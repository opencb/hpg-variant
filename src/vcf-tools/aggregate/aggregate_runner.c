/*
 * Copyright (c) 2012-2014 Cristina Yenyxe Gonzalez Garcia (EMBL-EBI)
 * Copyright (c) 2012-2014 Ignacio Medina (EMBL-EBI)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * hpg-variant is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-variant. If not, see <http://www.gnu.org/licenses/>.
 */

#include "aggregate_runner.h"


int run_aggregate(shared_options_data_t *shared_options_data, aggregate_options_data_t *options_data) {
    file_stats_t *file_stats = file_stats_new();
    
    // List that stores the batches of records filtered by each thread
    list_t *output_list[shared_options_data->num_threads];
    // List that stores which thread filtered the next batch to save
    list_t *next_token_list = malloc(sizeof(list_t));

    int ret_code;
    double start, stop, total;
    
    vcf_file_t *vcf_file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
    if (!vcf_file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }
    
    // Initialize variables related to the different threads
    for (int i = 0; i < shared_options_data->num_threads; i++) {
        output_list[i] = (list_t*) malloc(sizeof(list_t));
        list_init("input", 1, shared_options_data->num_threads * shared_options_data->batch_lines, output_list[i]);
    }
    list_init("next_token", shared_options_data->num_threads, INT_MAX, next_token_list);
    
    LOG_INFO("About to retrieve statistics from VCF file...\n");

#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            if (shared_options_data->batch_bytes > 0) {
                ret_code = vcf_parse_batches_in_bytes(shared_options_data->batch_bytes, vcf_file);
            } else if (shared_options_data->batch_lines > 0) {
                ret_code = vcf_parse_batches(shared_options_data->batch_lines, vcf_file);
            }

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            notify_end_parsing(vcf_file);
        }
        
#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            individual_t **individuals = NULL;
            khash_t(ids) *sample_ids = NULL;
            khash_t(str) *phenotype_ids = NULL;
            int num_phenotypes = 0;
            
            start = omp_get_wtime();
            
            int i = 0;
            vcf_batch_t *batch = NULL;
            while ((batch = fetch_vcf_batch(vcf_file)) != NULL) {
                LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->records->size, batch->records->capacity);

                // Divide the list of passed records in ranges of size defined in config file
                int num_chunks;
                int *chunk_sizes = NULL;
                array_list_t *input_records = batch->records;
                int *chunk_starts = create_chunks(input_records->size, 
                                                  ceil((float) shared_options_data->batch_lines / shared_options_data->num_threads), 
                                                  &num_chunks, &chunk_sizes);
                
                // OpenMP: Launch a thread for each range
                #pragma omp parallel for num_threads(shared_options_data->num_threads)
                for (int j = 0; j < num_chunks; j++) {
                    LOG_DEBUG_F("[%d] Stats invocation\n", omp_get_thread_num());
                    // Invoke variant stats
                    int index = omp_get_thread_num() % shared_options_data->num_threads;
                    ret_code = get_variants_stats((vcf_record_t**) (input_records->items + chunk_starts[j]),
                                                  chunk_sizes[j], individuals, sample_ids, num_phenotypes, output_list[index], file_stats); 
                }
                
                // Insert as many tokens as elements correspond to each thread,
                // linking the stats to some extra data from the variant to 
                // reconstruct the original line with extra INFO fields
                int variant_idx = 0;
                for (int t = 0; t < num_chunks; t++) {
                    for (int s = 0; s < chunk_sizes[t]; s++, variant_idx++) {
                        // Create auxiliary information
                        vcf_record_t *record = (vcf_record_t*) input_records->items[variant_idx];
                        variant_auxdata_t *aux = variant_auxdata_new(record);
                        // Link with variant stats
                        list_item_t *token_item = list_item_new(t, 0, aux);
                        list_insert_item(token_item, next_token_list);
                    }
                }
                
                free(chunk_starts);
                free(chunk_sizes);
                vcf_batch_free(batch);
                
                i++;
            }
            
            stop = omp_get_wtime();
            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
            
            // Decrease list writers count
            for (i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(next_token_list);
                list_decr_writers(output_list[i]);
            }
        }
        
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d writes the output\n", omp_get_thread_num());
            
            char *prefix_filename = calloc(strlen(shared_options_data->vcf_filename), sizeof(char));
            get_filename_from_path(shared_options_data->vcf_filename, prefix_filename);
            char default_filename[strlen(prefix_filename) + 7];
            sprintf(default_filename, "%s.aggregated", prefix_filename);
            
            char *aggregated_filename;
            FILE *output_fd = get_output_file(shared_options_data, default_filename, &aggregated_filename);
            LOG_INFO_F("Output filename = %s\n", aggregated_filename);
            
            // For each variant, generate a new line
            list_item_t *token_item = NULL;
            int header_written = 0;
            while ( token_item = list_remove_item(next_token_list) ) {
                // Write header and delimiter line
                if (!header_written) {
                    add_aggregator_header(vcf_file, options_data->overwrite);
                    write_vcf_header_nosamples(vcf_file, output_fd);
                    header_written = 1;
                    // TODO Write INFO headers for own "HPG_" annotations (if the "overwrite" flag is set)
                }
                
                variant_auxdata_t *aux_data = token_item->data_p; // Get the variant auxiliary data
                list_item_t *output_item = list_remove_item(output_list[token_item->id]); // Get the statistics
                assert(output_item);
                variant_stats_t *var_stats = output_item->data_p;
                
                // TODO Generate INFO field using statistics
                char *info = merge_info_and_stats(aux_data->info, var_stats, options_data->overwrite);
                
                // Combine with the rest of fields to get a vcf_record_t
                vcf_record_t *final_record = vcf_record_new();
                set_vcf_record_chromosome(var_stats->chromosome, strlen(var_stats->chromosome), final_record);
                set_vcf_record_position(var_stats->position, final_record);
                set_vcf_record_id(aux_data->id, strlen(aux_data->id), final_record);
                set_vcf_record_reference(var_stats->ref_allele, strlen(var_stats->ref_allele), final_record);
                set_vcf_record_alternate(var_stats->alt_alleles, strlen(var_stats->alt_alleles), final_record);
                set_vcf_record_quality(aux_data->quality, final_record);
                set_vcf_record_filter(aux_data->filter, strlen(aux_data->filter), final_record);
                set_vcf_record_info(info, strlen(info), final_record);
                
                // Write vcf_record_t to the file
                write_vcf_record(final_record, output_fd);
                
                // Free resources
                free(info);
                vcf_record_free(final_record);
                variant_stats_free(var_stats);
                variant_auxdata_free(aux_data);
                list_item_free(output_item);
                list_item_free(token_item);
            }

            fclose(output_fd);
            free(prefix_filename);
            free(aggregated_filename);
        }
    }
    
    free(file_stats);
    
    free(next_token_list);
    for (int i = 0; i < shared_options_data->num_threads; i++) {
        free(output_list[i]);
    }
    
    vcf_close(vcf_file);
    
    return 0;
}

void add_aggregator_header(vcf_file_t *vcf_file, int overwrite) {
    // TODO This should use the --config argument
    array_list_t *config_search_paths = get_configuration_search_paths(0, NULL);
    char *configuration_file = retrieve_config_file("vcf-info-fields.conf", config_search_paths);
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, configuration_file);
    assert(ret_code);
    
    // Add headers for fields created by HPG Variant
    if (!overwrite) {
        assert(add_aggregator_header_entry(config, "HPG_AC", vcf_file));
        assert(add_aggregator_header_entry(config, "HPG_AF", vcf_file));
        assert(add_aggregator_header_entry(config, "HPG_AN", vcf_file));
        assert(add_aggregator_header_entry(config, "HPG_GTC", vcf_file));
    } else {
        assert(add_aggregator_header_entry(config, "AC", vcf_file));
        assert(add_aggregator_header_entry(config, "AF", vcf_file));
        assert(add_aggregator_header_entry(config, "AN", vcf_file));
        assert(add_aggregator_header_entry(config, "HPG_GTC", vcf_file));
    }
}

int add_aggregator_header_entry(config_t *config, char *info_field_name, vcf_file_t *vcf_file) {
    const char *field_value;
    int ret_code = config_lookup_string(config, info_field_name, &field_value);
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("Information about subfield '%s' of INFO not found in configuration file\n", info_field_name);
        return 0;
    } else {
        vcf_header_entry_t* header_entry = vcf_header_entry_new();
        set_vcf_header_entry_name("INFO", 4, header_entry);
        add_vcf_header_entry_value(field_value, strlen(field_value), header_entry);
        return add_vcf_header_entry(header_entry, vcf_file);
    }
}

char *merge_info_and_stats(char *info, variant_stats_t *stats, int overwrite) {
    khash_t(info_fields) *fields_hash = kh_init(info_fields);
    
    // If the original INFO column was not empty, initialize hash map with its fields
    if (strcmp(".", info)) {
        int num_fields;
        char *dup_info = strdup(info);
        char **fields = split(dup_info, ";", &num_fields);

        for (int i = 0; i < num_fields; i++) {
            int num_subfields;
            char **subfields = split(fields[i], "=", &num_subfields);

            int ret = add_to_hash(fields_hash, subfields[0], (num_subfields > 1) ? subfields[1] : NULL);
            assert(ret);

            free(fields[i]);
            free(subfields);
        }
        free(fields);
        free(dup_info);
    }
    
    int ret;
    // Add some fields from statistics to hash map (AC, AF, AN, MAF)
    // AC: ALT alleles counts
    int *new_AC = stats->alleles_count + 1;
    char *AC_str = (char*) calloc (6 * stats->num_alleles, sizeof(char));
    for (int i = 0; i < stats->num_alleles - 1; i++) {
        sprintf(AC_str + strlen(AC_str), "%d,", new_AC[i]);
    }
    AC_str[strlen(AC_str)-1] = 0;
    ret = add_to_hash(fields_hash, "HPG_AC", AC_str); assert(ret);
    
    // AF: ALT alleles frequencies
    float *new_AF = stats->alleles_freq + 1;
    char *AF_str = (char*) calloc (6 * stats->num_alleles, sizeof(char));
    for (int i = 0; i < stats->num_alleles - 1; i++) {
        sprintf(AF_str + strlen(AF_str), "%.3f,", new_AF[i]);
    }
    AF_str[strlen(AF_str)-1] = 0;
    ret = add_to_hash(fields_hash, "HPG_AF", AF_str); assert(ret);
    
    // AN: Number of all alleles together
    int new_AN = 0; 
    for (int i = 0; i < stats->num_alleles; i++) {
        new_AN += stats->alleles_count[i];
    }
    char *AN_str = (char*) malloc (16 * sizeof(char));
    sprintf(AN_str, "%d", new_AN);
    ret = add_to_hash(fields_hash, "HPG_AN", AN_str); assert(ret);
    
//    // Minor allele frequency
//    float *maf = stats->maf;
//    char *maf_str = (char*) malloc (5 * sizeof(char)); // MAF format = x.abc
//    sprintf(maf_str, "%.3f", maf);
//    ret = add_to_hash(fields_hash, "HPG_MAF", maf_str); assert(ret);
    
    // GTC: Genotype counts
    char *GTC_str = report_variant_genotypes_stats(stats);
    ret = add_to_hash(fields_hash, "HPG_GTC", GTC_str); assert(ret);
    
    // Merge them considering the overwriting flag
    char *new_info = (char*) calloc (strlen(info) + 128, sizeof(char));
    for (int k = kh_begin(fields_hash); k < kh_end(fields_hash); k++) {
        if (kh_exist(fields_hash, k)) {
            char *key = kh_key(fields_hash, k);
            char *value = kh_value(fields_hash, k);
            assert(value);
            
            if (overwrite) {
                if (!strcmp(key, "AC") || !strcmp(key, "AF") || !strcmp(key, "AN")) {
                    // Ignore, they are going to be overwritten anyway
                } else if (!strcmp(key, "HPG_AC") || !strcmp(key, "HPG_AF") || !strcmp(key, "HPG_AN")) {
                    // Remove the HPG part
                    char *suffix = key + 4;
                    strcat(new_info, suffix);
                    strcat(new_info, "=");
                    strcat(new_info, value);
                    strcat(new_info, ";");
                } else {
                    // Non-overwritten field from the original INFO column
                    // or fields calculated by HPG Variant
                    strcat(new_info, key);
                    strcat(new_info, "=");
                    strcat(new_info, value);
                    strcat(new_info, ";");
                }
            } else {
                strcat(new_info, key);
                strcat(new_info, "=");
                strcat(new_info, value);
                strcat(new_info, ";");
            }
            
            if (!starts_with(key, "HPG_")) { free(key); }
            free(value);
        }
    }
    new_info[strlen(new_info)-1] = 0;
    
    kh_destroy(info_fields, fields_hash);
    return new_info;
}

int add_to_hash(kh_info_fields_t *hash, char *key, char *value) {
    int ret;
    khiter_t k = kh_put(info_fields, hash, key, &ret);
    if (!ret) kh_del(info_fields, hash, k);
    kh_value(hash, k) = value;
    
    return ret;
}

char* report_variant_genotypes_stats(variant_stats_t *var_stats) {
    // The field representation will be allele1/allele2:count, like 0/1:204
    // It also includes the count of missing genotypes, hence the +1
    char *genotype_count_field = (char*) calloc ((var_stats->num_alleles * var_stats->num_alleles + 1) * 16, sizeof(char));
    int num_chars_written = 0;
    int gt_count = 0;
//    float gt_freq = 0;
    
    for (int i = 0; i < var_stats->num_alleles; i++) {
        for (int j = i; j < var_stats->num_alleles; j++) {
            int idx1 = i * var_stats->num_alleles + j;
            if (i == j) {
                gt_count = var_stats->genotypes_count[idx1];
//                gt_freq = var_stats->genotypes_freq[idx1];
            } else {
                int idx2 = j * var_stats->num_alleles + i;
                gt_count = var_stats->genotypes_count[idx1] + var_stats->genotypes_count[idx2];
//                gt_freq = var_stats->genotypes_freq[idx1] + var_stats->genotypes_freq[idx2];
            }

            num_chars_written += sprintf(genotype_count_field + num_chars_written, "%d/%d:%d,", i, j, gt_count);
        }
    }
    
    num_chars_written += sprintf(genotype_count_field + num_chars_written, "./.:%d,", var_stats->missing_genotypes);
    
    genotype_count_field[num_chars_written] = 0;
    
    return genotype_count_field;
}

variant_auxdata_t* variant_auxdata_new(vcf_record_t *record) {
    variant_auxdata_t *aux = (variant_auxdata_t*) malloc (sizeof(variant_auxdata_t));
    aux->id = strndup(record->id, record->id_len);
    aux->quality = record->quality;
    aux->filter = strndup(record->filter, record->filter_len);
    aux->info = strndup(record->info, record->info_len);
    return aux;
}

void variant_auxdata_free(variant_auxdata_t *data) {
    free(data->filter);
    free(data->id);
    free(data->info);
    free(data);
}
