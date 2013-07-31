/*
 * Copyright (c) 2012-2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
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

#include "stats.h"


int run_stats(shared_options_data_t *shared_options_data, stats_options_data_t *options_data) {
    file_stats_t *file_stats = file_stats_new();
    sample_stats_t **sample_stats;
    
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
    
    ped_file_t *ped_file = NULL;
    if (shared_options_data->ped_filename) {
        ped_file = ped_open(shared_options_data->ped_filename);
        if (!ped_file) {
            LOG_FATAL("PED file does not exist!\n");
        }
        if(options_data->variable) {
            set_variable_field(options_data->variable, 0, ped_file);
        } else {
            set_variable_field("PHENO", 6, ped_file);
        }
        
        if(options_data->variable_groups) {
            int n, m;
            char *variable_groups = strdup(options_data->variable_groups);
            char **groups;
            char **phenos_in_group;
            groups = split(variable_groups, ":", &n);
            for(int i = 0; i < n; i++){
                phenos_in_group = split(groups[i], ",", &m);
                if(set_phenotype_group(phenos_in_group, m, ped_file) < 0) {
                    LOG_ERROR("Variable can't appear in two groups\n");
                    return DUPLICATED_VARIABLE;
                }
                free(phenos_in_group);
            }
            ped_file->accept_new_values = 0;
            
            free(variable_groups);
            free(groups);
        } else {
            ped_file->accept_new_values = 1;
        }
        if(options_data->phenotype) {
            int n;
            char* phenotypes = strdup(options_data->phenotype);
            char** pheno_values = split(phenotypes, ",", &n);
            if(n != 2) {
                LOG_ERROR("To handle case-control test, only two phenotypes are supported\n");
                return MORE_THAN_TWO_PHENOTYPES;
            } else {
                set_unaffected_phenotype(pheno_values[0],ped_file);
                set_affected_phenotype(pheno_values[1],ped_file);
            }
        } else {
            set_unaffected_phenotype("1", ped_file);
            set_affected_phenotype("2", ped_file);
        }
        
        LOG_INFO("About to read PED file...\n");
        // Read PED file before doing any processing
        ret_code = ped_read(ped_file);
        if (ret_code != 0) {
            LOG_FATAL_F("Can't read PED file: %s\n", ped_file->filename);
        }
        if(!ped_file->num_field) {
            LOG_ERROR_F("Can't find the specified field \"%s\" in file: %s \n", options_data->variable, ped_file->filename);
            return VARIABLE_FIELD_NOT_FOUND;
        }
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
            int num_phenotypes;
            
            start = omp_get_wtime();
            
            int i = 0;
            vcf_batch_t *batch = NULL;
            while ((batch = fetch_vcf_batch(vcf_file)) != NULL) {
                if (i == 0) {
                    sample_stats = malloc (get_num_vcf_samples(vcf_file) * sizeof(sample_stats_t*));
                    for (int j = 0; j < get_num_vcf_samples(vcf_file); j++) {
                        sample_stats[j] = sample_stats_new(array_list_get(j, vcf_file->samples_names));
                    }
                    
                    if (ped_file) {
                        // Create map to associate the position of individuals in the list of samples defined in the VCF file
                        sample_ids = associate_samples_and_positions(vcf_file);
                        // Sort individuals in PED as defined in the VCF file
                        individuals = sort_individuals(vcf_file, ped_file);
                        // Get the khash of the phenotypes in PED file
                        phenotype_ids = get_phenotypes(ped_file);
                        num_phenotypes = get_num_variables(ped_file);
                    }
                }
                
                if (i % 50 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->records->size, batch->records->capacity);
                }

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
                    // Invoke variant stats and/or sample stats when applies
                    if (options_data->variant_stats) {
                        int index = omp_get_thread_num() % shared_options_data->num_threads;
                        ret_code = get_variants_stats((vcf_record_t**) (input_records->items + chunk_starts[j]),
                                                      chunk_sizes[j], individuals, sample_ids,num_phenotypes, output_list[index], file_stats); 
                    }
                    
                    if (options_data->sample_stats) {
                        ret_code |= get_sample_stats((vcf_record_t**) (input_records->items + chunk_starts[j]), 
                                                      chunk_sizes[j], individuals, sample_ids, sample_stats, file_stats);
                    }
                }
                
                if (options_data->variant_stats) {
                    // Insert as many tokens as elements correspond to each thread
                    for (int t = 0; t < num_chunks; t++) {
                        for (int s = 0; s < chunk_sizes[t]; s++) {
                            list_item_t *token_item = list_item_new(t, 0, NULL);
                            list_insert_item(token_item, next_token_list);
                        }
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
            
            if (sample_ids) { kh_destroy(ids, sample_ids); }
            if (individuals) { free(individuals); }
        }
        
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d writes the output\n", omp_get_thread_num());
            
            char *stats_prefix = get_vcf_stats_filename_prefix(shared_options_data->vcf_filename, 
                                                               shared_options_data->output_filename, 
                                                               shared_options_data->output_directory);
            
            // File names and descriptors for output to plain text files
            char *stats_filename, *summary_filename, *phenotype_filename;
            FILE *stats_fd, *summary_fd, **phenotype_fd;
            
            char *stats_db_name;
            sqlite3 *db = NULL;
            khash_t(stats_chunks) *hash;
            
            khash_t(str) *phenotype_ids;
            int num_phenotypes;
            if(ped_file){
                phenotype_ids = get_phenotypes(ped_file);
                num_phenotypes = get_num_variables(ped_file);
            }
            
            if (options_data->save_db) {
                delete_files_by_extension(shared_options_data->output_directory, "db");
                stats_db_name = calloc(strlen(stats_prefix) + strlen(".db") + 2, sizeof(char));
                sprintf(stats_db_name, "%s.db", stats_prefix);
                create_stats_db(stats_db_name, VCF_CHUNKSIZE, create_vcf_query_fields, &db);
                hash = kh_init(stats_chunks);
            }
            
            // Write variant (and global) statistics
            if (options_data->variant_stats) {
                stats_filename = get_variant_stats_output_filename(stats_prefix);
                if (!(stats_fd = fopen(stats_filename, "w"))) {
                    LOG_FATAL_F("Can't open file for writing statistics of variants: %s\n", stats_filename);
                }
                
                //Open one file for each phenotype
                if(ped_file){
                    phenotype_fd = malloc(sizeof(FILE*)*num_phenotypes);
                    if(options_data->variable_groups){
                        int n;
                        char *variable_groups = strdup(options_data->variable_groups);
                        char ** names = split(variable_groups, ":", &n);
                        for(int i = 0; i < n; i++) {
                            phenotype_filename = get_variant_phenotype_stats_output_filename(stats_prefix, names[i]);
                            if(!(phenotype_fd[i] = fopen(phenotype_filename, "w"))) {
                                LOG_FATAL_F("Can't open file for writing statistics of variants per phenotype: %s\n", stats_filename);
                            }
                            free(phenotype_filename);
                        }
                        free(names);
                        free(variable_groups);
                    } else {
                 
                        for (khint_t i = kh_begin(phenotype_ids); i != kh_end(phenotype_ids); ++i) {
                            if (!kh_exist(phenotype_ids,i)) continue;
                            
                            phenotype_filename = get_variant_phenotype_stats_output_filename(stats_prefix, kh_key(phenotype_ids,i));
                            if(!(phenotype_fd[kh_val(phenotype_ids,i)] = fopen(phenotype_filename, "w"))) {
                                LOG_FATAL_F("Can't open file for writing statistics of variants per phenotype: %s\n", stats_filename);
                            }
                            free(phenotype_filename);
                        }
                    }
                }
                // Write header
                report_vcf_variant_stats_header(stats_fd);
                if(ped_file){
                    for(int i = 0; i < num_phenotypes; i++)
                        report_vcf_variant_phenotype_stats_header(phenotype_fd[i]);
                }
                
                // For each variant, generate a new line
                int avail_stats = 0;
                variant_stats_t *var_stats_batch[VCF_CHUNKSIZE];
                list_item_t *token_item = NULL, *output_item = NULL;
                while ( token_item = list_remove_item(next_token_list) ) {
                    output_item = list_remove_item(output_list[token_item->id]);
                    assert(output_item);
                    var_stats_batch[avail_stats] = output_item->data_p;
                    avail_stats++;
                    
                    // Run only when certain amount of stats is available
                    if (avail_stats >= VCF_CHUNKSIZE) {
                        report_vcf_variant_stats(stats_fd, db, hash, avail_stats, var_stats_batch);
                        
                        if(ped_file)
                            for(int i = 0; i < num_phenotypes; i++)
                                report_vcf_variant_phenotype_stats(phenotype_fd[i], avail_stats, var_stats_batch, i);

                        // Free all stats from the "batch"
                        for (int i = 0; i < avail_stats; i++) {
                            variant_stats_free(var_stats_batch[i]);
                        }
                        avail_stats = 0;
                    }
                    
                    // Free resources
                    list_item_free(output_item);
                    list_item_free(token_item);
                }
                
                if (avail_stats > 0) {
                    report_vcf_variant_stats(stats_fd, db, hash, avail_stats, var_stats_batch);
                    
                    if(ped_file)
                        for(int i = 0; i < num_phenotypes; i++)
                            report_vcf_variant_phenotype_stats(phenotype_fd[i], avail_stats, var_stats_batch, i);

                    // Free all stats from the "batch"
                    for (int i = 0; i < avail_stats; i++) {
                        variant_stats_free(var_stats_batch[i]);
                    }
                    avail_stats = 0;
                }
                
                // Write whole file stats (data only got when launching variant stats)
                summary_filename = get_vcf_file_stats_output_filename(stats_prefix);
                if (!(summary_fd = fopen(summary_filename, "w"))) {
                    LOG_FATAL_F("Can't open file for writing statistics summary: %s\n", summary_filename);
                }
                report_vcf_summary_stats(summary_fd, db, file_stats);
                
                free(stats_filename);
                free(summary_filename);
                
                // Close variant stats file
                if (stats_fd) { fclose(stats_fd); }
                if (summary_fd) { fclose(summary_fd); }
				if(ped_file){
		            for(int i = 0; i < num_phenotypes; i++)
		                if(phenotype_fd[i]) fclose(phenotype_fd[i]);
					free(phenotype_fd);
				}
            }
            
            // Write sample statistics
            if (options_data->sample_stats) {
                stats_filename = get_sample_stats_output_filename(stats_prefix);
                if (!(stats_fd = fopen(stats_filename, "w"))) {
                    LOG_FATAL_F("Can't open file for writing statistics of samples: %s\n", stats_filename);
                }
                
                report_vcf_sample_stats_header(stats_fd);
                report_vcf_sample_stats(stats_fd, NULL, vcf_file->samples_names->size, sample_stats);
                
                // Close sample stats file
                free(stats_filename);
                if (stats_fd) { fclose(stats_fd); }
            }
            
            free(stats_prefix);
            
            if (db) {
                insert_chunk_hash(VCF_CHUNKSIZE, hash, db);
                create_stats_index(create_vcf_index, db);
                close_stats_db(db, hash);
            }
            
        }
    }
    
    for (int i = 0; i < get_num_vcf_samples(vcf_file); i++) {
        sample_stats_free(sample_stats[i]);
    }
    free(sample_stats);
    free(file_stats);
    
    free(next_token_list);
    for (int i = 0; i < shared_options_data->num_threads; i++) {
        free(output_list[i]);
    }
    
    vcf_close(vcf_file);
    if (ped_file) { ped_close(ped_file, 1,1); }
    
    return 0;
}
