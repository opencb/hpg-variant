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

#include "merge.h"

static int *get_format_indices_per_file(vcf_record_file_link **position_in_files, int position_occurrences, 
                                        vcf_file_t **files, int num_files, array_list_t *format_fields);
static inline void concat_sample_allele_tolerant(int allele, int len, char *sample, cp_hashtable* alleles_table, char *reference, char **split_alternates);
static inline void concat_sample_allele_strict(int allele, int len, char *sample, cp_hashtable* alleles_table, char **split_alternates);
static char *get_empty_sample(int num_format_fields, int gt_pos, merge_options_data_t *options);


int merge_vcf_headers(vcf_file_t** files, int num_files, merge_options_data_t* options, list_t* output_list) {
    cp_hashtable *filter_entries = cp_hashtable_create_by_option(COLLECTION_MODE_PLAIN, 16, 
                                                                 cp_hash_istring, (cp_compare_fn) strcmp, 
                                                                 NULL, free,
                                                                 NULL, free);
    cp_hashtable *format_entries = cp_hashtable_create_by_option(COLLECTION_MODE_PLAIN, 16, 
                                                                 cp_hash_istring, (cp_compare_fn) strcmp, 
                                                                 NULL, free,
                                                                 NULL, free);
    cp_hashtable *info_entries = cp_hashtable_create_by_option(COLLECTION_MODE_PLAIN, 16, 
                                                                 cp_hash_istring, (cp_compare_fn) strcmp, 
                                                                 NULL, free,
                                                                 NULL, free);
    cp_hashtable *other_entries = cp_hashtable_create_by_option(COLLECTION_MODE_PLAIN, 16, 
                                                                 cp_hash_istring, (cp_compare_fn) strcmp, 
                                                                 NULL, free,
                                                                 NULL, free);
    
    vcf_file_t *file;
    vcf_header_entry_t *in_entry, *out_entry, *aux;
    
    for (int i = 0; i < num_files; i++) {
        file = files[i];
        for (int j = 0; j < file->header_entries->size; j++) {
            out_entry = NULL;
            in_entry = array_list_get(j, file->header_entries);
//             printf("1) value = %s\n", array_list_get(0, in_entry->values));
            // Merge FILTER, FORMAT and other entries that are __not__ INFO ones
            if (!strncmp("FILTER", in_entry->name, in_entry->name_len)) {
                aux = cp_hashtable_get(filter_entries, array_list_get(0, in_entry->values));
                if (!aux) {
                    cp_hashtable_put(filter_entries, array_list_get(0, in_entry->values), in_entry);
                    out_entry = in_entry;
                }
            } else if (!strncmp("FORMAT", in_entry->name, in_entry->name_len)) {
                aux = cp_hashtable_get(format_entries, array_list_get(0, in_entry->values));
                if (!aux) {
                    cp_hashtable_put(format_entries, array_list_get(0, in_entry->values), in_entry);
                    out_entry = in_entry;
                }
            } else if (!strncmp("INFO", in_entry->name, in_entry->name_len)) {
                aux = cp_hashtable_get(info_entries, array_list_get(0, in_entry->values));
                if (!aux) {
                    cp_hashtable_put(info_entries, array_list_get(0, in_entry->values), in_entry);
                    out_entry = in_entry;
                }
            } else {
                aux = cp_hashtable_get(other_entries, array_list_get(0, in_entry->values));
                if (!aux) {
                    cp_hashtable_put(other_entries, array_list_get(0, in_entry->values), in_entry);
                    out_entry = in_entry;
                }
            }
            
            // Insert into the list of entries of the output file
            if (out_entry) {
                list_item_t *item = list_item_new(i, MERGED_HEADER, out_entry);
                list_insert_item(item, output_list);
            }
        }
    }
    
    // Create INFO entries
    char *info_field;
    char *field_value;
    
    char *configuration_file = retrieve_config_file("vcf-info-fields.conf", options->config_search_paths);
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, configuration_file);
    assert(ret_code);

    for (int i = 0; i < options->num_info_fields; i++) {
        info_field = options->info_fields[i];
        // TODO search info field in config file
        ret_code = config_lookup_string(config, info_field, &field_value);
        if (ret_code == CONFIG_FALSE) {
            LOG_ERROR_F("Information about subfield '%s' of INFO not found in configuration file\n", info_field);
        } else {
            aux = cp_hashtable_get(info_entries, field_value);
            if (!aux) {
//                 printf("2) value = %s\n", field_value);
                out_entry = vcf_header_entry_new();
                set_vcf_header_entry_name("INFO", 4, out_entry);
                add_vcf_header_entry_value(field_value, strlen(field_value), out_entry);
                
                list_item_t *item = list_item_new(i, MERGED_HEADER, out_entry);
                list_insert_item(item, output_list);
            }
        }
    }
    
    config_destroy(config);
    free(config);
    
    cp_hashtable_destroy(filter_entries);
    cp_hashtable_destroy(format_entries);
    cp_hashtable_destroy(info_entries);
    cp_hashtable_destroy(other_entries);
}


array_list_t *merge_vcf_sample_names(vcf_file_t **files, int num_files) {
    khash_t(names) *samples = kh_init(names);
    // Use a khash for checking repeated samples
    for (int i = 0; i < num_files; i++) {
        for (int j = 0; j < files[i]->samples_names->size; j++) {
            char *sample_name = array_list_get(j, files[i]->samples_names);
            khiter_t iter = kh_get(names, samples, sample_name);
            if (iter != kh_end(samples)) {
                kh_destroy(names, samples);
                LOG_ERROR_F("Sample %s appears more than once.\n", sample_name);
                return NULL;
            } else {
                int ret;
                iter = kh_put(names, samples, strdup(sample_name), &ret);
            }
        }
    }
    kh_destroy(names, samples);
    
    // If no errors were found, return the list of all sample names
    array_list_t *sample_names = array_list_new(get_num_vcf_samples(files[0]) * 2, 2, COLLECTION_MODE_ASYNCHRONIZED);
    for (int i = 0; i < num_files; i++) {
        array_list_insert_all(files[i]->samples_names->items, files[i]->samples_names->size, sample_names);
    }
    return sample_names;
}


int merge_vcf_records(array_list_t **records_by_position, int num_positions, vcf_file_t **files, int num_files, merge_options_data_t *options, list_t *output_list) {
    int info;
    vcf_record_t *merged;
    for (int i = 0; i < num_positions; i++) {
        merged = merge_position((vcf_record_file_link **) records_by_position[i]->items, records_by_position[i]->size, files, num_files, options, &info);
        
        if (merged) {
            list_item_t *item = list_item_new(i, MERGED_RECORD, merged);
            list_insert_item(item, output_list);
        }
    }
    
    printf("Merging %d variants...\n", num_positions);
    
    return 0;
}

vcf_record_t *merge_position(vcf_record_file_link **position_in_files, int position_occurrences, 
                             vcf_file_t **files, int num_files, merge_options_data_t *options, int *err_code) {
    vcf_record_t *result = vcf_record_new();
    
    // Check consistency among chromosome, position and reference
    for (int i = 0; i < position_occurrences-1; i++) {
        assert(position_in_files[i]);
        assert(position_in_files[i+1]);
        assert(position_in_files[i]->record->chromosome);
        assert(position_in_files[i+1]->record->chromosome);
        assert(position_in_files[i]->record->chromosome_len);
        
        if (strncmp(position_in_files[i]->record->chromosome, position_in_files[i+1]->record->chromosome, position_in_files[i]->record->chromosome_len)) {
            LOG_ERROR_F("Positions %.*s:%ld and %.*s:%ld can't be merged: Discordant chromosome\n", 
                        position_in_files[i]->record->chromosome_len, position_in_files[i]->record->chromosome, position_in_files[i]->record->position,
                        position_in_files[i+1]->record->chromosome_len, position_in_files[i+1]->record->chromosome, position_in_files[i+1]->record->position);
            *err_code = DISCORDANT_CHROMOSOME;
            return NULL;
        } else if (position_in_files[i]->record->position != position_in_files[i+1]->record->position) {
            LOG_ERROR_F("Positions %.*s:%ld and %.*s:%ld can't be merged: Discordant position\n", 
                        position_in_files[i]->record->chromosome_len, position_in_files[i]->record->chromosome, position_in_files[i]->record->position,
                        position_in_files[i+1]->record->chromosome_len, position_in_files[i+1]->record->chromosome, position_in_files[i+1]->record->position);
            *err_code = DISCORDANT_POSITION;
            return NULL;
        } else if (strncmp(position_in_files[i]->record->reference, position_in_files[i+1]->record->reference, position_in_files[i]->record->reference_len)) {
            // If the strict reference flag is set, both reference alleles must be the same
            // Otherwise, accept them if their first nucleotide is the same
            if (options->strict_reference || position_in_files[i]->record->reference[0] != position_in_files[i+1]->record->reference[0]) {
                LOG_ERROR_F("Position %.*s:%ld can't be merged: Discordant reference alleles (%.*s, %.*s)\n", 
                            position_in_files[i]->record->chromosome_len, position_in_files[i]->record->chromosome, position_in_files[i]->record->position,
                            position_in_files[i]->record->reference_len, position_in_files[i]->record->reference, 
                            position_in_files[i+1]->record->reference_len, position_in_files[i+1]->record->reference);
                *err_code = DISCORDANT_REFERENCE;
                return NULL;
            }
        }
    }
    
    // Once checked, copy constant fields
    set_vcf_record_chromosome(strndup(position_in_files[0]->record->chromosome, position_in_files[0]->record->chromosome_len), 
                              position_in_files[0]->record->chromosome_len, result);
    set_vcf_record_position(position_in_files[0]->record->position, result);
    set_vcf_record_reference(strndup(position_in_files[0]->record->reference, position_in_files[0]->record->reference_len),
                             position_in_files[0]->record->reference_len, result);
    
    // Get first non-dot ID
    // TODO what can we do when having several ID?
    char *id = merge_id_field(position_in_files, position_occurrences);
    set_vcf_record_id(id, strlen(id), result);
    
    // Calculate weighted mean of the quality
    set_vcf_record_quality(merge_quality_field(position_in_files, position_occurrences), result);
    
    // Concatenate alternates and set their order (used later to assign samples' alleles number)
    cp_hashtable *alleles_table = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP, 8, 
                                                                cp_hash_istring, 
                                                                (cp_compare_fn) strcasecmp,
                                                                NULL, free,
                                                                NULL, free);
    char *alternate = merge_alternate_field(position_in_files, position_occurrences, alleles_table);
    set_vcf_record_alternate(alternate, strlen(alternate), result);
    
    // Concatenate failed filters
    char *filter = merge_filter_field(position_in_files, position_occurrences);
    set_vcf_record_filter(filter, strlen(filter), result);
    
    // Get the union of all FORMAT fields with their corresponding position
    // Include INFO and FILTER fields in FORMAT when required by the user
    array_list_t *format_fields = array_list_new(16, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    format_fields->compare_fn = strcasecmp;
    char *format = merge_format_field(position_in_files, position_occurrences, options, format_fields);
    set_vcf_record_format(format, strlen(format), result);
    
    // Get indexes for the format in each file
    // TODO Illustrate with an example
    int *format_indices = get_format_indices_per_file(position_in_files, position_occurrences, files, num_files, format_fields);
    
    // Create the text for empty samples
    int filter_pos, info_pos, gt_pos;
    char *dupaux = strndup(result->format, result->format_len);
    filter_pos = get_field_position_in_format("SFT", dupaux);
    free(dupaux);
    
    dupaux = strndup(result->format, result->format_len);
    info_pos = get_field_position_in_format("IN", dupaux);
    free(dupaux);
    
    dupaux = strndup(result->format, result->format_len);
    gt_pos = get_field_position_in_format("GT", dupaux);
    free(dupaux);
    
    char *empty_sample = get_empty_sample(format_fields->size, gt_pos, options);
//     printf("empty sample = %s\n", empty_sample);

    // Generate samples using the reordered FORMAT fields and the new alleles' numerical values
    // Include INFO and FILTER fields from the original file when required by the user
    array_list_free(result->samples, NULL);
    result->samples = merge_samples(position_in_files, position_occurrences, files, num_files, alleles_table, 
                                    format_fields, format_indices, empty_sample, gt_pos, filter_pos, info_pos, options);

    // Merge INFO field
    char *info = merge_info_field(position_in_files, position_occurrences, options->info_fields, options->num_info_fields,
                                  result, alleles_table, empty_sample);
    set_vcf_record_info(info, strlen(info), result);

    free(empty_sample);
    free(format_indices);
    array_list_free(format_fields, free);
    cp_hashtable_destroy(alleles_table);

    return result;
}


/* ******************************
 *    Field merging functions   *
 * ******************************/

char* merge_id_field(vcf_record_file_link** position_in_files, int position_occurrences) {
    vcf_record_t *input;
    char *result = NULL;
    
    for (int i = 0; i < position_occurrences; i++) {
        input = position_in_files[i]->record;
        if (strncmp(".", input->id, input->id_len)) {
            result = strndup(input->id, input->id_len);
            break;
        }
    }
    
    if (result == NULL) {
        result = strdup(".");
    }
    
    return result;
}


char* merge_alternate_field(vcf_record_file_link** position_in_files, int position_occurrences, cp_hashtable* alleles_table) {
    vcf_file_t *file;
    vcf_record_t *input;
    size_t max_len = 0, concat_len = 0;
    char *reference = strndup(position_in_files[0]->record->reference, position_in_files[0]->record->reference_len);
    char *alternates = NULL, *cur_alternate = NULL, *aux = NULL;
    char *pos_alternates = NULL, **split_pos_alternates = NULL;
    int cur_alternate_len, num_pos_alternates;
    
    int cur_index = 0;
    int *allele_index = (int*) calloc (1, sizeof(int)); *allele_index = cur_index;
    cp_hashtable_put(alleles_table, reference, allele_index);
    
    for (int i = 0; i < position_occurrences; i++) {
        file = position_in_files[i]->file;
        input = position_in_files[i]->record;

        // Check reference allele
        cur_alternate = strndup(input->reference, input->reference_len);
        cur_alternate_len = input->reference_len;
        if (strcmp(reference, cur_alternate) && !cp_hashtable_contains(alleles_table, cur_alternate)) {
            if (!alternates) {
                alternates = strdup(cur_alternate);
                max_len = cur_alternate_len;
            } else {
                char *aux = realloc(alternates, max_len + cur_alternate_len + 2);
                if (aux) {
                    strncat(aux, ",", 1);
        //             printf("1) cur_alternate = %.*s\naux = %s\n---------\n", concat_len, cur_alternate, aux);
                    strncat(aux, cur_alternate, cur_alternate_len);
                    alternates = aux;
        //             printf("2) alternates = %s\n---------\n", alternates);
                } else {
                    LOG_FATAL_F("Can't allocate memory for alternate alleles in position %.*s:%ld\n",
                                input->chromosome_len, input->chromosome, input->position);
                }

                max_len += cur_alternate_len + 1;
            }

            // In case the allele is not in the hashtable, insert it with a new index
            allele_index = (int*) calloc (1, sizeof(int));
            *allele_index = ++cur_index;
            cp_hashtable_put(alleles_table, cur_alternate, allele_index);
        } else {
            free(cur_alternate);
        }
        
        // Check alternate alleles
        pos_alternates = strndup(input->alternate, input->alternate_len);
        split_pos_alternates = split(pos_alternates, ",", &num_pos_alternates);

        for (int j = 0; j < num_pos_alternates; j++) {
            cur_alternate = split_pos_alternates[j];
            cur_alternate_len = strlen(cur_alternate);

            if (!cp_hashtable_contains(alleles_table, cur_alternate)) {
                if (!alternates) {
                    alternates = strdup(cur_alternate);
                    max_len = cur_alternate_len;
                } else {
                    char *aux = realloc(alternates, max_len + cur_alternate_len + 2);
                    if (aux) {
                        strncat(aux, ",", 1);
            //             printf("1) cur_alternate = %.*s\naux = %s\n---------\n", concat_len, cur_alternate, aux);
                        strncat(aux, cur_alternate, cur_alternate_len);
                        alternates = aux;
            //             printf("2) alternates = %s\n---------\n", alternates);
                    } else {
                        LOG_FATAL_F("Can't allocate memory for alternate alleles in position %.*s:%ld\n",
                                    input->chromosome_len, input->chromosome, input->position);
                    }

                    max_len += cur_alternate_len + 1;
                }
                
                // In case the allele is not in the hashtable, insert it with a new index
                allele_index = (int*) calloc (1, sizeof(int));
                *allele_index = ++cur_index;
                cp_hashtable_put(alleles_table, cur_alternate, allele_index);
            } else {
                free(cur_alternate);
            }
        }

        free(pos_alternates);
        free(split_pos_alternates);
    }

    return alternates;
}


float merge_quality_field(vcf_record_file_link** position_in_files, int position_occurrences) {
    vcf_file_t *file;
    vcf_record_t *input;
    float accum_quality = 0.0f;
    int total_samples = 0;
    int file_num_samples;
    float result;
    
    for (int i = 0; i < position_occurrences; i++) {
        file = position_in_files[i]->file;
        input = position_in_files[i]->record;
        file_num_samples = get_num_vcf_samples(file);
        
        if (input->quality > 0) {
            accum_quality += input->quality * file_num_samples;
        }
        total_samples += file_num_samples;
    }
    
    if (total_samples > 0) {
        result = accum_quality / total_samples;
    } else {
        result = -1;
    }
    
    return result;
}


char* merge_filter_field(vcf_record_file_link** position_in_files, int position_occurrences) {
    char *result;   
    vcf_record_t *input;
    int filter_text_len = 0;
    array_list_t *failed_filters = array_list_new(8, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    failed_filters->compare_fn = strcasecmp;
    
    // Flags
    int pass_found = 0;
    int miss_found = 0;
    
    // Register all failed filters
    for (int i = 0; i < position_occurrences; i++) {
        input = position_in_files[i]->record;
        
        if (!strncmp(input->filter, "PASS", input->filter_len)) {
            pass_found = 1;
        } else if (!strncmp(input->filter, ".", input->filter_len)) {
            miss_found = 1;
        } else {
            char *filter_field = strndup(input->filter, input->filter_len);
            int num_filters;
            char **all_filters = split(filter_field, ";", &num_filters);
            
            for (int j = 0; j < num_filters; j++) {
                char *filter = all_filters[j];
                if (!array_list_contains(filter, failed_filters)) {
                    array_list_insert(filter, failed_filters);
                    filter_text_len += strlen(filter) + 1; // concat field + ","
                } else {
                    free(filter);
                }
            }
            free(all_filters);
            free(filter_field);
        }
    }
    
    // Should there be any failed filters, write them as result
    if (failed_filters->size > 0) {
        result = calloc (filter_text_len + 1, sizeof(char));
        strcat(result, (char*) array_list_get(0, failed_filters));
        for (int i = 1; i < failed_filters->size; i++) {
            strncat(result, ";", 1);
            strcat(result, (char*) array_list_get(i, failed_filters));
        }
    } else {
        if (pass_found) {
            result = strndup("PASS", 4);
        } else {
            result = strndup(".", 1);
        }
    }
    
    array_list_free(failed_filters, free);
    
    return result;
}


char *merge_info_field(vcf_record_file_link **position_in_files, int position_occurrences, char **info_fields, int num_fields,
                       vcf_record_t *output_record, cp_hashtable *alleles, char *empty_sample) {
    if (num_fields == 0) {
        return strndup(".", 1);
    }

    size_t len = 0;
    size_t max_len = 128;
    char *result = calloc (max_len, sizeof(char));
    
    list_t *stats_list = malloc (sizeof(list_t));
    list_init("stats", 1, INT_MAX, stats_list);
    file_stats_t *file_stats = file_stats_new();
    variant_stats_t *variant_stats = NULL;
    int dp = 0, mq0 = 0;
    double mq = 0;
    int calculate_stats = 0, calculate_dp = 0, calculate_mq = 0;
    
    // Precalculate auxiliary values
    for (int i = 0; i < num_fields; i++) {
        if (!strncmp(info_fields[i], "AC", 2) ||   // allele count in genotypes, for each ALT allele
            !strncmp(info_fields[i], "AF", 2)) {
            calculate_stats = 1;
            continue;
        }
        if (!strncmp(info_fields[i], "DP", 2) ||    // combined depth across samples
            !strncmp(info_fields[i], "QD", 2)) {
            calculate_dp = 1;
            continue;
        }
        if (!strncmp(info_fields[i], "MQ0", 3) ||    // Number of MAPQ == 0 reads covering this record
            !strncmp(info_fields[i], "MQ", 2)) {
            calculate_mq = 1;
            continue;
        }
    }
    precalculate_aux_values_for_annotation(calculate_stats, calculate_dp, calculate_mq, output_record,
                                           &variant_stats, file_stats, stats_list, &dp, &mq0, &mq);
    
    for (int i = 0; i < num_fields; i++) {
        if (len >= max_len - 32) {
            char *aux = realloc(result, max_len + 128);
            if (aux) {
                result = aux;
                max_len += 128;
            } else {
                LOG_FATAL("Can't allocate memory for file merging\n");
            }
        }
        
        // Composition of the INFO field
        
        char *field = NULL;
        size_t field_len = 0;
        
        if (!strncmp(info_fields[i], "AC", 2)) {   // allele count in genotypes, for each ALT allele
            field = get_annotation_allele_count(variant_stats, &field_len);
            
        } else if (!strncmp(info_fields[i], "AF", 2)) {   // allele frequency for each ALT allele
            field = get_annotation_allele_freq(variant_stats, &field_len);
            
        } else if (!strncmp(info_fields[i], "AN", 2)) {    // total number of alleles in called genotypes
            field = get_annotation_allele_number(variant_stats, &field_len);
            
        } else if (!strncmp(info_fields[i], "DB", 2)) {    // dbSNP membership
            for (int j = 0; j < position_occurrences; j++) {
                if (strstr(position_in_files[j]->record->info, "DB")) {
                    strncat(result, "DB;", 3);
                    len += 3;
                    break;
                }
            }
            
        } else if (!strncmp(info_fields[i], "DP", 2)) {    // combined depth across samples
            field = get_annotation_read_depth(dp, &field_len);
            
        } else if (!strncmp(info_fields[i], "H2", 2)) {    // membership in hapmap2
            for (int j = 0; j < position_occurrences; j++) {
                if (strstr(position_in_files[j]->record->info, "H2")) {
                    strncat(result, "H2;", 3);
                    len += 3;
                    break;
                }
            }
            
        } else if (!strncmp(info_fields[i], "H3", 2)) {    // membership in hapmap3
            for (int j = 0; j < position_occurrences; j++) {
                if (strstr(position_in_files[j]->record->info, "H3")) {
                    strncat(result, "H3;", 3);
                    len += 3;
                    break;
                }
            }
            
        } else if (!strncmp(info_fields[i], "MQ0", 3)) {   // Number of MAPQ == 0 reads covering this record
            field = get_annotation_mapping_quality_zero(mq0, &field_len);
            
        } else if (!strncmp(info_fields[i], "MQ", 2)) {    // RMS mapping quality
            field = get_annotation_mapping_quality(mq, &field_len);
            
        } else if (!strncmp(info_fields[i], "NS", 2)) {    // Number of samples with data
            field = get_annotation_non_missing_samples(output_record, empty_sample, &field_len);
            
        } else if (!strncmp(info_fields[i], "QD", 2)) {    // quality by depth (GATK)
            field = get_annotation_quality_by_depth(output_record, dp, &field_len);
            
        } else if (!strncmp(info_fields[i], "SOMATIC", 7)) {   // the record is a somatic mutation, for cancer genomics
            for (int j = 0; j < position_occurrences; j++) {
                if (strstr(position_in_files[j]->record->info, "SOMATIC")) {
                    strncat(result, "SOMATIC;", 8);
                    len += 8;
                    break;
                }
            }
            
        } else if (!strncmp(info_fields[i], "VALIDATED", 9)) { // validated by follow-up experiment
            for (int j = 0; j < position_occurrences; j++) {
                if (strstr(position_in_files[j]->record->info, "VALIDATED")) {
                    strncat(result, "VALIDATED;", 10);
                    len += 10;
                    break;
                }
            }
        }
        
        if (field) {
            strncat(result, field, field_len);
            strncat(result, ";", 1);
            len += field_len + 1;
            free(field);
        }
    }
    
    if (result[len-1] == ';') {
        result[len-1] = '\0';
    } else {
        result[len] = '\0';
    }
    
    if(variant_stats) {
        variant_stats_free(variant_stats);
    }
    file_stats_free(file_stats);
    free(stats_list);
    
    return result;
}


char* merge_format_field(vcf_record_file_link** position_in_files, int position_occurrences, merge_options_data_t *options, array_list_t* format_fields) {
    char *result;   
    vcf_record_t *input;
    int format_text_len = 0;
    
    // Split FORMAT of the input record and register the non-previously inserted ones
    for (int i = 0; i < position_occurrences; i++) {
        input = position_in_files[i]->record;
        int num_fields;
        char *dupformat = strndup(input->format, input->format_len);
        char **fields = split(dupformat, ":", &num_fields);
        for (int j = 0; j < num_fields; j++) {
            if (!array_list_contains(fields[j], format_fields)) {
                array_list_insert(fields[j], format_fields);
                format_text_len += strlen(fields[j]) + 1; // concat field + ":"
            } else {
                free(fields[j]);
            }
        }
        
        free(dupformat);
        free(fields);
    }
    
    if (options->copy_filter) {
        array_list_insert(strndup("SFT", 3), format_fields);
        format_text_len += options->copy_info ? 4 : 3; // If INFO is copied, another ":" is needed
    }
    if (options->copy_info) {
        array_list_insert(strndup("IN", 2), format_fields);
        format_text_len += 2;
    }
    
    result = calloc (format_text_len + 1, sizeof(char));
    strcat(result, (char*) array_list_get(0, format_fields));
    for (int i = 1; i < format_fields->size; i++) {
        strncat(result, ":", 1);
        strcat(result, (char*) array_list_get(i, format_fields));
    }
    
    return result;
}


array_list_t* merge_samples(vcf_record_file_link** position_in_files, int position_occurrences, vcf_file_t **files, int num_files, 
                            cp_hashtable* alleles_table, array_list_t* format_fields, int *format_indices, char *empty_sample,
                            int gt_pos, int filter_pos, int info_pos, merge_options_data_t *options) {
    int empty_sample_len = strlen(empty_sample);
    array_list_t *result = array_list_new(num_files * 2, 1.5, COLLECTION_MODE_ASYNCHRONIZED);
    
    for (int i = 0; i < num_files; i++) {
        vcf_record_file_link *link = NULL;
        for (int j = 0; j < position_occurrences; j++) {
            if (!strcmp(position_in_files[j]->file->filename, files[i]->filename)) {
                link = position_in_files[j];
                break;
            }
        }
        
        if (link) {
            vcf_record_t *record = link->record;
            int num_alternates, num_sample_fields;
            char *dupreference = strndup(record->reference, record->reference_len);
            char *dupalternates = strndup(record->alternate, record->alternate_len);
            char **split_alternates = split(dupalternates, ",", &num_alternates);

            char sample[1024];
            int len;
            for (int j = 0; j < get_num_vcf_samples(files[i]); j++) {
                memset(sample, 0, 256 * sizeof(char));
                len = 0;
                
                char *dupsample = strdup(array_list_get(j, record->samples));
                char **split_sample = split(dupsample, ":", &num_sample_fields);

                for (int k = 0; k < format_fields->size; k++) {
                    int idx = format_indices[i*format_fields->size + k];
//                     printf("k = %d\tidx = %d\n", k, format_indices[i*format_fields->size + k]);
                    if (idx < 0) {  // Missing
                        if (k == gt_pos) {
                            strncat(sample, "./.", 3);
                            len += 3;
                        } else if (k == filter_pos) {
                            strncat(sample, record->filter, record->filter_len);
                            len += record->filter_len;
                        } else if (k == info_pos) {
                            strncat(sample, record->info, record->info_len);
                            len += record->info_len;
                        } else {
                            strncat(sample, ".", 1);
                            len++;
                        }
                    } else {    // Not-missing
                        if (k == gt_pos) {
                            // Manage genotypes (in multiallelic variants, allele indices must be recalculated)
                            char *dupsplit = strdup(split_sample[idx]);
                            int allele1, allele2;
                            int allele_ret = get_alleles(dupsplit, 0, &allele1, &allele2);
                            free(dupsplit);
                            // TODO Haploid chromosomes such as X?
                            if (allele_ret == ALL_ALLELES_MISSING) {
                                strncat(sample, "./.", 3);
                                len += 3;
                            } else {
                                if (options->strict_reference) {
                                    concat_sample_allele_strict(allele1, len, sample, alleles_table, split_alternates);
                                    strncat(sample, split_sample[idx] + 1, 1);
                                    concat_sample_allele_strict(allele2, len, sample, alleles_table, split_alternates);
                                    len += 3;
                                } else {
                                    concat_sample_allele_tolerant(allele1, len, sample, alleles_table, dupreference, split_alternates);
                                    len++;
                                    strncat(sample, split_sample[idx] + 1, 1);
                                    len++;
                                    concat_sample_allele_tolerant(allele2, len, sample, alleles_table, dupreference, split_alternates);
                                    len++;
                                }
                            }
                        } else if (k == filter_pos) {
                            strncat(sample, record->filter, record->filter_len);
                            len += record->filter_len;
                        } else if (k == info_pos) {
                            strncat(sample, record->info, record->info_len);
                            len += record->info_len;
                        } else {
//                             printf("sample = %s\n", array_list_get(j, record->samples));
                            strcat(sample, split_sample[idx]);
                            len += strlen(split_sample[idx]);
                        }
                    }
                    
                    
                    if (k < format_fields->size - 1) {
                        strncat(sample, ":", 1);
                        len++;
                    }
                }
                
                sample[len] = '\0';
                array_list_insert(strndup(sample, len), result);
                
                for (int j = 0; j < num_sample_fields; j++) {
                    free(split_sample[j]);
                }
                free(split_sample);
                free(dupsample);

            }

            for (int j = 0; j < num_alternates; j++) {
                free(split_alternates[j]);
            }
            free(split_alternates);
            free(dupalternates);
            free(dupreference);

        } else {
            // If the file has no samples in that position, fill with empty samples
            for (int j = 0; j < get_num_vcf_samples(files[i]); j++) {
                array_list_insert(strndup(empty_sample, empty_sample_len), result);
            }
        }
    }
    
    return result;
}


/* ******************************
 *      Auxiliary functions     *
 * ******************************/

static int *get_format_indices_per_file(vcf_record_file_link **position_in_files, int position_occurrences, 
                                        vcf_file_t **files, int num_files, array_list_t *format_fields) {
    int *indices = malloc (format_fields->size * num_files * sizeof(size_t));
    int in_file = 0;
    vcf_record_t *record;
    for (int i = 0; i < num_files; i++) {
        in_file = 0;
        
        for (int j = 0; j < position_occurrences; j++) {
            // If the position is in this file, then get the record format
            assert(files[i]);
            assert(position_in_files[j]);
            if (!strcmp(position_in_files[j]->file->filename, files[i]->filename)) {
                in_file = 1;
                record = position_in_files[j]->record;
                int num_fields;
                char *dupformat = strndup(record->format, record->format_len);
                char **fields = split(dupformat, ":", &num_fields);
                
                for (int k = 0; k < format_fields->size; k++) {
                    indices[i*format_fields->size + k] = -1;
                    
                    for (int m = 0; m < num_fields; m++) {
                        if (!strcmp(array_list_get(k, format_fields), fields[m])) {
                            indices[i*format_fields->size + k] = m;
                            break;
                        }
                    }
                }
                
                for (int m = 0; m < num_fields; m++) {
                    free(fields[m]);
                }
                free(fields);
                free(dupformat);
                
                break;
            }
        }
        
        // If the position is not in the file, set all positions as -1
        if (!in_file) {
            for (int k = 0; k < format_fields->size; k++) {
                indices[i*format_fields->size + k] = -1;
            }
        }
    }
    
    return indices;
}

static inline void concat_sample_allele_tolerant(int allele, int len, char *sample, cp_hashtable* alleles_table, char *reference, char **split_alternates) {
    if (allele < 0) {
        strncat(sample, ".", 1);
    } else if (allele == 0) {
        int aux_idx = *((int*) cp_hashtable_get(alleles_table, reference));
        sprintf(sample + len, "%d", aux_idx);
    } else {
        int aux_idx = *((int*) cp_hashtable_get(alleles_table, split_alternates[allele-1]));
        sprintf(sample + len, "%d", aux_idx);
    }
}

static inline void concat_sample_allele_strict(int allele, int len, char *sample, cp_hashtable* alleles_table, char **split_alternates) {
    if (allele < 0) {
        strncat(sample, ".", 1);
    } else if (allele == 0) {
        strncat(sample, "0", 1);
    } else {
        int aux_idx = *((int*) cp_hashtable_get(alleles_table, split_alternates[allele-1]));
        sprintf(sample + len, "%d", aux_idx);
    }
}

static char *get_empty_sample(int num_format_fields, int gt_pos, merge_options_data_t *options) {
    assert(options);
    int sample_len = num_format_fields * 2 + 2; // Each field + ':' = 2 chars, except for GT which is 1 char more
    char *sample = (char*) calloc (sample_len, sizeof(char));
    for (int j = 0; j < num_format_fields; j++) {
        if (j > 0) {
            strncat(sample, ":", 1);
        }
        if (j != gt_pos) {
            strncat(sample, ".", 1);
        } else {
            if (options->missing_mode == MISSING) {
                strncat(sample, "./.", 3);
            } else if (options->missing_mode == REFERENCE) {
                strncat(sample, "0/0", 3);
            }
        }
    }
    return sample;
}

vcf_record_file_link *vcf_record_file_link_new(vcf_record_t *record, vcf_file_t *file) {
    vcf_record_file_link *link = malloc(sizeof(vcf_record_file_link));
    link->record = record;
    link->file = file;
    return link;
}

void vcf_record_file_link_free(vcf_record_file_link *link) {
    assert(link);
    vcf_record_free_deep(link->record);
    free(link);
}
