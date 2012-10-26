/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
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
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, "vcf-info-fields.cfg");
    
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
    array_list_t *samples = array_list_new(get_num_vcf_samples(files[0]) * 2, 2, COLLECTION_MODE_ASYNCHRONIZED);
    for (int i = 0; i < num_files; i++) {
        array_list_insert_all(files[i]->samples_names->items, files[i]->samples_names->size, samples);
    }
    return samples;
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
    vcf_file_t *file;
    vcf_record_t *input;
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
            LOG_ERROR_F("Position %.*s:%ld can't be merged: Discordant reference alleles (%.*s, %.*s)\n", 
                        position_in_files[i]->record->chromosome_len, position_in_files[i]->record->chromosome, position_in_files[i]->record->position,
                        position_in_files[i]->record->reference_len, position_in_files[i]->record->reference, 
                        position_in_files[i+1]->record->reference_len, position_in_files[i+1]->record->reference);
            *err_code = DISCORDANT_REFERENCE;
            return NULL;
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
    char *alternates = NULL, *cur_alternate = NULL, *aux = NULL;
    char *pos_alternates = NULL, **split_pos_alternates = NULL;
    int num_pos_alternates;
    
    int cur_index = 0;
    int *allele_index = (int*) calloc (1, sizeof(int)); *allele_index = cur_index;
    cp_hashtable_put(alleles_table,
                     strndup(position_in_files[0]->record->reference, position_in_files[0]->record->reference_len), 
                     allele_index);
    
    for (int i = 0; i < position_occurrences; i++) {
        file = position_in_files[i]->file;
        input = position_in_files[i]->record;

    	position = input->position;
        pos_alternates = strndup(input->alternate, input->alternate_len);
        split_pos_alternates = split(pos_alternates, ",", &num_pos_alternates);

        for (int j = 0; j < num_pos_alternates; j++) {
        	cur_alternate = split_pos_alternates[j];

			if (!alternates) {
				alternates = strdup(cur_alternate); // Need to be dup because it will be inserted (and destroyed) in the alleles table
				max_len = input->alternate_len;

				allele_index = (int*) calloc (1, sizeof(int));
				*allele_index = ++cur_index;
				cp_hashtable_put(alleles_table, cur_alternate, allele_index);
			} else {
				if (!cp_hashtable_contains(alleles_table, cur_alternate)) {
					concat_len = strlen(cur_alternate);
					aux = realloc(alternates, max_len + concat_len + 2);
					if (aux) {
						// Concatenate alternate value to the existing list
						strncat(aux, ",", 1);
	//                     printf("1) cur_alternate = %.*s\naux = %s\n---------\n", concat_len, cur_alternate, aux);
						strncat(aux, cur_alternate, concat_len);
						alternates = aux;
						max_len += concat_len + 1;
	//                     printf("2) alternates = %s\n---------\n", alternates);

						// In case the allele is not in the hashtable, insert it with a new index
						allele_index = (int*) calloc (1, sizeof(int));
						*allele_index = ++cur_index;
						cp_hashtable_put(alleles_table, cur_alternate, allele_index);
					} else {
						LOG_FATAL_F("Can't allocate memory for alternate alleles in position %s:%ld\n",
									input->chromosome, input->position);
					}
				}
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
            char *filter = strndup(input->filter, input->filter_len);
            if (!array_list_contains(filter, failed_filters)) {
                array_list_insert(filter, failed_filters);
                filter_text_len += strlen(filter) + 1; // concat field + ","
            }
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
    size_t len = 0;
    size_t max_len = 128;
    char *result = calloc (max_len, sizeof(char));
    char *aux;
    
    list_t *stats_list = malloc (sizeof(list_t));
    list_init("stats", 1, INT_MAX, stats_list);
    file_stats_t *file_stats = file_stats_new();
    variant_stats_t *variant_stats = NULL;
    int dp = 0, mq0 = 0;
    int dp_checked = 0, mq_checked = 0, stats_checked = 0;
    double mq = 0;
    
    for (int i = 0; i < num_fields; i++) {
//         printf("\n----------\nfield = %s\nresult = %s\n--------\n", info_fields[i], result);
        if (len >= max_len - 32) {
            aux = realloc(result, max_len + 128);
            if (aux) {
                result = aux;
                max_len += 128;
            } else {
                LOG_FATAL("Can't allocate memory for file merging");
            }
        }
        
        // Conditional precalculations
        // TODO will it be faster to calculate these once and for all or keep asking strncmp?
        
        if (!stats_checked && 
            (!strncmp(info_fields[i], "AC", 2) ||   // allele count in genotypes, for each ALT allele
             !strncmp(info_fields[i], "AF", 2))) {  // allele frequency for each ALT allele
            get_variants_stats(&output_record, 1, stats_list, file_stats);
            variant_stats = list_remove_item(stats_list)->data_p;
        }
        
        if (!dp_checked && 
            (!strncmp(info_fields[i], "DP", 2) ||    // combined depth across samples
             !strncmp(info_fields[i], "QD", 2))) {   // quality by depth (GATK)
            int dp_pos = -1;
            for (int j = 0; j < output_record->samples->size; j++) {
                dp_pos = get_field_position_in_format("DP", strndup(output_record->format, output_record->format_len));
                if (dp_pos >= 0) {
                    dp += atoi(get_field_value_in_sample(strdup((char*) array_list_get(j, output_record->samples)), dp_pos));
                }
            }
            dp_checked = 1;
        }
        
        if (!mq_checked &&
            (!strncmp(info_fields[i], "MQ0", 3) ||    // Number of MAPQ == 0 reads covering this record
             !strncmp(info_fields[i], "MQ", 2))) {    // RMS mapping quality
            int mq_pos;
            int cur_gq;
            for (int j = 0; j < output_record->samples->size; j++) {
                mq_pos = get_field_position_in_format("GQ", strndup(output_record->format, output_record->format_len));
                if (mq_pos < 0) {
                    continue;
                }
                
                cur_gq = atoi(get_field_value_in_sample(strdup((char*) array_list_get(j, output_record->samples)), mq_pos));
//                 printf("sample = %s\tmq_pos = %d\tvalue = %d\n", array_list_get(j, record->samples), mq_pos,
//                        atoi(get_field_value_in_sample(strdup((char*) array_list_get(j, record->samples)), mq_pos)));
                if (cur_gq == 0) {
                    mq0++;
                } else {
                    mq += cur_gq * cur_gq;
                }
            }
            mq = sqrt(mq / output_record->samples->size);
            mq_checked = 1;
        }
        
        // Composition of the INFO field
        
        if (!strncmp(info_fields[i], "AC", 2)) {   // allele count in genotypes, for each ALT allele
            strncat(result, "AC=", 3);
            len += 3;
            for (int j = 1; j < variant_stats->num_alleles; j++) {
                if (j < variant_stats->num_alleles - 1) {
                    sprintf(result+len, "%d,", variant_stats->alleles_count[j]);
                } else {
                    sprintf(result+len, "%d;", variant_stats->alleles_count[j]);
                }
                len = strlen(result);
            }
            
        } else if (!strncmp(info_fields[i], "AF", 2)) {   // allele frequency for each ALT allele
            // For each ALT, AF = alt_freq / (1 - ref_freq)
            strncat(result, "AF=", 3);
            len += 3;
            for (int j = 1; j < variant_stats->num_alleles; j++) {
                if (j < variant_stats->num_alleles - 1) {
                    sprintf(result+len, "%.3f,", variant_stats->alleles_freq[j] / (1 - variant_stats->alleles_freq[0]));
                } else {
                    sprintf(result+len, "%.3f;", variant_stats->alleles_freq[j] / (1 - variant_stats->alleles_freq[0]));
                }
                len = strlen(result);
            }
            
        } else if (!strncmp(info_fields[i], "AN", 2)) {    // total number of alleles in called genotypes
            sprintf(result+len, "AN=%ld;", cp_hashtable_count(alleles));
            len = strlen(result);
            
        } else if (!strncmp(info_fields[i], "DB", 2)) {    // dbSNP membership
            for (int j = 0; j < position_occurrences; j++) {
                if (strstr(position_in_files[j]->record->info, "DB")) {
                    strncat(result, "DB;", 3);
                    len += 3;
                    break;
                }
            }
            
        } else if (!strncmp(info_fields[i], "DP", 2)) {    // combined depth across samples
            sprintf(result+len, "DP=%d;", dp);
            len = strlen(result);
            
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
            sprintf(result+len, "MQ0=%d;", mq0);
            len = strlen(result);
            
        } else if (!strncmp(info_fields[i], "MQ", 2)) {    // RMS mapping quality
            sprintf(result+len, "MQ=%.3f;", mq);
            len = strlen(result);
            
        } else if (!strncmp(info_fields[i], "NS", 2)) {    // Number of samples with data
            int ns = 0;
            for (int j = 0; j < output_record->samples->size; j++) {
                if (strcmp(array_list_get(j, output_record->samples), empty_sample)) {
                    ns++;
                }
            }
            sprintf(result+len, "NS=%d;", ns);
            len = strlen(result);
            
        } else if (!strncmp(info_fields[i], "QD", 2)) {    // quality by depth (GATK)
            if (output_record->quality < 0) {
                strncat(result, "QD=.;", 1);
                len += 5;
            } else if (output_record->quality > 0) {
                sprintf(result+len, "QD=%.3f;", output_record->quality / dp);
                len = strlen(result);
            } else {
                strncat(result, "QD=0.000;", 1);
                len += 5;
            }
            
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
                if (strstr(position_in_files[j]->record->info, "SOMATIC")) {
                    strncat(result, "VALIDATED;", 10);
                    len += 10;
                    break;
                }
            }
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
            }
        }
        
        free(dupformat);
        free(fields);
    }
    
    if (options->copy_filter) {
        array_list_insert(strndup("SFT", 3), format_fields);
        format_text_len += 3;
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
    int aux_idx;
    
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
                            if (allele_ret == 3) {
                                strncat(sample, "./.", 3);
                                len += 3;
                            } else {
                                if (allele1 < 0) {
                                    strncat(sample, ".", 1);
                                } else if (allele1 == 0) {
                                    strncat(sample, "0", 1);
                                } else {
//                                     printf("alternate = %.*s\nidx = %d\n", 
//                                            record->alternate_len, record->alternate, 
//                                            *((int*) cp_hashtable_get(alleles_table, record->alternate)));
                                    aux_idx = *((int*) cp_hashtable_get(alleles_table, split_alternates[allele1-1]));
                                    sprintf(sample + len, "%d", aux_idx);
                                }
                                len++;
                                
                                strncat(sample, split_sample[idx] + 1, 1);
                                len++;
                                
                                if (allele2 < 0) {
                                    strncat(sample, ".", 1);
                                } else if (allele2 == 0) {
                                    strncat(sample, "0", 1);
                                } else {
                                    aux_idx = *((int*) cp_hashtable_get(alleles_table, split_alternates[allele2-1]));
                                    sprintf(sample + len, "%d", aux_idx);
                                }
                                len++;
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

int *get_format_indices_per_file(vcf_record_file_link **position_in_files, int position_occurrences, 
                                    vcf_file_t **files, int num_files, array_list_t *format_fields) {
    int *indices = malloc (format_fields->size * num_files * sizeof(size_t));
    int in_file = 0;
    vcf_record_t *record;
    char *split_format;
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

char *get_empty_sample(int num_format_fields, int gt_pos, merge_options_data_t *options) {
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
