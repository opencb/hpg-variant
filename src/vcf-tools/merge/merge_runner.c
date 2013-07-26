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

#include "merge_runner.h"

#define TREE_LIMIT         (shared_options_data->batch_lines)

static int num_chromosomes;
static char **chromosome_order;

static void free_merge_tree(kh_pos_t* positions_read);


int run_merge(shared_options_data_t *shared_options_data, merge_options_data_t *options_data) {
    if (options_data->num_files == 1) {
        LOG_INFO("Just one VCF file specified, no need to merge\n");
        return 0;
    }
    
    list_t *read_list[options_data->num_files];
    memset(read_list, 0, options_data->num_files * sizeof(list_t*));
    list_t *output_header_list = (list_t*) malloc (sizeof(list_t));
    list_init("headers", shared_options_data->num_threads, INT_MAX, output_header_list);
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, shared_options_data->max_batches * shared_options_data->batch_lines, output_list);
    list_t *merge_tokens_list = (list_t*) malloc (sizeof(list_t));
    list_init("tokens", 1, INT_MAX, merge_tokens_list);
    
    int ret_code = 0;
    double start, stop, total;
    vcf_file_t *files[options_data->num_files];
    memset(files, 0, options_data->num_files * sizeof(vcf_file_t*));
    
    // Initialize variables related to the different files
    for (int i = 0; i < options_data->num_files; i++) {
        files[i] = vcf_open(options_data->input_files[i], shared_options_data->max_batches);
        if (!files[i]) {
            LOG_FATAL_F("VCF file %s does not exist!\n", options_data->input_files[i]);
        }
        
        read_list[i] = (list_t*) malloc(sizeof(list_t));
        list_init("text", 1, shared_options_data->max_batches, read_list[i]);
    }
    
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }

    chromosome_order = get_chromosome_order(shared_options_data->host_url, shared_options_data->species,
                                            shared_options_data->version, &num_chromosomes);
    
#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            ret_code = vcf_multiread_batches(read_list, shared_options_data->batch_lines, files, options_data->num_files);

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) {
                LOG_ERROR_F("Error %d while reading VCF files\n", ret_code);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        }
        
#pragma omp section
        {
            // Enable nested parallelism
            omp_set_nested(1);
            
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            int num_eof_found = 0;
            int eof_found[options_data->num_files];
            memset(eof_found, 0, options_data->num_files * sizeof(int));
            
            list_item_t *items[options_data->num_files];
            memset(items, 0, options_data->num_files * sizeof(list_item_t*));
            char *texts[options_data->num_files];
            memset(texts, 0, options_data->num_files * sizeof(char*));
            
            // Positions read and to merge
            khash_t(pos) *positions_read = kh_init(pos);
            // Last chromosome and positions read from each file (to block file reading when necessary)
            char *last_chromosome_read[options_data->num_files];
            long last_position_read[options_data->num_files];
            memset(last_chromosome_read, 0, options_data->num_files * sizeof(char*));
            memset(last_position_read, 0, options_data->num_files * sizeof(long));
            // Last chromosome and position merged
            char *max_chromosome_merged = NULL;
            long max_position_merged = LONG_MAX;
            // Whether each file was ahead the last merged position
            int file_ahead_last_merged[options_data->num_files];
            memset(file_ahead_last_merged, 0, options_data->num_files * sizeof(int));
            // Whether the headers of all files have been merged
            int header_merged = 0;
            // Token returned from the merging process: notifies that the output 
            // can be written to file and the number of records to write
            int token = 0;
            
            
            start = omp_get_wtime();

            while (num_eof_found < options_data->num_files) {
                /* Process:
                 * - N threads getting batches of VCF records and inserting them in a data structure. The common minimum 
                 * position of each group of batches will also be stored.
                 * - If the data structure reaches certain size or the end of a chromosome, merge positions prior to the 
                 * last minimum registered.
                 */
                
                // Getting text elements in a critical region guarantees that each thread gets variants in positions in the same range
                for (int i = 0; i < options_data->num_files; i++) {
                    if (eof_found[i] || file_ahead_last_merged[i] > 0) {
                        continue;
                    }
                    
                    items[i] = list_remove_item(read_list[i]);
                    if (items[i] == NULL || !strcmp(items[i]->data_p, "")) {
                        LOG_INFO_F("[%d] EOF found in file %s\n", omp_get_thread_num(), options_data->input_files[i]);
                        // Mark as finished
                        eof_found[i] = 1;
                        num_eof_found++;
                        // Set last chromosome and position to the maximum possible so they do not disturb pending files
                        last_chromosome_read[i] = chromosome_order[num_chromosomes-1];
                        last_position_read[i] = LONG_MAX;
                        
                        if(items[i] != NULL && !strcmp(items[i]->data_p, "")) {
                            free(items[i]->data_p);
                            list_item_free(items[i]);
                            LOG_DEBUG_F("[%d] Text batch freed\n", omp_get_thread_num());
                        } else {
                            LOG_DEBUG_F("[%d] No need to free text batch\n", omp_get_thread_num());
                        }
                        
                        continue;
                    }
                    
                    assert(items[i]->data_p != NULL);
                    texts[i] = items[i]->data_p;
                    
//                     printf("[%d] text batch from file %d\tcontents = '%s'\n", omp_get_thread_num(), i, texts[i]);
                }
                
                for (int i = 0; i < options_data->num_files; i++) {
                    if (eof_found[i] || file_ahead_last_merged[i] > 0) {
                        continue;
                    }
                    
                    char *text_begin = texts[i];
                    char *text_end = text_begin + strlen(text_begin);
                    assert(text_end != NULL);
                    
//                     printf("batch = '%.*s'\n", text_end - text_begin, text_begin);
                    
                    // Get VCF batches from text batches
                    vcf_reader_status *status = vcf_reader_status_new(shared_options_data->batch_lines, 0);
                    ret_code = run_vcf_parser(text_begin, text_end, shared_options_data->batch_lines, files[i], status);
                    
                    if (ret_code) {
                        // TODO stop?
                        LOG_ERROR_F("Error %d while reading the file %s\n", ret_code, files[i]->filename);
                        continue;
                    }

//                     printf("batches = %d\n", files[i]->record_batches->length);
                    vcf_batch_t *batch = fetch_vcf_batch_non_blocking(files[i]);
                    if (!batch) {
                        continue;
                    }
                    
                    // Insert records into hashtable
                    for (int j = 0; j < batch->records->size; j++) {
                        vcf_record_t *record = vcf_record_copy(array_list_get(j, batch->records));
                        vcf_record_file_link *link = vcf_record_file_link_new(record, files[i]);
                        char key[64];
                        compose_key_value(record->chromosome, record->position, key);
                        int ret = insert_position_read(key, link, positions_read);
                        assert(ret);
                    }
                    
                    vcf_record_t *current_record = (vcf_record_t*) array_list_get(batch->records->size - 1, batch->records);
                    //int updated = calculate_merge_interval(current_record, &max_chromosome_merged, &max_position_merged, chromosome_order, num_chromosomes);
                    last_chromosome_read[i] = strndup(current_record->chromosome, current_record->chromosome_len);
                    last_position_read[i] = current_record->position;
                    
                    // Free batch and its contents
                    vcf_reader_status_free(status);
                    vcf_batch_free(batch);
                    list_item_free(items[i]);
                }
                
                // Calculate least common position to merge
                for (int i = 0; i < options_data->num_files; i++) {
                    char *file_chromosome = last_chromosome_read[i];
                    long file_position = last_position_read[i];
                    if (max_chromosome_merged == NULL) {
                        // Max merged chrom:position not set, assign without any other consideration
                        max_chromosome_merged = strdup(file_chromosome);
                        max_position_merged = file_position;
                        //return 1;
                    } else {
                        char *current_chromosome = file_chromosome;
                        long unsigned int current_position = file_position;

                //         printf("current = %s:%ld\tmax = %s:%ld\n", current_chromosome, current_position, *max_chromosome_merged, *max_position_merged);
                        int chrom_comparison = compare_chromosomes(current_chromosome, max_chromosome_merged, chromosome_order, num_chromosomes);
                        int position_comparison = compare_positions(current_position, max_position_merged);

                        // Max merged chrom:position is greater than the last one in this batch
                        if (chrom_comparison < 0 || (chrom_comparison == 0 && position_comparison < 0)) {
                            max_chromosome_merged = strdup(current_chromosome);
                            max_position_merged = current_position;
                            //return 1;
                        }
                    }
                }
                
                if (num_eof_found == options_data->num_files) {
                    max_chromosome_merged = chromosome_order[num_chromosomes-1];
                    max_position_merged = LONG_MAX;
                }
                
                LOG_DEBUG_F("Last position is (%s, %ld). Files ahead: ", max_chromosome_merged, max_position_merged);
                // Check which files are ahead the last position merged and will not be read in the next iteration
                for (int i = 0; i < options_data->num_files; i++) {
                    char *file_chromosome = last_chromosome_read[i];
                    long file_position = last_position_read[i];
                    
                    int chrom_comparison = compare_chromosomes(file_chromosome, max_chromosome_merged, chromosome_order, num_chromosomes);
                    int position_comparison = compare_positions(file_position, max_position_merged);
                    
                    if (chrom_comparison < 0 || (chrom_comparison == 0 && position_comparison <= 0)) {
                        file_ahead_last_merged[i] = 0;
                    } else {
                        file_ahead_last_merged[i] = 1;
                        LOG_DEBUG_F("%d (%s, %ld), ", i, file_chromosome, file_position);
                    }
                }
                LOG_DEBUG("\n");
                
                // Merge headers, if not previously done
                if (!header_merged) {
                    // Check correction of input file headers
                    array_list_t *sample_names = merge_vcf_sample_names(files, options_data->num_files);
                    if (!sample_names) {
                        // Avoid as many leaks as possible
                        array_list_free(sample_names, NULL);
                        free_merge_tree(positions_read);
                        free(max_chromosome_merged);
                        
                        LOG_FATAL("Files can not be merged!\n");
                    } else {
                        array_list_free(sample_names, NULL);
                    }
                    
                    // Run the merge itself
                    merge_vcf_headers(files, options_data->num_files, options_data, output_header_list);
                    header_merged = 1;
                    
                    // Decrease list writers count
                    for (int i = 0; i < shared_options_data->num_threads; i++) {
                        list_decr_writers(output_header_list);
                    }
                }
                
                // If the data structure reaches certain size or the end of a chromosome, 
                // merge positions prior to the last minimum registered
                if (num_eof_found < options_data->num_files && kh_size(positions_read) > TREE_LIMIT) {
                    LOG_INFO_F("Merging until position %s:%ld\n", max_chromosome_merged, max_position_merged);
                    token = merge_interval(positions_read, max_chromosome_merged, max_position_merged, chromosome_order, num_chromosomes,
                                   	   	   files, shared_options_data, options_data, output_list);
                }
                // When reaching EOF for all files, merge the remaining entries
                else if (num_eof_found == options_data->num_files && kh_size(positions_read) > 0) {
                    LOG_INFO_F("Merging remaining positions (last = %s:%ld)\n", chromosome_order[num_chromosomes - 1], LONG_MAX);
                    token = merge_remaining_interval(positions_read, files, shared_options_data, options_data, output_list);
                }
                
                if (token) {
                    int *token_ptr = malloc (sizeof(int)); *token_ptr = token;
                    list_item_t *item = list_item_new(1, 0, token_ptr);
                    list_insert_item(item, merge_tokens_list);
                }

                // Set variables ready for next iteration of the algorithm
                if (max_chromosome_merged) {
                    free(max_chromosome_merged);
                }
            	token = 0;
                max_chromosome_merged = NULL;
                max_position_merged = LONG_MAX;
            }
            
            kh_destroy(pos, positions_read);
            
            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Decrease list writers count
            for (int i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
            list_decr_writers(merge_tokens_list);
        }
        
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d writes the output\n", omp_get_thread_num());
    
            start = omp_get_wtime();

            // Create file streams for results
            char aux_filename[32]; memset(aux_filename, 0, 32 * sizeof(char));
            sprintf(aux_filename, "merge_from_%d_files.vcf", options_data->num_files);
            
            char *merge_filename;
            FILE *merge_fd = get_output_file(shared_options_data, aux_filename, &merge_filename);
            LOG_INFO_F("Output filename = %s\n", merge_filename);
            free(merge_filename);
            
            list_item_t *item1 = NULL, *item2 = NULL;
            vcf_header_entry_t *entry;
            vcf_record_t *record;
            int *num_records;
            
            // Write headers
            while ((item1 = list_remove_item(output_header_list))) {
                entry = item1->data_p;
                write_vcf_header_entry(entry, merge_fd);
                list_item_free(item1);
            }
            
            // Write delimiter
            array_list_t *sample_names = merge_vcf_sample_names(files, options_data->num_files);
            write_vcf_delimiter_from_samples((char**) sample_names->items, sample_names->size, merge_fd);
            
            // Write records
            // When a token is present, it means a set of batches has been merged. The token contains the number of records merged.
            // In this case, the records must be sorted by chromosome and position, and written afterwards.
            while ((item1 = list_remove_item(merge_tokens_list))) {
                num_records = item1->data_p;
                vcf_record_t *records[*num_records];
                for (int i = 0; i < *num_records; i++) {
                    item2 = list_remove_item(output_list);
                    if (!item2) {
                        break;
                    }

                    records[i] = item2->data_p;
                    list_item_free(item2);
                }

                // Sort records
                qsort(records, *num_records, sizeof(vcf_record_t*), record_cmp);

                // Write and free sorted records
                for (int i = 0; i < *num_records; i++) {
                    record = records[i];
                    write_vcf_record(record, merge_fd);
                    vcf_record_free_deep(record);
                }

                free(num_records);
                list_item_free(item1);
            }
            
            // Close file
            array_list_free(sample_names, NULL);
            if (merge_fd != NULL) { fclose(merge_fd); }
            
            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        }
    }

    // Free variables related to the different files
    for (int i = 0; i < options_data->num_files; i++) {
        if(files[i]) { vcf_close(files[i]); }
        if(read_list[i]) { free(read_list[i]); }
    }
    free(output_list);
    free(output_header_list);
    free(merge_tokens_list);
    
    return ret_code;
}



int insert_position_read(char key[64], vcf_record_file_link* link, kh_pos_t* positions_read) {
    int ret;
    array_list_t *records_in_position;
    khiter_t iter = kh_get(pos, positions_read, key);
    if (iter != kh_end(positions_read)) {
        records_in_position = kh_value(positions_read, iter);
        ret = array_list_insert(link, records_in_position);
    } else {
        records_in_position = array_list_new(8, 1.5, COLLECTION_MODE_SYNCHRONIZED);
        ret = array_list_insert(link, records_in_position);
        iter = kh_put(pos, positions_read, strdup(key), &ret);
        if (ret) {
            kh_value(positions_read, iter) = records_in_position;
        }
    }
    
    return ret;
}


int calculate_merge_interval(vcf_record_t* current_record, char** max_chromosome_merged, long unsigned int* max_position_merged,
                              char **chromosome_order, int num_chromosomes) {
    if (*max_chromosome_merged == NULL) {
        // Max merged chrom:position not set, assign without any other consideration
        *max_chromosome_merged = strndup(current_record->chromosome, current_record->chromosome_len);
        *max_position_merged = current_record->position;
        return 1;
    } else {
        char *current_chromosome = strndup(current_record->chromosome, current_record->chromosome_len);
        long unsigned int current_position = current_record->position;
        
//         printf("current = %s:%ld\tmax = %s:%ld\n", current_chromosome, current_position, *max_chromosome_merged, *max_position_merged);
        int chrom_comparison = compare_chromosomes(current_chromosome, *max_chromosome_merged, chromosome_order, num_chromosomes);
        int position_comparison = compare_positions(current_position, *max_position_merged);
        
        // Max merged chrom:position is greater than the last one in this batch
        if (chrom_comparison < 0 || (chrom_comparison == 0 && position_comparison < 0)) {
            *max_chromosome_merged = current_chromosome;
            *max_position_merged = current_position;
            return 1;
        } else {
            assert(current_chromosome);
            free(current_chromosome);
            return 0;
        }
    }
}


int merge_interval(kh_pos_t* positions_read, char *max_chromosome_merged, unsigned long max_position_merged,
                    char **chromosome_order, int num_chromosomes, vcf_file_t **files, 
                    shared_options_data_t *shared_options_data, merge_options_data_t *options_data, list_t *output_list) {
    int num_entries = 0;

    #pragma omp parallel for num_threads(shared_options_data->num_threads) reduction(+:num_entries)
    for (int k = kh_begin(positions_read); k < kh_end(positions_read); k++) {
        if (kh_exist(positions_read, k)) {
            array_list_t *records_in_position = kh_value(positions_read, k);
            assert(records_in_position);
            
            vcf_record_t *record = ((vcf_record_file_link*) array_list_get(0, records_in_position))->record;
            vcf_record_file_link **links = NULL;
            int num_links = 0;
            
            // Remove positions prior to the last chromosome:position to merge
            int cmp_chrom = compare_chromosomes(record->chromosome, max_chromosome_merged, chromosome_order, num_chromosomes);
            if (cmp_chrom < 0 || (cmp_chrom == 0 && compare_positions(record->position, max_position_merged) <= 0)) {
                links = (vcf_record_file_link**) records_in_position->items;
                num_links = records_in_position->size;
            }
            
            // Launch merge
            if (num_links > 0) {
//                 printf("links[0] = %s:%ld in file %s\n", links[0]->record->chromosome, links[0]->record->position, links[0]->file->filename);
                int err_code = 0;
                vcf_record_t *merged = merge_position(links, num_links, files, options_data->num_files, options_data, &err_code);
                
                if (!err_code) {
                    list_item_t *item = list_item_new(k, MERGED_RECORD, merged);
                    list_insert_item(item, output_list);
                    num_entries += 1;
                }
                
                // Free empty nodes (lists of records in the same position)
                free(kh_key(positions_read, k));
                array_list_free(records_in_position, vcf_record_file_link_free);
                kh_del(pos, positions_read, k);
            }
        } // End kh_exist
    }

    return num_entries;
}


int merge_remaining_interval(kh_pos_t* positions_read, vcf_file_t **files, shared_options_data_t *shared_options_data,
                              merge_options_data_t *options_data, list_t *output_list) {
    int num_entries = 0;

    #pragma omp parallel for num_threads(shared_options_data->num_threads) reduction(+:num_entries)
    for (int k = kh_begin(positions_read); k < kh_end(positions_read); k++) {
        if (kh_exist(positions_read, k)) {
            array_list_t *records_in_position = kh_value(positions_read, k);
            assert(records_in_position);
            
            // Launch merge
            int err_code = 0;
            vcf_record_t *merged = merge_position((vcf_record_file_link **) records_in_position->items, records_in_position->size, 
                                                  files, options_data->num_files, options_data, &err_code);
            
            if (!err_code) {
                list_item_t *item = list_item_new(k, MERGED_RECORD, merged);
                list_insert_item(item, output_list);
                num_entries += 1;
            }
            
            // Free empty nodes (lists of records in the same position)
            free(kh_key(positions_read, k));
            array_list_free(records_in_position, vcf_record_file_link_free);
            kh_del(pos, positions_read, k);
        }
    }

    return num_entries;
}



static void compose_key_value(const char *chromosome, const long position, char *key) {
    assert(key);
    strcpy(key, chromosome);
    sprintf(key+strlen(key), "_%ld", position);
}

static int record_cmp(const void *data1, const void *data2) {
	vcf_record_t **record1 = (vcf_record_t **) data1;
	vcf_record_t **record2 = (vcf_record_t **) data2;

	int cmp = 0;
	for (int i = 0; i < num_chromosomes; i++) {
		assert(chromosome_order[i]);
		if (!strncmp(chromosome_order[i], (*record1)->chromosome, (*record1)->chromosome_len)) {
			if (!strncmp(chromosome_order[i], (*record2)->chromosome, (*record2)->chromosome_len)) {
				cmp = 0;
				break;
			}
			cmp = -1;
			break;
		} else if (!strncmp(chromosome_order[i], (*record2)->chromosome, (*record2)->chromosome_len)) {
			cmp = 1;
			break;
		}
	}

	if (!cmp) {
		cmp = (*record1)->position - (*record2)->position;
	}

	return cmp;
}

static void free_merge_tree(kh_pos_t* positions_read) {
    for (int k = kh_begin(positions_read); k < kh_end(positions_read); k++) {
        if (kh_exist(positions_read, k)) {
            array_list_free(kh_value(positions_read, k), vcf_record_file_link_free);
            free(kh_key(positions_read, k));
            kh_del(pos, positions_read, k);
        }
    }
}