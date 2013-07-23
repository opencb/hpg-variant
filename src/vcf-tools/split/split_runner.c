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

#include "split_runner.h"


int run_split(shared_options_data_t *shared_options_data, split_options_data_t *options_data) {
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, MIN(10, shared_options_data->max_batches) * shared_options_data->batch_lines, output_list);
    cp_hashtable *output_files = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                                               50,
                                                               cp_hash_istring,
                                                               (cp_compare_fn) strcasecmp,
                                                               NULL,
                                                               (cp_destructor_fn) free_file_key,
                                                               NULL,
                                                               (cp_destructor_fn) free_file_descriptor
                                                              );
    
    int ret_code = 0;
    double start, stop, total;
    vcf_file_t *file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
    
    if (!file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }
    
#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            if (shared_options_data->batch_bytes > 0) {
                ret_code = vcf_parse_batches_in_bytes(shared_options_data->batch_bytes, file);
            } else if (shared_options_data->batch_lines > 0) {
                ret_code = vcf_parse_batches(shared_options_data->batch_lines, file);
            }

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) {
                LOG_ERROR_F("Error %d while reading the file %s\n", ret_code, file->filename);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            notify_end_parsing(file);
        }
        
#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            omp_set_num_threads(shared_options_data->num_threads);
            
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            start = omp_get_wtime();

            int i = 0;
            vcf_batch_t *batch = NULL;
            while ((batch = fetch_vcf_batch(file)) != NULL) {
//                 vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
                array_list_t *input_records = batch->records;

                if (i % 50 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->records->size, batch->records->capacity);
                }

                // Divide the list of passed records in ranges of size defined in config file
                int num_chunks;
                int *chunk_sizes;
                int *chunk_starts = create_chunks(input_records->size, 
                                                  ceil((float) shared_options_data->batch_lines / shared_options_data->num_threads), 
                                                  &num_chunks, &chunk_sizes);
                
                // OpenMP: Launch a thread for each range
                #pragma omp parallel for
                for (int j = 0; j < num_chunks; j++) {
                    LOG_DEBUG_F("[%d] Split invocation\n", omp_get_thread_num());
                    if (options_data->criterion == SPLIT_CHROMOSOME) {
                        ret_code = split_by_chromosome((vcf_record_t**) (input_records->items + chunk_starts[j]), 
                                                       chunk_sizes[j],
                                                       output_list);
                    } else if (options_data->criterion == SPLIT_COVERAGE) {
                        ret_code = split_by_coverage((vcf_record_t**) (input_records->items + chunk_starts[j]), chunk_sizes[j], 
						     options_data->intervals, options_data->num_intervals, output_list);
                    }
                }
//                 if (i % 50 == 0) { LOG_INFO_F("*** %dth split invocation finished\n", i); }
                
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
                list_decr_writers(output_list);
            }
        }
        
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d writes the output\n", omp_get_thread_num());
    
            start = omp_get_wtime();

            // Create file streams for results
            list_item_t* item = NULL;
            split_result_t *result;
            FILE *split_fd = NULL;
            char split_filename[1024];
            char input_filename[256];
            get_filename_from_path(shared_options_data->vcf_filename, input_filename);
            
            while ((item = list_remove_item(output_list)) != NULL) {
                result = item->data_p;
                
                sprintf(split_filename, "%s/%s_%s", shared_options_data->output_directory, result->split_name, input_filename);
//                 printf("Split filename = '%s'\n", split_filename);
                
                split_fd = cp_hashtable_get(output_files, result->split_name);
                if (!split_fd) {
                    // If this is the first line to be written to the file, create the file descriptor...
                    split_fd = fopen(split_filename, "w");
                    cp_hashtable_put(output_files, strdup(result->split_name), split_fd);
                    // ...and insert the header
                    write_vcf_header(file, split_fd);
                }
                
                // Write line into the file
                write_vcf_record(result->record, split_fd);
                vcf_record_free_deep(result->record);
                
                free_split_result(result);
                list_item_free(item);
            }
            
            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        }
    }

    free_output(output_files);
    free(output_list);
    vcf_close(file);
    
    return ret_code;
}



static int initialize_output(cp_hashtable **output_files) {
    // Initialize collections of file descriptors and summary counters
    *output_files = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                                 50,
                                                 cp_hash_istring,
                                                 (cp_compare_fn) strcasecmp,
                                                 NULL,
                                                 (cp_destructor_fn) free_file_key,
                                                 NULL,
                                                 (cp_destructor_fn) free_file_descriptor
                                                );
                    
    return 0;
}

static void free_output(cp_hashtable *output_files) {
    // Free file descriptors and summary counters
    cp_hashtable_destroy(output_files);
}

static void free_file_key(char *key) {
    LOG_DEBUG_F("Free file key: %s\n", key);
    free(key);
}

static void free_file_descriptor(FILE *fd) {
    LOG_DEBUG("Free file descriptor\n");
    fclose(fd);
}
