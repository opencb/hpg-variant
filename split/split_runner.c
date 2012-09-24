#include "split_runner.h"


int run_split(shared_options_data_t *shared_options_data, split_options_data_t *options_data) {
//     list_t *read_list = (list_t*) malloc(sizeof(list_t));
//     list_init("batches", 1, shared_options_data->max_batches, read_list);
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, MIN(10, shared_options_data->max_batches) * shared_options_data->batch_lines, output_list);
//     cp_hashtable **output_files = (cp_hashtable**) malloc (sizeof(cp_hashtable*));
//     initialize_output(output_files);
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
                ret_code = vcf_parse_batches_in_bytes(shared_options_data->batch_bytes, file, 1);
            } else if (shared_options_data->batch_lines > 0) {
                ret_code = vcf_parse_batches(shared_options_data->batch_lines, file, 1);
            }

//             ret_code = vcf_parse_batches(read_list, shared_options_data->batch_lines, file, 1);

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) {
                LOG_ERROR_F("Error %d while reading the file %s\n", ret_code, file->filename);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            notify_end_reading(file);
//             list_decr_writers(read_list);
        }
        
#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            omp_set_num_threads(shared_options_data->num_threads);
            
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            start = omp_get_wtime();

            int i = 0;
//             list_item_t* item = NULL;
//             while ((item = list_remove_item(read_list)) != NULL) {
            vcf_batch_t *batch = NULL;
            while ((batch = fetch_vcf_batch(file)) != NULL) {
//                 vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
                array_list_t *input_records = batch->records;

                if (i % 100 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->records->size, batch->records->capacity);
                }

                // Divide the list of passed records in ranges of size defined in config file
                int num_chunks;
                int *chunk_sizes;
                int *chunk_starts = create_chunks(input_records->size, shared_options_data->entries_per_thread, &num_chunks, &chunk_sizes);
                
                // OpenMP: Launch a thread for each range
                #pragma omp parallel for
                for (int j = 0; j < num_chunks; j++) {
                    LOG_DEBUG_F("[%d] Split invocation\n", omp_get_thread_num());
                    if (options_data->criterion == CHROMOSOME) {
                        ret_code = split_by_chromosome((vcf_record_t**) (input_records->items + chunk_starts[j]), 
                                                       chunk_sizes[j],
                                                       output_list);
                    } else if (options_data->criterion == GENE) {
                        ret_code = 0;
                    }
                }
//                 if (i % 50 == 0) { LOG_INFO_F("*** %dth split invocation finished\n", i); }
                
                free(chunk_starts);
                // Can't free batch contents because they will be written to file by another thread
//                 free(item->data_p);
//                 list_item_free(item);
                
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
            int dirname_len = strlen(shared_options_data->output_directory);
            
            list_item_t* item = NULL;
            split_result_t *split;
            FILE *split_fd = NULL;
            char split_filename[1024];
            char input_filename[256];
            get_filename_from_path(shared_options_data->vcf_filename, input_filename);
            
            while ((item = list_remove_item(output_list)) != NULL) {
                split = item->data_p;
                
                memset(split_filename, 0, 1024 * sizeof(char));
                sprintf(split_filename, "%s/%s_%s", shared_options_data->output_directory, split->split_name, input_filename);
//                 sprintf(split_filename, "%s/%s.vcf", shared_options_data->output_directory, split->split_name, shared_options_data->vcf_filename);
                
//                 printf("Split filename = '%s'\n", split_filename);
                
                split_fd = cp_hashtable_get(output_files, split->split_name);
                if (!split_fd) {
                    // TODO If its the first line to write into the file, create file and include the header
                    split_fd = fopen(split_filename, "w");
                    cp_hashtable_put(output_files, split->split_name, split_fd);
                    
                    write_vcf_file(file, split_fd);
                }
                
                // TODO write line into the file
                write_vcf_record(split->record, split_fd);
                vcf_record_free(split->record);
            }
            
            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        }
    }

    free_output(output_files);
//     free(read_list);
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
