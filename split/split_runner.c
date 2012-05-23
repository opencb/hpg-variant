#include "split_runner.h"


int run_split(global_options_data_t *global_options_data, split_options_data_t *options_data) {
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("batches", 1, options_data->max_batches, read_list);
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", options_data->num_threads, MIN(10, options_data->max_batches) * options_data->batch_size, output_list);
    cp_hashtable **output_files = (cp_hashtable**) malloc (sizeof(cp_hashtable*));
    initialize_output(output_files);
    
    int ret_code = 0;
    double start, stop, total;
    vcf_file_t *file = vcf_open(global_options_data->vcf_filename);
    
    if (!file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ret_code = create_directory(global_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", global_options_data->output_directory);
    }
    
#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            ret_code = vcf_read_batches(read_list, options_data->batch_size, file, 0);

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) {
                LOG_ERROR_F("Error %d while reading the file %s\n", ret_code, file->filename);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            list_decr_writers(read_list);
        }
        
#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            omp_set_num_threads(options_data->num_threads);
            
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            start = omp_get_wtime();

            int i = 0;
            list_item_t* item = NULL;
            while ((item = list_remove_item(read_list)) != NULL) {
                vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
                list_t *input_records = batch;

                if (i % 50 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->length, batch->max_length);
                }

                // Divide the list of passed records in ranges of size defined in config file
                int max_chunk_size = options_data->variants_per_thread;
                int num_chunks;
                list_item_t **chunk_starts = create_chunks(input_records, max_chunk_size, &num_chunks);
                
                // OpenMP: Launch a thread for each range
                #pragma omp parallel for
                for (int j = 0; j < num_chunks; j++) {
                    LOG_DEBUG_F("[%d] Split invocation\n", omp_get_thread_num());
                    if (options_data->criterion == CHROMOSOME) {
                        ret_code = split_by_chromosome(chunk_starts[j], max_chunk_size, output_list);
                    } else if (options_data->criterion == GENE) {
                        ret_code = 0;
                    }
                }
                if (i % 25 == 0) { LOG_INFO_F("*** %dth split invocation finished\n", i); }
                
                free(chunk_starts);
                vcf_batch_free(item->data_p);
                list_item_free(item);
                
                i++;
            }

            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Decrease list writers count
            for (i = 0; i < options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }
        
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d writes the output\n", omp_get_thread_num());
            
            char *split_filename;
            FILE *split_fd;
    
            // Create file streams for results)
            int dirname_len = strlen(global_options_data->output_directory);
            int filename_len;
            if (global_options_data->output_filename == NULL || strlen(global_options_data->output_filename) == 0) {
                filename_len = strlen("split-tool-output");
            
                split_filename = (char*) calloc ((dirname_len + filename_len + 1), sizeof(char));
                strncat(split_filename, global_options_data->output_directory, dirname_len);
                strncat(split_filename, "split-tool-output", filename_len);
            } else {
                filename_len = strlen(global_options_data->output_filename);
            
                split_filename = (char*) calloc ((dirname_len + filename_len + 1), sizeof(char));
                strncat(split_filename, global_options_data->output_directory, dirname_len);
                strncat(split_filename, global_options_data->output_filename, filename_len);
            }
            
            LOG_DEBUG_F("split output filename = %s\n", split_filename);
            split_fd = fopen(split_filename, "w");
            free(split_filename);
            LOG_DEBUG("File streams created\n");
            
            list_item_t* item = NULL;
            FILE *fd = NULL;
            while ((item = list_remove_item(output_list)) != NULL) {
                split = item->data_p;
                
                // TODO If its the first line to write into the file, include the header
                
                // TODO write line into the file
            }
            
            // Close file
            if (split_fd != NULL) { fclose(split_fd); }
        }
    }

//     write_summary_file();
//     write_result_file(global_options_data, options_data);

    ret_code = free_output(*output_files);
    free(read_list);
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
