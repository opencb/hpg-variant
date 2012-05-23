#include "split.h"

// // Output file descriptors (must be kept between function calls)
// static cp_hashtable *output_files = NULL;
// static FILE *all_variants_file = NULL;
// static FILE *summary_file = NULL;
// // Consequence type counters (for summary, must be kept between function calls)
// static cp_hashtable *summary_count = NULL;
// 
// // Line buffers and their maximum size (one per thread)
// static char **line;
// static char **output_line;
// static int *max_line_size;
// 
// // Output directory (non-accessible directly from CURL callback function)
// static char *output_directory;
// static size_t output_directory_len;
// 
// static int batch_num;
// 
// static list_t *output_list;

int run_split(global_options_data_t *global_options_data, split_options_data_t *options_data) {
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("batches", 1, options_data->max_batches, read_list);
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", options_data->num_threads, MIN(10, options_data->max_batches) * options_data->batch_size, output_list);

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
    
//     // Remove all .txt files in folder
//     ret_code = delete_files_by_extension(global_options_data->output_directory, "txt");
//     if (ret_code != 0) {
//         return ret_code;
//     }
//     output_directory = global_options_data->output_directory;
//     output_directory_len = strlen(output_directory);
//     
//     // Initialize environment for connecting to the web service
//     ret_code = init_http_environment(0);
//     if (ret_code != 0) {
//         return ret_code;
//     }
//     
//     // Initialize collections of file descriptors and summary counters
//     ret_code = initialize_ws_output(options_data->num_threads, global_options_data->output_directory);
//     if (ret_code != 0) {
//         return ret_code;
//     }
//  
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
            FILE *passed_file = NULL, *failed_file = NULL;
            
//             filter_t **filters = NULL;
//             int num_filters = 0;
//             if (options_data->chain != NULL) {
//                 filters = sort_filter_chain(options_data->chain, &num_filters);
//             }
    
            start = omp_get_wtime();

            int i = 0;
            list_item_t* item = NULL;
            
//             if (global_options_data->output_filename != NULL && 
//                 strlen(global_options_data->output_filename) > 0) {
//                 int dirname_len = strlen(global_options_data->output_directory);
//                 int filename_len = strlen(global_options_data->output_filename);
//             
//                 char *passed_filename = (char*) calloc ((dirname_len + filename_len + 1), sizeof(char));
//                 strncat(passed_filename, global_options_data->output_directory, dirname_len);
//                 strncat(passed_filename, global_options_data->output_filename, filename_len);
//                 passed_file = fopen(passed_filename, "w");
//             
//                 char *failed_filename = (char*) calloc ((dirname_len + filename_len + 10), sizeof(char));
//                 strncat(failed_filename, global_options_data->output_directory, dirname_len);
//                 strncat(failed_filename, global_options_data->output_filename, filename_len);
//                 strncat(failed_filename, ".filtered", 9);
//                 failed_file = fopen(failed_filename, "w");
//                 
//                 LOG_DEBUG_F("passed filename = %s\nfailed filename = %s\n", passed_filename, failed_filename);
//                 
//                 free(passed_filename);
//                 free(failed_filename);
//             }
//             
//             // TODO doesn't work (segfault)
// /*            if (!passed_file) {
//                 passed_file = stdout;
//             }
//             if (!failed_file) {
//                 failed_file = stderr;
//             }
//             
//             LOG_DEBUG("File streams created\n");*/
//             
//             // Write file format, header entries and delimiter
//             if (passed_file != NULL) { vcf_write_to_file(file, passed_file); }
//             if (failed_file != NULL) { vcf_write_to_file(file, failed_file); }
// 
//             LOG_DEBUG("VCF header written\n");
//             
            while ((item = list_remove_item(read_list)) != NULL) {
                vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
//                 list_t *input_records = batch;
//                 list_t *passed_records = NULL, *failed_records = NULL;

                if (i % 20 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->length, batch->max_length);
                }
// 
//                 if (filters == NULL) {
//                     passed_records = input_records;
//                 } else {
//                     failed_records = (list_t*) malloc(sizeof(list_t));
//                     list_init("failed_records", 1, INT_MAX, failed_records);
//                     passed_records = run_filter_chain(input_records, failed_records, filters, num_filters);
//                 }
// 
//                 // Write records that passed to a separate file, and query the WS with them as args
//                 if (passed_records->length > 0) {
//                     // Divide the list of passed records in ranges of size defined in config file
//                     int max_chunk_size = options_data->variants_per_request;
//                     int num_chunks;
//                     list_item_t **chunk_starts = create_chunks(passed_records, max_chunk_size, &num_chunks);
//                     
//                     // OpenMP: Launch a thread for each range
//                     #pragma omp parallel for
//                     for (int j = 0; j < num_chunks; j++) {
//                         LOG_DEBUG_F("[%d] WS invocation\n", omp_get_thread_num());
//                         ret_code = invoke_split_ws(url, chunk_starts[j], max_chunk_size);
//                     }
//                     free(chunk_starts);
//                     
//                     LOG_INFO_F("*** %dth web service invocation finished\n", i);
//                     
//                     if (ret_code) {
//                         LOG_FATAL_F("Effect web service error: %s\n", get_last_http_error(ret_code));
//                         break;
//                     }
//                 }
//                 
//                 // Write records that passed and failed to separate files
//                 if (passed_file != NULL && failed_file != NULL) {
//                     if (passed_records != NULL && passed_records->length > 0) {
//                         write_batch(passed_records, passed_file);
//                     }
//                     if (failed_records != NULL && failed_records->length > 0) {
//                         write_batch(failed_records, failed_file);
//                     }
//                 }
//                 
//                 // Free items in both lists (not their internal data)
//                 if (passed_records != input_records) {
//                     LOG_DEBUG_F("[Batch %d] %zu passed records\n", i, passed_records->length);
//                     list_free_deep(passed_records, NULL);
//                 }
//                 if (failed_records) {
//                     LOG_DEBUG_F("[Batch %d] %zu failed records\n", i, failed_records->length);
//                     list_free_deep(failed_records, NULL);
//                 }
                // Free batch and its contents
                vcf_batch_free(item->data_p);
                list_item_free(item);
                
                i++;
            }

            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
// 
//             // Free resources
//             if (passed_file) { fclose(passed_file); }
//             if (failed_file) { fclose(failed_file); }
//             
//             // Free filters
//             for (i = 0; i < num_filters; i++) {
//                 filter_t *filter = filters[i];
//                 filter->free_func(filter);
//             }
//             free(filters);
//             
            // Decrease list writers count
            for (i = 0; i < options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }
        
#pragma omp section
        {
            // Thread which writes the results to all_variants, summary and one file per consequence type
            list_item_t* item = NULL;
            char *line;
            FILE *fd = NULL;
            while ((item = list_remove_item(output_list)) != NULL) {
                line = item->data_p;
//                 
//                 // Write entry in the consequence type file
//                 fd = cp_hashtable_get(output_files, &(item->type));
//                 int ret = fprintf(fd, "%s\n", line);
//                 if (ret < 0) {
//                     LOG_ERROR_F("Error writing to file: '%s'\n", line);
//                 } /*else {
//                     fflush(fd);
//                 }*/
//                 
//                 // Write in all_variants
//                 ret = fprintf(all_variants_file, "%s\n", line);
//                 if (ret < 0) {
//                     LOG_ERROR_F("Error writing to all_variants: '%s'\n", line);
//                 } /*else {
//                     fflush(all_variants_file);
//                 }*/
//                 
                free(line);
                list_item_free(item);
            }
            
        }
    }

//     write_summary_file();
//     write_result_file(global_options_data, options_data);
// 
//     ret_code = free_ws_output(options_data->num_threads);
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
