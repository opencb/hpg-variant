#include "filter.h"

int run_filter(shared_options_data_t *shared_options_data, filter_options_data_t *options_data) {
    int ret_code;
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

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            notify_end_reading(file);
        }
        
#pragma omp section
        {
            filter_t **filters = NULL;
            int num_filters = 0;
            if (options_data->chain != NULL) {
                filters = sort_filter_chain(options_data->chain, &num_filters);
            }
    
            FILE *passed_file = NULL, *failed_file = NULL;
            char *passed_filename, *failed_filename;
//             get_output_files(shared_options_data, &passed_file, &failed_file);
            
            if (shared_options_data->output_filename == NULL || strlen(shared_options_data->output_filename) == 0) {
                int dirname_len = strlen(shared_options_data->output_directory);
                int filename_len = strlen("filter-tool-output");
                passed_filename = (char*) calloc ((dirname_len + filename_len + 15), sizeof(char));
                sprintf(passed_filename, "%s/%s.vcf.filtered", shared_options_data->output_directory, "filter-tool-output");
                failed_filename = (char*) calloc ((dirname_len + filename_len + 15), sizeof(char));
                sprintf(failed_filename, "%s/%s.vcf.rejected", shared_options_data->output_directory, "filter-tool-output");
                
            } else {
                int dirname_len = strlen(shared_options_data->output_directory);
                int filename_len = strlen(shared_options_data->output_filename);
                passed_filename = (char*) calloc ((dirname_len + filename_len + 11), sizeof(char));
                sprintf(passed_filename, "%s/%s.filtered", shared_options_data->output_directory, shared_options_data->output_filename);            
                failed_filename = (char*) calloc ((dirname_len + filename_len + 11), sizeof(char));
                sprintf(failed_filename, "%s/%s.rejected", shared_options_data->output_directory, shared_options_data->output_filename);
            }
            
            LOG_DEBUG_F("passed filename = %s\nfailed filename = %s\n", passed_filename, failed_filename);
            passed_file = fopen(passed_filename, "w");
            if (options_data->save_rejected) {
                failed_file = fopen(failed_filename, "w");
            }
            free(passed_filename);
            free(failed_filename);
            LOG_DEBUG("File streams created\n");
            
            start = omp_get_wtime();

            int i = 0;
            vcf_batch_t *batch = NULL;
            while ((batch = fetch_vcf_batch(file)) != NULL) {
                if (i == 0) {
                    // Write file format, header entries and delimiter
                    write_vcf_header(file, passed_file);
                    if (options_data->save_rejected) {
                        write_vcf_header(file, failed_file);
                    }

                    LOG_DEBUG("VCF header written created\n");
                }
                
                array_list_t *input_records = batch->records;
                array_list_t *passed_records, *failed_records;

                if (i % 100 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->records->size, batch->records->capacity);
                }

                if (filters == NULL) {
                    passed_records = input_records;
                } else {
                    failed_records = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
                    passed_records = run_filter_chain(input_records, failed_records, filters, num_filters);
                }

                // Write records that passed and failed to 2 new separated files
                if (passed_records != NULL && passed_records->size > 0) {
                    LOG_DEBUG_F("[batch %d] %zu passed records\n", i, passed_records->size);
                #pragma omp critical 
                    {
                        for (int r = 0; r < passed_records->size; r++) {
                            write_vcf_record(passed_records->items[r], passed_file);
                        }
//                         write_batch(passed_records, passed_file);
                    }
                }
                
                if (options_data->save_rejected && failed_records != NULL && failed_records->size > 0) {
                    LOG_DEBUG_F("[batch %d] %zu failed records\n", i, failed_records->size);
                #pragma omp critical 
                    {
                        for (int r = 0; r < failed_records->size; r++) {
                            write_vcf_record(failed_records->items[r], failed_file);
                        }
//                         write_batch(failed_records, failed_file);
                    }
                }
                
                // Free batch and its contents
                vcf_batch_free(batch);
                
                // Free items in both lists (not their internal data)
                if (passed_records != input_records) {
                    array_list_free(passed_records, NULL);
                }
                if (failed_records) {
                    array_list_free(failed_records, NULL);
                }
                
                i++;
            }

            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Free resources
            if (passed_file) { fclose(passed_file); }
            if (failed_file) { fclose(failed_file); }
            
            free_filters(filters, num_filters);
        }
    }
    
    vcf_close(file);
    
    return 0;
}
