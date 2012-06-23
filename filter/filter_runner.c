#include "filter.h"

int run_filter(shared_options_data_t *shared_options_data, filter_options_data_t *options_data) {
    
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("batches", 1, shared_options_data->max_batches, read_list);

    int ret_code;
    double start, stop, total;
    vcf_file_t *file = vcf_open(shared_options_data->vcf_filename);
    
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

            ret_code = vcf_read_batches(read_list, shared_options_data->batch_size, file, 1);

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            list_decr_writers(read_list);
        }
        
#pragma omp section
        {
            char *passed_filename, *failed_filename;
            FILE *passed_file, *failed_file;
            
            filter_t **filters = NULL;
            int num_filters = 0;
            if (options_data->chain != NULL) {
                filters = sort_filter_chain(options_data->chain, &num_filters);
            }
    
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
            
            LOG_INFO_F("passed filename = %s\nfailed filename = %s\n", passed_filename, failed_filename);
            passed_file = fopen(passed_filename, "w");
            failed_file = fopen(failed_filename, "w");
            free(passed_filename);
            free(failed_filename);
            LOG_DEBUG("File streams created\n");
            
            start = omp_get_wtime();

            int i = 0;
            list_item_t* item = NULL;
            while ((item = list_remove_item(read_list)) != NULL) {
                if (i == 0) {
                    // Write file format, header entries and delimiter
                    vcf_write_to_file(file, passed_file);
                    vcf_write_to_file(file, failed_file);

                    LOG_DEBUG("VCF header written created\n");
                }
                
                vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
                list_t *input_records = batch;
                list_t *passed_records, *failed_records;

                if (i % 100 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->length, batch->max_length);
                }

                if (filters != NULL) {
                    failed_records = (list_t*) malloc(sizeof(list_t));
                    list_init("failed_records", 1, INT_MAX, failed_records);
                    passed_records = run_filter_chain(input_records, failed_records, filters, num_filters);
                } else {
                    passed_records = input_records;
                }

                // Write records that passed and failed to 2 new separated files
                if (passed_records != NULL && passed_records->length > 0) {
                    LOG_DEBUG_F("[batch %d] %zu passed records\n", i, passed_records->length);
                    write_batch(passed_records, passed_file);
                }
                
                if (failed_records != NULL && failed_records->length > 0) {
                    LOG_DEBUG_F("[batch %d] %zu failed records\n", i, failed_records->length);
                    write_batch(failed_records, failed_file);
                }
                
                vcf_batch_free(item->data_p);
                list_item_free(item);
                
                if (passed_records && passed_records != input_records) { free(passed_records); }
                if (failed_records) { free(failed_records); }
                
                i++;
            }

            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Free resources
            if (passed_file) { fclose(passed_file); }
            if (failed_file) { fclose(failed_file); }
            
            free(filters);
        }
    }
    
    vcf_close(file);
    
    free(read_list);
    
    return 0;
}
