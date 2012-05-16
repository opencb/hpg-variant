#include "stats.h"

int run_stats(global_options_data_t *global_options_data, stats_options_data_t *options_data) {
    
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("batches", 1, options_data->max_batches, read_list);

    int ret_code;
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

            ret_code = vcf_read_batches(read_list, options_data->batch_size, file, 1);

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            list_decr_writers(read_list);
        }
        
#pragma omp section
        {
            char *stats_filename;
            FILE *stats_file;
    
            if (global_options_data->output_filename == NULL || 
                strlen(global_options_data->output_filename) == 0) {
                int dirname_len = strlen(global_options_data->output_directory);
                int filename_len = strlen("stats-tool-output.vcf");
            
                stats_filename = (char*) calloc ((dirname_len + filename_len + 1), sizeof(char));
                strncat(stats_filename, global_options_data->output_directory, dirname_len);
                strncat(stats_filename, "stats-tool-output.vcf", filename_len);
            } else {
                int dirname_len = strlen(global_options_data->output_directory);
                int filename_len = strlen(global_options_data->output_filename);
            
                stats_filename = (char*) calloc ((dirname_len + filename_len + 1), sizeof(char));
                strncat(stats_filename, global_options_data->output_directory, dirname_len);
                strncat(stats_filename, global_options_data->output_filename, filename_len);
            }
            
            LOG_DEBUG_F("stats filename = %s\n", stats_filename);
            stats_file = fopen(stats_filename, "w");
            free(stats_filename);
            LOG_DEBUG("File streams created\n");
            
            start = omp_get_wtime();
            
            int i = 0;
            list_item_t* item = NULL;
            while ((item = list_remove_item(read_list)) != NULL) {
                vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
                list_t *input_records = batch;

                if (i % 200 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->length, batch->max_length);
                }

                
//                 // Write records that passed and failed to 2 new separated files
//                 if (passed_records != NULL && passed_records->length > 0) {
//                     LOG_DEBUG_F("[batch %d] %zu passed records\n", i, passed_records->length);
//                     write_batch(passed_records, stats_file);
//                 }
//                 
//                 if (failed_records != NULL && failed_records->length > 0) {
//                     LOG_DEBUG_F("[batch %d] %zu failed records\n", i, failed_records->length);
//                     write_batch(failed_records, failed_file);
//                 }
                
                vcf_batch_free(item->data_p);
                list_item_free(item);
                
                i++;
            }

            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Free resources
            if (stats_file != NULL) { fclose(stats_file); }
        }
    }
    
    vcf_close(file);
    
    free(read_list);
    
    return 0;
}

