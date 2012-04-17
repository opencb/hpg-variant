#include "vcf_tools_runner.h"

int execute_vcf_tools(global_options_data_t *global_options_data, vcf_tools_options_data_t *options_data) {
    
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("batches", 1, options_data->max_batches, read_list);

    int ret_code;
    double start, stop, total;
    vcf_file_t *file = vcf_open(global_options_data->vcf_filename);
    
#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            dprintf("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            ret_code = vcf_read_batches(read_list, options_data->batch_size, file, 1);

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) dprintf("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code);

            bprintf("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            bprintf("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            list_decr_writers(read_list);
        }
        
#pragma omp section
        {
            FILE *passed_file, *failed_file;
            
            filter_t **filters = NULL;
            int num_filters = 0;
            if (options_data->chain != NULL) {
                filters = sort_filter_chain(options_data->chain, &num_filters);
            }
    
            start = omp_get_wtime();

            int i = 0;
            list_item_t* item = NULL;
            
            if (global_options_data->output_filename != NULL && 
                strlen(global_options_data->output_filename) > 0) {
                passed_file = fopen(global_options_data->output_filename, "w");
                int filename_len = strlen(global_options_data->output_filename);
                char *failed_filename = (char*) malloc ((filename_len+10) * sizeof(char));
                strncpy(failed_filename, global_options_data->output_filename, filename_len);
                strncpy(failed_filename + filename_len, ".filtered", 9);
                failed_filename[filename_len+10] = '\0';
                failed_file = fopen(failed_filename, "w");
            }
            
            // TODO doesn't work (segfault)
            if (!passed_file) {
                passed_file = stdout;
            }
            if (!failed_file) {
                failed_file = stderr;
            }
            
            dprintf("File streams created\n");
            
            // Write file format, header entries and delimiter
            vcf_write_to_file(file, passed_file);
            vcf_write_to_file(file, failed_file);

            dprintf("VCF header written created\n");
            
            while ((item = list_remove_item(read_list)) != NULL) {
                vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
                list_t *input_records = batch;
                list_t *passed_records, *failed_records;

                if (i % 200 == 0) {
                    dprintf("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->length, batch->max_length);
                }

                if (filters != NULL) {
                    failed_records = (list_t*) malloc(sizeof(list_t));
                    list_init("failed_records", 1, INT_MAX, failed_records);
                    passed_records = run_filter_chain(input_records, failed_records, filters, num_filters);
                }

                
//                 if(is_split()) {
//                         split(passed_records, batch);
//                         
//                 }
//                 if(is_sort()) {
//                         sort(passed_records);
//                         
//                 }
//                 if(is_stats()) {
//                         stats(passed_records);
//                 }
//                 if(is_merge()) {
//                         merge(passed_records, ...);
//                 }
//                 if(is_validate()) {
//                         validate(passed_records);
//                 }
                

                // Write records that passed and failed to 2 new separated files
                if (passed_records != NULL && passed_records->length > 0) {
                    dprintf("[batch %d] %zu passed records\n", i, passed_records->length);
                    write_batch(passed_records, passed_file);
                }
                
                if (failed_records != NULL && failed_records->length > 0) {
                    dprintf("[batch %d] %zu failed records\n", i, failed_records->length);
                    write_batch(failed_records, failed_file);
                }
                
                vcf_batch_free(item->data_p);
                list_item_free(item);
                
                if (passed_records) free(passed_records);
                if (failed_records) free(failed_records);
                
                i++;
            }

            stop = omp_get_wtime();

            total = stop - start;

            bprintf("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            bprintf("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Free resources
            fclose(passed_file);
            fclose(failed_file);
            
            free(filters);
        }
    }
    
    vcf_close(file);
    
    free(read_list);
    
    return 0;
}
