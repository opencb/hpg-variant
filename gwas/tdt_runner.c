#include "tdt_runner.h"

// int permute = 0;

int run_tdt_test(global_options_data_t* global_options_data, gwas_options_data_t* options_data) {
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("batches", 1, options_data->max_batches, read_list);
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", options_data->num_threads, options_data->max_batches * options_data->batch_size, output_list);

    int ret_code = 0;
    double start, stop, total;
    vcf_file_t *file = vcf_open(global_options_data->vcf_filename);
    ped_file_t *ped_file = ped_open(global_options_data->ped_filename);
    size_t output_directory_len = strlen(global_options_data->output_directory);
    
    LOG_INFO("About to read PED file...\n");
    // Read PED file before doing any proccessing
    ret_code = ped_read(ped_file);
    if (ret_code != 0) {
        LOG_FATAL_F("Can't read PED file: %s\n", ped_file->filename);
    }
    
    // Try to create the directory where the output files will be stored
    ret_code = create_directory(global_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", global_options_data->output_directory);
    }
    
    LOG_INFO("About to perform TDT test...\n");

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

            if (ret_code) {
                LOG_FATAL_F("Error %d while reading the file %s\n", ret_code, file->filename);
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
            cp_hashtable *sample_ids = NULL;
            
            // Create chain of filters for the VCF file
            filter_t **filters = NULL;
            int num_filters = 0;
            if (options_data->chain != NULL) {
                filters = sort_filter_chain(options_data->chain, &num_filters);
            }
    
            start = omp_get_wtime();

            if (global_options_data->output_filename != NULL && 
                strlen(global_options_data->output_filename) > 0) {
                int dirname_len = strlen(global_options_data->output_directory);
                int filename_len = strlen(global_options_data->output_filename);
            
                char *passed_filename = (char*) calloc (dirname_len + filename_len + 11, sizeof(char));
                sprintf(passed_filename, "%s/%s.filtered", global_options_data->output_directory, global_options_data->output_filename);
                passed_file = fopen(passed_filename, "w");
            
                char *failed_filename = (char*) calloc (dirname_len + filename_len + 2, sizeof(char));
                sprintf(failed_filename, "%s/%s", global_options_data->output_directory, global_options_data->output_filename);
                failed_file = fopen(failed_filename, "w");
                
                LOG_DEBUG_F("passed filename = %s\nfailed filename = %s\n", passed_filename, failed_filename);
                
                free(passed_filename);
                free(failed_filename);
            }
            
            int i = 0;
            list_item_t *item = NULL;
            while ((item = list_remove_item(read_list)) != NULL) {
                if (i == 0) {
                    // Create map to associate the position of individuals in the list of samples defined in the VCF file
                    sample_ids = associate_samples_and_positions(file);
                    
                    // Write file format, header entries and delimiter
                    if (passed_file != NULL) { vcf_write_to_file(file, passed_file); }
                    if (failed_file != NULL) { vcf_write_to_file(file, failed_file); }
                    
                    LOG_DEBUG("VCF header written\n");
                }
                
                vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
                list_t *input_records = batch;
                list_t *passed_records = NULL, *failed_records = NULL;

                if (i % 20 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->length, batch->max_length);
                }

                if (filters == NULL) {
                    passed_records = input_records;
                } else {
                    failed_records = (list_t*) malloc(sizeof(list_t));
                    list_init("failed_records", 1, INT_MAX, failed_records);
                    passed_records = run_filter_chain(input_records, failed_records, filters, num_filters);
                }

                // Launch TDT test over records that passed the filters
                if (passed_records->length > 0) {
                    // Divide the list of passed records in ranges of size defined in config file
                    int max_chunk_size = options_data->variants_per_request;
                    int num_chunks;
                    list_item_t **chunk_starts = create_chunks(passed_records, max_chunk_size, &num_chunks);
                    
                    // OpenMP: Launch a thread for each range
                    #pragma omp parallel for
                    for (int j = 0; j < num_chunks; j++) {
                        LOG_DEBUG_F("[%d] Test execution\n", omp_get_thread_num());
                        ret_code = tdt_test(ped_file, chunk_starts[j], max_chunk_size, sample_ids, output_list);
                    }
                    free(chunk_starts);
                    
                    if (i % 10 == 0) { 
                        LOG_INFO_F("*** %dth TDT execution finished\n", i);
                    }
                    
                    if (ret_code) {
//                         LOG_FATAL_F("TDT error: %s\n", get_last_http_error(ret_code));
                        break;
                    }
                }
                
                // Write records that passed and failed to separate files
                if (passed_file != NULL && failed_file != NULL) {
                    if (passed_records != NULL && passed_records->length > 0) {
                        write_batch(passed_records, passed_file);
                    }
                    if (failed_records != NULL && failed_records->length > 0) {
                        write_batch(failed_records, failed_file);
                    }
                }
                
                // Free items in both lists (not their internal data)
                if (passed_records != input_records) {
                    LOG_DEBUG_F("[Batch %d] %zu passed records\n", i, passed_records->length);
                    list_free_deep(passed_records, NULL);
                }
                if (failed_records) {
                    LOG_DEBUG_F("[Batch %d] %zu failed records\n", i, failed_records->length);
                    list_free_deep(failed_records, NULL);
                }
                // Free batch and its contents
                vcf_batch_free(item->data_p);
                list_item_free(item);
                
                i++;
            }

            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Free resources
            if (sample_ids) { cp_hashtable_destroy(sample_ids); }
            
            // Free filters
            for (i = 0; i < num_filters; i++) {
                filter_t *filter = filters[i];
                filter->free_func(filter);
            }
            free(filters);
            
            // Decrease list writers count
            for (i = 0; i < options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }

#pragma omp section
        {
            // Thread which writes the results to the output file
            FILE *fd = NULL;    // TODO check if output file is defined
            char *path = NULL, *filename = NULL;
            size_t filename_len = 0;
            
            // Set whole path to the output file
            if (global_options_data->output_filename != NULL && 
                strlen(global_options_data->output_filename) > 0) {
                filename_len = strlen(global_options_data->output_filename);
                filename = global_options_data->output_filename;
            } else {
                filename_len = strlen("hpg-variant.tdt");
                filename = strdup("hpg-variant.tdt");
            }
            path = (char*) calloc ((output_directory_len + filename_len + 2), sizeof(char));
            sprintf(path, "%s/%s", global_options_data->output_directory, filename);
            fd = fopen(path, "w");
            
            LOG_INFO_F("TDT output filename = %s\n", path);
            free(filename);
            free(path);
            
            // Write data: header + one line per variant
            list_item_t* item = NULL;
            tdt_result_t *result;
            fprintf(fd, "CHR           BP       A1      A2         T       U           OR           CHISQ         P-VALUE\n");
            while ((item = list_remove_item(output_list)) != NULL) {
                result = item->data_p;
                
                fprintf(fd, "%s\t%8ld\t%s\t%s\t%3d\t%3d\t%6f\t%6f\t%6f\n",
                       result->chromosome, result->position, result->reference, result->alternate, 
                       result->t1, result->t2, result->odds_ratio, result->chi_square, result->p_value);
                
                tdt_result_free(result);
                list_item_free(item);
            }
            
            fclose(fd);
        }
    }
    
    free(read_list);
    free(output_list);
    vcf_close(file);
    // TODO delete conflicts among frees
    ped_close(ped_file, 0);
        
    return ret_code;
}


cp_hashtable* associate_samples_and_positions(vcf_file_t* file) {
    LOG_DEBUG_F("** %zu sample names read\n", file->samples_names->length);
    list_t *sample_names = file->samples_names;
    cp_hashtable *sample_ids = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                                             sample_names->length * 2,
                                                             cp_hash_string,
                                                             (cp_compare_fn) strcasecmp,
                                                             NULL,
                                                             NULL,
                                                             NULL,
                                                             (cp_destructor_fn) free
                                                            );
    
    list_item_t *sample_item = sample_names->first_p;
    int *index;
    char *name;
    for (int i = 0; i < sample_names->length && sample_item != NULL; i++) {
        name = sample_item->data_p;
        index = (int*) malloc (sizeof(int)); *index = i;
        cp_hashtable_put(sample_ids, (char*) sample_item->data_p, index);
        
        sample_item = sample_item->next_p;
    }
    
//     char **keys = (char**) cp_hashtable_get_keys(sample_names);
//     int num_keys = cp_hashtable_count(sample_names);
//     for (int i = 0; i < num_keys; i++) {
//         printf("%s\t%d\n", keys[i], *((int*) cp_hashtable_get(sample_ids, keys[i])));
//     }
    
    return sample_ids;
}
