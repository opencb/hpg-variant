#include "assoc_runner.h"

int run_association_test(shared_options_data_t* shared_options_data, gwas_options_data_t* options_data) {
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("text", 1, shared_options_data->max_batches, read_list);
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, INT_MAX, output_list);

    int ret_code = 0;
    vcf_file_t *file = vcf_open(shared_options_data->vcf_filename);
    if (!file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ped_file_t *ped_file = ped_open(shared_options_data->ped_filename);
    if (!ped_file) {
        LOG_FATAL("PED file does not exist!\n");
    }
    
    size_t output_directory_len = strlen(shared_options_data->output_directory);
    
    LOG_INFO("About to read PED file...\n");
    // Read PED file before doing any proccessing
    ret_code = ped_read(ped_file);
    if (ret_code != 0) {
        LOG_FATAL_F("Can't read PED file: %s\n", ped_file->filename);
    }
               
    // Try to create the directory where the output files will be stored
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }
    
    LOG_INFO("About to perform basic association test...\n");

#pragma omp parallel sections num_threads(3) private (ret_code)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            double start = omp_get_wtime();

            ret_code = 0;
            if (shared_options_data->batch_bytes > 0) {
                ret_code = vcf_read_batches_in_bytes(read_list, shared_options_data->batch_bytes, file);
            } else if (shared_options_data->batch_lines > 0) {
                ret_code = vcf_read_batches(read_list, shared_options_data->batch_lines, file);
            }

            double stop = omp_get_wtime();
            double total = stop - start;

            if (ret_code) {
                LOG_FATAL_F("Error %d while reading the file %s\n", ret_code, file->filename);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            list_decr_writers(read_list);
        }

#pragma omp section
        {
            // Enable nested parallelism
            omp_set_nested(1);
            
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            volatile int initialization_done = 0;
            cp_hashtable *sample_ids = NULL;
            
            // Create chain of filters for the VCF file
            filter_t **filters = NULL;
            int num_filters = 0;
            if (shared_options_data->chain != NULL) {
                filters = sort_filter_chain(shared_options_data->chain, &num_filters);
            }
            FILE *passed_file = NULL, *failed_file = NULL;
            get_output_files(shared_options_data, &passed_file, &failed_file);
    
            double start = omp_get_wtime();
            
            double *factorial_logarithms = NULL;
            
// #pragma omp parallel num_threads(shared_options_data->num_threads) shared(initialization_done, sample_ids, factorial_logarithms, filters)
#pragma omp parallel num_threads(shared_options_data->num_threads) shared(initialization_done, factorial_logarithms, filters)
            {
            family_t **families = (family_t**) cp_hashtable_get_values(ped_file->families);
            int num_families = get_num_families(ped_file);
            
            list_t *vcf_batches_list = (list_t*) malloc(sizeof(list_t));
            list_init("batches", 1, shared_options_data->max_batches, vcf_batches_list);
             
            int i = 0;
            list_item_t *item = NULL;
            while ((item = list_remove_item(read_list)) != NULL) {
                char *text_begin = item->data_p;
                char *text_end = text_begin + strlen(text_begin);
                
                assert(text_end != NULL);
                assert(vcf_batches_list != NULL);
                
                vcf_reader_status *status = vcf_reader_status_new(shared_options_data->batch_lines, 1, 1);
                execute_vcf_ragel_machine(text_begin, text_end, vcf_batches_list, shared_options_data->batch_lines, file, status);
                
                // Initialize structures needed for TDT and write headers of output files
                if (!initialization_done) {
                    sample_ids = associate_samples_and_positions(file);
# pragma omp critical
                {
                    // Guarantee that just one thread performs this operation
                    if (!initialization_done) {
                        // Create map to associate the position of individuals in the list of samples defined in the VCF file
//                         sample_ids = associate_samples_and_positions(file);
                        
                        // Write file format, header entries and delimiter
                        if (passed_file != NULL) { vcf_write_to_file(file, passed_file); }
                        if (failed_file != NULL) { vcf_write_to_file(file, failed_file); }
                        
                        LOG_DEBUG("VCF header written\n");
                        
                        if (options_data->task == FISHER) {
                            factorial_logarithms = init_logarithm_array(file->num_samples * 10);
                        }
                        
                        initialization_done = 1;
                    }
                }
                }
                
                list_item_t *batch_item = list_remove_item(vcf_batches_list);
                vcf_batch_t *batch = batch_item->data_p;
                
                array_list_t *input_records = batch->records;
                array_list_t *passed_records = NULL, *failed_records = NULL;

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

                // Launch TDT test over records that passed the filters
                int num_variants = MIN(shared_options_data->batch_lines, passed_records->size);
                if (passed_records->size > 0) {
//                     LOG_DEBUG_F("[%d] Test execution\n", omp_get_thread_num());
                    // TODO is this if neccessary? factorial_logarithms will be null in chi-square
//                     if (options_data->task == ASSOCIATION_BASIC) {
//                         assoc_test(options_data->task, (vcf_record_t**) passed_records->items, num_variants, 
//                                    families, num_families, sample_ids, NULL, output_list);
//                     } else if (options_data->task == FISHER) {
                        assoc_test(options_data->task, (vcf_record_t**) passed_records->items, num_variants, 
                                   families, num_families, sample_ids, factorial_logarithms, output_list);
//                     }
                }
                
                // Write records that passed and failed to separate files
                if (passed_file != NULL && failed_file != NULL) {
                    if (passed_records != NULL && passed_records->size > 0) {
                #pragma omp critical 
                    {
                        for (int r = 0; r < passed_records->size; r++) {
                            write_record(passed_records->items[r], passed_file);
                        }
//                         write_batch(passed_records, passed_file);
                    }
                    }
                    if (failed_records != NULL && failed_records->size > 0) {
                #pragma omp critical 
                    {
                        for (int r = 0; r < passed_records->size; r++) {
                            write_record(failed_records->items[r], failed_file);
                        }
//                         write_batch(failed_records, failed_file);
                    }
                    }
                }
                
                // Free items in both lists (not their internal data)
                if (passed_records != input_records) {
//                     LOG_DEBUG_F("[Batch %d] %zu passed records\n", i, passed_records->size);
                    array_list_free(passed_records, NULL);
                }
                if (failed_records) {
//                     LOG_DEBUG_F("[Batch %d] %zu failed records\n", i, failed_records->size);
                    array_list_free(failed_records, NULL);
                }
                
                // Free batch and its contents
                vcf_reader_status_free(status);
                vcf_batch_free(batch);
                list_item_free(batch_item);
//                 free(item->data_p);
                list_item_free(item);
                
                i++;
            }  
            }

            double stop = omp_get_wtime();
            double total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Free resources
            if (sample_ids) { cp_hashtable_destroy(sample_ids); }
            
            // Free filters
            for (int i = 0; i < num_filters; i++) {
                filter_t *filter = filters[i];
                filter->free_func(filter);
            }
            free(filters);
            
            // Decrease list writers count
            for (int i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }

#pragma omp section
        {
            // Thread which writes the results to the output file
            FILE *fd = NULL;
            char *path = NULL, *filename = NULL;
            size_t filename_len = 0;
            
            // Set whole path to the output file
            if (options_data->task == ASSOCIATION_BASIC) {
                filename_len = strlen("hpg-variant.assoc");
                filename = strdup("hpg-variant.assoc");
            } else if (options_data->task == FISHER) {
                filename_len = strlen("hpg-variant.fisher");
                filename = strdup("hpg-variant.fisher");
            }
            path = (char*) calloc ((output_directory_len + filename_len + 2), sizeof(char));
            sprintf(path, "%s/%s", shared_options_data->output_directory, filename);
            fd = fopen(path, "w");
            
            LOG_INFO_F("Association test output filename = %s\n", path);
            free(filename);
            free(path);
            
            double start = omp_get_wtime();
            
            // Write data: header + one line per variant
            list_item_t* item = NULL;
            
            if (options_data->task == ASSOCIATION_BASIC) {
                assoc_basic_result_t *result;
                fprintf(fd, "#CHR          BP       A1      C_A1    C_U1         F_A1            F_U1       A2      C_A2    C_U2         F_A2            F_U2              OR           CHISQ         P-VALUE\n");
                while ((item = list_remove_item(output_list)) != NULL) {
                    result = item->data_p;
                    
                    fprintf(fd, "%s\t%8ld\t%s\t%3d\t%3d\t%6f\t%6f\t%s\t%3d\t%3d\t%6f\t%6f\t%6f\t%6f\t%6f\n",
                            result->chromosome, result->position, 
                            result->reference, result->affected1, result->unaffected1, 
                            (double) result->affected1 / (result->affected1 + result->affected2), 
                            (double) result->unaffected1 / (result->unaffected1 + result->unaffected2), 
                            result->alternate, result->affected2, result->unaffected2, 
                            (double) result->affected2 / (result->affected1 + result->affected2), 
                            (double) result->unaffected2 / (result->unaffected1 + result->unaffected2), 
                            result->odds_ratio, result->chi_square, result->p_value);
                    
                    assoc_basic_result_free(result);
                    list_item_free(item);
                }
            } else if (options_data->task == FISHER) {
                assoc_fisher_result_t *result;
                fprintf(fd, "#CHR          BP       A1      C_A1    C_U1         F_A1            F_U1       A2      C_A2    C_U2         F_A2            F_U2              OR         P-VALUE\n");
                while ((item = list_remove_item(output_list)) != NULL) {
                    result = item->data_p;
                    
                    fprintf(fd, "%s\t%8ld\t%s\t%3d\t%3d\t%6f\t%6f\t%s\t%3d\t%3d\t%6f\t%6f\t%6f\t%6f\n",//\t%f\n",
                            result->chromosome, result->position, 
                            result->reference, result->affected1, result->unaffected1, 
                            (double) result->affected1 / (result->affected1 + result->affected2), 
                            (double) result->unaffected1 / (result->unaffected1 + result->unaffected2), 
                            result->alternate, result->affected2, result->unaffected2, 
                            (double) result->affected2 / (result->affected1 + result->affected2), 
                            (double) result->unaffected2 / (result->unaffected1 + result->unaffected2), 
                            result->odds_ratio, result->p_value);
                    
                    assoc_fisher_result_free(result);
                    list_item_free(item);
                }
            }
            
            fclose(fd);
            
            double stop = omp_get_wtime();
            double total = stop - start;

            LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
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
    LOG_DEBUG_F("** %zu sample names read\n", file->samples_names->size);
    array_list_t *sample_names = file->samples_names;
    cp_hashtable *sample_ids = cp_hashtable_create(sample_names->size * 2,
                                                   cp_hash_string,
                                                   (cp_compare_fn) strcasecmp
                                                  );
    
    int *index;
    char *name;
    for (int i = 0; i < sample_names->size; i++) {
        name = sample_names->items[i];
        index = (int*) malloc (sizeof(int)); *index = i;
        cp_hashtable_put(sample_ids, name, index);
    }
//     char **keys = (char**) cp_hashtable_get_keys(sample_names);
//     int num_keys = cp_hashtable_count(sample_names);
//     for (int i = 0; i < num_keys; i++) {
//         printf("%s\t%d\n", keys[i], *((int*) cp_hashtable_get(sample_ids, keys[i])));
//     }
    
    return sample_ids;
}
