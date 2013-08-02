/*
 * Copyright (c) 2012-2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
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

#include "filter_runner.h"

int run_filter(shared_options_data_t *shared_options_data, filter_options_data_t *options_data) {
    int ret_code;
    double start, stop, total;
    FILE *passed_file = NULL, *failed_file = NULL;
    // List that stores the batches of records filtered by each thread
    list_t *passed_list[shared_options_data->num_threads];
    // List that stores which thread filtered the next batch to save
    list_t *next_token_list = malloc(sizeof(list_t));
    
    vcf_file_t *vcf_file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
    if (!vcf_file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ped_file_t *ped_file = NULL;
    if (shared_options_data->ped_filename) {
        ped_file = ped_open(shared_options_data->ped_filename);
        if (!ped_file) {
            LOG_FATAL("PED file does not exist!\n");
        }
        LOG_INFO("About to read PED file...\n");
        // Read PED file before doing any processing
        ret_code = ped_read(ped_file);
        if (ret_code != 0) {
            LOG_FATAL_F("Can't read PED file: %s\n", ped_file->filename);
        }
    }
    
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }
    
    // Initialize variables related to the different threads
    for (int i = 0; i < shared_options_data->num_threads; i++) {
        passed_list[i] = (list_t*) malloc(sizeof(list_t));
        list_init("input", 1, shared_options_data->max_batches, passed_list[i]);
    }
    list_init("next_token", shared_options_data->num_threads, shared_options_data->max_batches, next_token_list);
    
    get_filtering_output_files(shared_options_data, &passed_file, &failed_file);
    if (!options_data->save_rejected) {
        fclose(failed_file);
    }
    LOG_DEBUG("Output files created\n");
    
#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            if (shared_options_data->batch_bytes > 0) {
                ret_code = vcf_parse_batches_in_bytes(shared_options_data->batch_bytes, vcf_file);
            } else if (shared_options_data->batch_lines > 0) {
                ret_code = vcf_parse_batches(shared_options_data->batch_lines, vcf_file);
            }

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            notify_end_parsing(vcf_file);
        }
        
#pragma omp section
        {
            // Enable nested parallelism
            omp_set_nested(1);
            
            filter_t **filters = NULL;
            int num_filters = 0;
            if (shared_options_data->chain != NULL) {
                filters = sort_filter_chain(shared_options_data->chain, &num_filters);
            }
            
            volatile int initialization_done = 0;
            individual_t **individuals = NULL;
            khash_t(ids) *sample_ids = NULL;
            
            start = omp_get_wtime();

            volatile int i = 0;
            
#pragma omp parallel num_threads(shared_options_data->num_threads) shared(i)
            {
            vcf_batch_t *batch = NULL;
            int index = omp_get_thread_num() % shared_options_data->num_threads;
            
            while ( batch = fetch_vcf_batch(vcf_file) ) {
                // Initialize structures needed for filtering
                if (!initialization_done) {
#pragma omp critical
                {
                    // Guarantee that just one thread performs this operation
                    if (!initialization_done) {
                        // Add headers associated to the defined filters
                        vcf_header_entry_t **filter_headers = get_filters_as_vcf_headers(filters, num_filters);
                        for (int j = 0; j < num_filters; j++) {
                            add_vcf_header_entry(filter_headers[j], vcf_file);
                        }

                        // Write file format, header entries and delimiter
                        write_vcf_header(vcf_file, passed_file);
                        if (options_data->save_rejected) {
                            write_vcf_header(vcf_file, failed_file);
                        }

                        LOG_DEBUG("VCF headers written\n");

                        if (ped_file) {
                            // Create map to associate the position of individuals in the list of samples defined in the VCF file
                            sample_ids = associate_samples_and_positions(vcf_file);
                            // Sort individuals in PED as defined in the VCF file
                            individuals = sort_individuals(vcf_file, ped_file);
                        }
                        
                        initialization_done = 1;
                    }
                }
                }
                
                array_list_t *input_records = batch->records;
                array_list_t *passed_records = NULL, *failed_records = NULL;

                if (i % 100 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->records->size, batch->records->capacity);
                }
                
                #pragma omp atomic
                i++;
                
                int num_variables = ped_file? get_num_variables(ped_file): 0;
                if (filters == NULL) {
                    passed_records = input_records;
                } else {
                    failed_records = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
                    passed_records = run_filter_chain(input_records, failed_records, individuals, sample_ids,num_variables, filters, num_filters);
                }
                
                filter_temp_output_t *output = filter_temp_output_new(batch, passed_records, failed_records);
                list_item_t *output_item = list_item_new(i, 0, output);
                list_insert_item(output_item, passed_list[index]);
                
                list_item_t *token_item = list_item_new(index, 0, NULL);
                list_insert_item(token_item, next_token_list);
                
            }
            }

            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Notify end of operations
            for (int i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(next_token_list);
                list_decr_writers(passed_list[i]);
            }
            
            if (sample_ids) { kh_destroy(ids, sample_ids); }
            free(individuals);
            
            free_filters(filters, num_filters);
        }
        
#pragma omp section
        {
            int i = 0;
            list_item_t *token_item = NULL, *output_item = NULL;
            while ( token_item = list_remove_item(next_token_list) ) {
                output_item = list_remove_item(passed_list[token_item->id]);
                filter_temp_output_t *output = output_item->data_p;
                
                assert(output);
                
                // Write records that passed and failed to 2 new separated files
                if (output->passed_records != NULL && output->passed_records->size > 0) {
                    LOG_DEBUG_F("[batch %d] %zu passed records\n", i, output->passed_records->size);
                    for (int r = 0; r < output->passed_records->size; r++) {
                        write_vcf_record(output->passed_records->items[r], passed_file);
                    }
                }

                if (options_data->save_rejected && output->failed_records != NULL && output->failed_records->size > 0) {
                    LOG_DEBUG_F("[batch %d] %zu failed records\n", i, output->failed_records->size);
                    for (int r = 0; r < output->failed_records->size; r++) {
                        write_vcf_record(output->failed_records->items[r], failed_file);
                    }
                }

                // Free batch and its contents
                filter_temp_output_free(output);
                list_item_free(output_item);
                list_item_free(token_item);
                
                i++;
            }
              
        }
    }
    
    // Close files
    if (passed_file) {
        fclose(passed_file);
    }
    if (options_data->save_rejected && failed_file) {
        fclose(failed_file);
    }

    vcf_close(vcf_file);
    if (ped_file) { ped_close(ped_file, 1, 1); }
    
    return 0;
}


/* ******************************
 *           Auxiliary          *
 * ******************************/

filter_temp_output_t *filter_temp_output_new(vcf_batch_t *batch, array_list_t *passed_records, array_list_t *failed_records) {
    filter_temp_output_t *ret = malloc (sizeof(filter_temp_output_t));
    ret->batch = batch;
    ret->passed_records = passed_records;
    ret->failed_records = failed_records;
    return ret;
}

void filter_temp_output_free(filter_temp_output_t *temp) {
    assert(temp);
    array_list_free(temp->passed_records, NULL);
    array_list_free(temp->failed_records, NULL);
    vcf_batch_free(temp->batch);
}
