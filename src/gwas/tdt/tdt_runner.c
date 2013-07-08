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

#include "tdt_runner.h"

int run_tdt_test(shared_options_data_t* shared_options_data) {
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, INT_MAX, output_list);

    int ret_code = 0;
    vcf_file_t *vcf_file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
    if (!vcf_file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ped_file_t *ped_file = ped_open(shared_options_data->ped_filename);
    if (!ped_file) {
        LOG_FATAL("PED file does not exist!\n");
    }
    
    LOG_INFO("About to read PED file...\n");
    // Read PED file before doing any processing
    ret_code = ped_read(ped_file);
    if (ret_code != 0) {
        LOG_FATAL_F("Can't read PED file: %s\n", ped_file->filename);
    }
    
    // Try to create the directory where the output files will be stored
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }
    
    LOG_INFO("About to perform TDT test...\n");

#pragma omp parallel sections private(ret_code)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Level %d: number of threads in the team - %d\n", 0, omp_get_num_threads());
            
            double start = omp_get_wtime();

            ret_code = vcf_read(vcf_file, 0,
                                (shared_options_data->batch_bytes > 0) ? shared_options_data->batch_bytes : shared_options_data->batch_lines,
                                shared_options_data->batch_bytes <= 0);

            double stop = omp_get_wtime();

            if (ret_code) {
                LOG_FATAL_F("Error %d while reading the file %s\n", ret_code, vcf_file->filename);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), stop - start);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), (stop - start) * 1000);

            notify_end_reading(vcf_file);
        }

#pragma omp section
        {
            LOG_DEBUG_F("Level %d: number of threads in the team - %d\n", 10, omp_get_num_threads());
            
            // Enable nested parallelism
            omp_set_nested(1);
            
            volatile int initialization_done = 0;
            // Pedigree information
            family_t **families = (family_t**) cp_hashtable_get_values(ped_file->families);
            int num_families = get_num_families(ped_file);
            individual_t **individuals = NULL;
            khash_t(ids) *sample_ids = NULL;
            
            // Create chain of filters for the VCF file
            filter_t **filters = NULL;
            int num_filters = 0;
            if (shared_options_data->chain != NULL) {
                filters = sort_filter_chain(shared_options_data->chain, &num_filters);
            }
            FILE *passed_file = NULL, *failed_file = NULL;
            get_filtering_output_files(shared_options_data, &passed_file, &failed_file);
    
            double start = omp_get_wtime();
            
            int i = 0;
//#pragma omp parallel num_threads(shared_options_data->num_threads) shared(initialization_done, families, sample_ids, filters)
#pragma omp parallel num_threads(shared_options_data->num_threads)
            {
            LOG_DEBUG_F("Level %d: number of threads in the team - %d\n", 11, omp_get_num_threads());
            
            char *text_begin, *text_end;
            vcf_reader_status *status;
            while(text_begin = fetch_vcf_text_batch(vcf_file)) {
                text_end = text_begin + strlen(text_begin);
                if (text_begin == text_end) { // EOF
                    free(text_begin);
                    break;
                }
                
#pragma omp critical 
                {
                    status = vcf_reader_status_new(shared_options_data->batch_lines, i);
                    i++;
                }
                
                if (shared_options_data->batch_bytes > 0) {
                    ret_code = run_vcf_parser(text_begin, text_end, 0, vcf_file, status);
                } else if (shared_options_data->batch_lines > 0) {
                    ret_code = run_vcf_parser(text_begin, text_end, shared_options_data->batch_lines, vcf_file, status);
                }
                
                if (ret_code > 0) {
                    LOG_FATAL_F("Error %d while parsing the file %s\n", ret_code, vcf_file->filename);
                }
                
                // Initialize structures needed for TDT and write headers of output files
                if (!initialization_done && vcf_file->samples_names->size > 0) {
#pragma omp critical
                {
                    // Guarantee that just one thread performs this operation
                    if (!initialization_done) {
                        // Create map to associate the position of individuals in the list of samples defined in the VCF file
                        sample_ids = associate_samples_and_positions(vcf_file);
                        // Sort individuals in PED as defined in the VCF file
                        individuals = sort_individuals(vcf_file, ped_file);
                        
                        // Add headers associated to the defined filters
                        vcf_header_entry_t **filter_headers = get_filters_as_vcf_headers(filters, num_filters);
                        for (int j = 0; j < num_filters; j++) {
                            add_vcf_header_entry(filter_headers[j], vcf_file);
                        }
                        
                        // Write file format, header entries and delimiter
                        if (passed_file != NULL) { write_vcf_header(vcf_file, passed_file); }
                        if (failed_file != NULL) { write_vcf_header(vcf_file, failed_file); }
                        
                        LOG_DEBUG("VCF header written\n");
                        
                        initialization_done = 1;
                    }
                }
                }
                
                // If it has not been initialized it means that header is not fully read
                if (!initialization_done) {
                    continue;
                }
                
                vcf_batch_t *batch = fetch_vcf_batch(vcf_file);

                if (i % 100 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->records->size, batch->records->capacity);
                }

                // Launch TDT test over records that passed the filters
                array_list_t *failed_records = NULL;
                assert(batch);
                assert(batch->records);
                array_list_t *passed_records = filter_records(filters, num_filters, individuals, sample_ids, batch->records, &failed_records);
                if (passed_records->size > 0) {
                    ret_code = tdt_test((vcf_record_t**) passed_records->items, passed_records->size, families, num_families, sample_ids, output_list);
                    if (ret_code) {
                        LOG_FATAL_F("[%d] Error in execution #%d of TDT\n", omp_get_thread_num(), i);
                    }
                }
                
                // Write records that passed and failed filters to separate files, and free them
                write_filtering_output_files(passed_records, failed_records, passed_file, failed_file);
                free_filtered_records(passed_records, failed_records, batch->records);
                
                // Free batch and its contents
                vcf_reader_status_free(status);
                vcf_batch_free(batch);
            }
            
            notify_end_parsing(vcf_file);
            }

            double stop = omp_get_wtime();
            
            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), stop - start);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), (stop - start) * 1000);

            // Free resources
            if (filters) {
                for (int i = 0; i < num_filters; i++) {
                    filter_t *filter = filters[i];
                    filter->free_func(filter);
                }
                free(filters);
            }
            
            if (sample_ids) { kh_destroy(ids, sample_ids); }
            if (individuals) { free(individuals); }
            free(families);
            
            // Decrease list writers count
            for (int i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }

#pragma omp section
        {
            // Thread which writes the results to the output file
            LOG_DEBUG_F("Level %d: number of threads in the team - %d\n", 20, omp_get_num_threads());
            
            // Get the file descriptor
            char *path;
            FILE *fd = get_output_file(shared_options_data, "hpg-variant.tdt", &path);
            LOG_INFO_F("TDT output filename = %s\n", path);
            
            double start = omp_get_wtime();
            
            // Write data: header + one line per variant
            write_output_header(fd);
            write_output_body(output_list, fd);
            
            fclose(fd);
            
            // Sort resulting file
            char *cmd = calloc (40 + strlen(path) * 4, sizeof(char));
            sprintf(cmd, "sort -k1,1h -k2,2n %s > %s.tmp && mv %s.tmp %s", path, path, path, path);
            
            int sort_ret = system(cmd);
            if (sort_ret) {
                LOG_WARN("TDT results could not be sorted by chromosome and position, will be shown unsorted\n");
            }
            
            free(cmd);
            free(path);
            
            double stop = omp_get_wtime();

            LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), stop - start);
            LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), (stop - start) * 1000);

        }
    }
    
    free(output_list);
    vcf_close(vcf_file);
    ped_close(ped_file, 1);
    
    return ret_code;
}


/* *******************
 * Output generation *
 * *******************/

void write_output_header(FILE *fd) {
    assert(fd);
    fprintf(fd, "#CHR         POS               ID      A1      A2         T       U           OR           CHISQ         P-VALUE\n");
}

void write_output_body(list_t* output_list, FILE *fd) {
    assert(fd);
    list_item_t* item = NULL;
    while (item = list_remove_item(output_list)) {
        tdt_result_t *result = item->data_p;
        
        fprintf(fd, "%s\t%8ld\t%s\t%s\t%s\t%3d\t%3d\t%6f\t%6f\t%6f\n",
                result->chromosome, result->position, result->id, result->reference, result->alternate, 
                result->t1, result->t2, result->odds_ratio, result->chi_square, result->p_value);
        
        tdt_result_free(result);
        list_item_free(item);
    }
}

