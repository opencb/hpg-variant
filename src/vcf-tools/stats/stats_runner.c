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

#include "stats.h"


int run_stats(shared_options_data_t *shared_options_data, stats_options_data_t *options_data) {
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, MIN(10, shared_options_data->max_batches) * shared_options_data->batch_lines, output_list);
    file_stats_t *file_stats = file_stats_new();
    sample_stats_t **sample_stats;

    int ret_code;
    double start, stop, total;
    
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
    
    LOG_INFO("About to retrieve statistics from VCF file...\n");

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
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            individual_t **individuals;
            khash_t(ids) *sample_ids = NULL;
            
            start = omp_get_wtime();
            
            int i = 0;
            vcf_batch_t *batch = NULL;
            while ((batch = fetch_vcf_batch(vcf_file)) != NULL) {
                if (i == 0) {
                    sample_stats = malloc (get_num_vcf_samples(vcf_file) * sizeof(sample_stats_t*));
                    for (int j = 0; j < get_num_vcf_samples(vcf_file); j++) {
                        sample_stats[j] = sample_stats_new(array_list_get(j, vcf_file->samples_names));
                    }
                    
                    if (ped_file) {
                        // Create map to associate the position of individuals in the list of samples defined in the VCF file
                        sample_ids = associate_samples_and_positions(vcf_file);
                        // Sort individuals in PED as defined in the VCF file
                        individuals = sort_individuals(vcf_file, ped_file);
                    }
                }
                
                if (i % 50 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->records->size, batch->records->capacity);
                }

                // Divide the list of passed records in ranges of size defined in config file
                int num_chunks;
                int *chunk_sizes;
                array_list_t *input_records = batch->records;
                int *chunk_starts = create_chunks(input_records->size, shared_options_data->entries_per_thread, &num_chunks, &chunk_sizes);
                
                // OpenMP: Launch a thread for each range
                #pragma omp parallel for num_threads(shared_options_data->num_threads)
                for (int j = 0; j < num_chunks; j++) {
                    LOG_DEBUG_F("[%d] Stats invocation\n", omp_get_thread_num());
                    // Invoke variant stats and/or sample stats when applies
                    if (options_data->variant_stats) {
                        ret_code = get_variants_stats((vcf_record_t**) (input_records->items + chunk_starts[j]), 
                                                      chunk_sizes[j], individuals, sample_ids, output_list, file_stats);
                    }
                    if (options_data->sample_stats) {
                        ret_code |= get_sample_stats((vcf_record_t**) (input_records->items + chunk_starts[j]), 
                                                      chunk_sizes[j], individuals, sample_ids, sample_stats, file_stats);
                    }
                }
                
                free(chunk_starts);
                free(chunk_sizes);
                vcf_batch_free(batch);
                
                i++;
            }
            
            // Write sample statistics
            char *stats_filename;
            FILE *stats_fd;
            
            if (options_data->sample_stats) {
                stats_filename = get_sample_stats_output_filename(shared_options_data->vcf_filename, 
                                                                  shared_options_data->output_filename, 
                                                                  shared_options_data->output_directory);
                if (!(stats_fd = fopen(stats_filename, "w"))) {
                    LOG_FATAL_F("Can't open file for writing statistics of samples: %s\n", stats_filename);
                }
                
                report_sample_variant_stats_header(stats_fd);
                report_sample_stats(stats_fd, NULL, vcf_file->samples_names->size, sample_stats);
                
                // Close sample stats file
                free(stats_filename);
                if (stats_fd != NULL) { fclose(stats_fd); }
            }
            
            stop = omp_get_wtime();
            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
            
            // Decrease list writers count
            for (i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
            
            if (sample_ids) { kh_destroy(ids, sample_ids); }
            if (individuals) { free(individuals); }
        }
        
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d writes the output\n", omp_get_thread_num());
            char *stats_filename, *summary_filename;
            FILE *stats_fd, *summary_fd;
    
            // Write variant statistics
            if (options_data->variant_stats) {
                stats_filename = get_variant_stats_output_filename(shared_options_data->vcf_filename, 
                                                                   shared_options_data->output_filename, 
                                                                   shared_options_data->output_directory);
                if (!(stats_fd = fopen(stats_filename, "w"))) {
                    LOG_FATAL_F("Can't open file for writing statistics of variants: %s\n", stats_filename);
                }
                
                // Write header
                report_vcf_variant_stats_header(stats_fd);
                
                // For each variant, generate a new line
                list_item_t* item = NULL;
                variant_stats_t *var_stats;
                while ((item = list_remove_item(output_list))) {
                    var_stats = item->data_p;
                    
                    report_vcf_variant_stats(stats_fd, NULL, var_stats);
                    
                    // Free resources
                    variant_stats_free(var_stats);
                    list_item_free(item);
                }
                
                // Write whole file stats (data only got when launching variant stats)
                summary_filename = get_vcf_file_stats_output_filename(shared_options_data->vcf_filename, 
                                                                      shared_options_data->output_filename, 
                                                                      shared_options_data->output_directory);
                if (!(summary_fd = fopen(summary_filename, "w"))) {
                    LOG_FATAL_F("Can't open file for writing statistics summary: %s\n", summary_filename);
                }
                report_vcf_summary_stats(summary_fd, NULL, file_stats);
                
                free(stats_filename);
                free(summary_filename);
                
                // Close variant stats file
                if (stats_fd) { fclose(stats_fd); }
                if (summary_fd && summary_fd != stdout) { fclose(summary_fd); }
            }
        }
        
    }
    
    vcf_close(vcf_file);
    if (ped_file) { ped_close(ped_file, 1); }
    free(sample_stats);
    free(file_stats);
//     free(read_list);
    free(output_list);
    
    return 0;
}
