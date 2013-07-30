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

#include "effect_runner.h"


int run_effect(char **urls, shared_options_data_t *shared_options_data, effect_options_data_t *options_data) {
    int ret_code = 0;
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
    
    char *output_directory = shared_options_data->output_directory;
    size_t output_directory_len = strlen(output_directory);
    
    ret_code = create_directory(output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", output_directory);
    }
    
    // Remove all .txt files in folder
    ret_code = delete_files_by_extension(output_directory, "txt");
    if (ret_code != 0) {
        return ret_code;
    }
    
    // Initialize environment for connecting to the web service
    ret_code = init_http_environment(0);
    if (ret_code != 0) {
        return ret_code;
    }
    
    // Output file descriptors
    static cp_hashtable *output_files = NULL;
    // Lines of the output data in the main .txt files
    static list_t *output_list = NULL;
    // Consequence type counters (for summary, must be kept between web service calls)
    static cp_hashtable *summary_count = NULL;
    // Gene list (for genes-with-variants, must be kept between web service calls)
    static cp_hashtable *gene_list = NULL;

    // Initialize collections of file descriptors and summary counters
    ret_code = initialize_output_files(output_directory, output_directory_len, &output_files);
    if (ret_code != 0) {
        return ret_code;
    }
    initialize_output_data_structures(shared_options_data, &output_list, &summary_count, &gene_list);
    initialize_ws_buffers(shared_options_data->num_threads);
    
    // Create job.status file
    char job_status_filename[output_directory_len + 10];
    sprintf(job_status_filename, "%s/job.status", output_directory);
    FILE *job_status = new_job_status_file(job_status_filename);
    if (!job_status) {
        LOG_FATAL("Can't create job status file\n");
    } else {
        update_job_status_file(0, job_status);
    }
    
 
#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            
            start = omp_get_wtime();
            
            ret_code = vcf_read(vcf_file, 1,
                                (shared_options_data->batch_bytes > 0) ? shared_options_data->batch_bytes : shared_options_data->batch_lines,
                                shared_options_data->batch_bytes <= 0);

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) {
                LOG_ERROR_F("Error %d while reading the file %s\n", ret_code, vcf_file->filename);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            notify_end_parsing(vcf_file);
        }
        
#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            // Filters and files for filtering output
            filter_t **filters = NULL;
            int num_filters = 0;
            if (shared_options_data->chain != NULL) {
                filters = sort_filter_chain(shared_options_data->chain, &num_filters);
            }
            FILE *passed_file = NULL, *failed_file = NULL, *non_processed_file = NULL;
            get_filtering_output_files(shared_options_data, &passed_file, &failed_file);
            
            // Pedigree information (used in some filters)
            individual_t **individuals;
            khash_t(ids) *sample_ids = NULL;
            
            // Filename structure outdir/vcfname.errors
            char *prefix_filename = calloc(strlen(shared_options_data->vcf_filename), sizeof(char));
            get_filename_from_path(shared_options_data->vcf_filename, prefix_filename);
            char *non_processed_filename = malloc((strlen(shared_options_data->output_directory) + strlen(prefix_filename) + 9) * sizeof(char));
            sprintf(non_processed_filename, "%s/%s.errors", shared_options_data->output_directory, prefix_filename);
            non_processed_file = fopen(non_processed_filename, "w");
            free(non_processed_filename);
            
            // Maximum size processed by each thread (never allow more than 1000 variants per query)
            if (shared_options_data->batch_lines > 0) {
                shared_options_data->entries_per_thread = MIN(MAX_VARIANTS_PER_QUERY, 
                            ceil((float) shared_options_data->batch_lines / shared_options_data->num_threads));
            } else {
                shared_options_data->entries_per_thread = MAX_VARIANTS_PER_QUERY;
            }
            LOG_DEBUG_F("entries-per-thread = %d\n", shared_options_data->entries_per_thread);
    
            int i = 0;
            vcf_batch_t *batch = NULL;
            int ret_ws_0 = 0, ret_ws_1 = 0, ret_ws_2 = 0;
            
            start = omp_get_wtime();

            while (batch = fetch_vcf_batch(vcf_file)) {
                if (i == 0) {
                    // Add headers associated to the defined filters
                    vcf_header_entry_t **filter_headers = get_filters_as_vcf_headers(filters, num_filters);
                    for (int j = 0; j < num_filters; j++) {
                        add_vcf_header_entry(filter_headers[j], vcf_file);
                    }
                        
                    // Write file format, header entries and delimiter
                    if (passed_file != NULL) { write_vcf_header(vcf_file, passed_file); }
                    if (failed_file != NULL) { write_vcf_header(vcf_file, failed_file); }
                    if (non_processed_file != NULL) { write_vcf_header(vcf_file, non_processed_file); }
                    
                    LOG_DEBUG("VCF header written\n");
                    
                    if (ped_file) {
                        // Create map to associate the position of individuals in the list of samples defined in the VCF file
                        sample_ids = associate_samples_and_positions(vcf_file);
                        // Sort individuals in PED as defined in the VCF file
                        individuals = sort_individuals(vcf_file, ped_file);
                    }
                }
                
//                     printf("batch loaded = '%.*s'\n", 50, batch->text);
//                     printf("batch text len = %zu\n", strlen(batch->text));

//                 if (i % 10 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->records->size, batch->records->capacity);
//                 }

                int reconnections = 0;
                int max_reconnections = 3; // TODO allow to configure?

                // Write records that passed to a separate file, and query the WS with them as args
                array_list_t *failed_records = NULL;
                int num_variables = ped_file? get_num_variables(ped_file): 0;
                array_list_t *passed_records = filter_records(filters, num_filters, individuals, sample_ids, num_variables, batch->records, &failed_records);
                if (passed_records->size > 0) {
                    // Divide the list of passed records in ranges of size defined in config file
                    int num_chunks;
                    int *chunk_sizes;
                    int *chunk_starts = create_chunks(passed_records->size, shared_options_data->entries_per_thread, &num_chunks, &chunk_sizes);
                    
                    do {
                        // OpenMP: Launch a thread for each range
                        #pragma omp parallel for num_threads(shared_options_data->num_threads)
                        for (int j = 0; j < num_chunks; j++) {
                            int tid = omp_get_thread_num();
                            LOG_DEBUG_F("[%d] WS invocation\n", tid);
                            LOG_DEBUG_F("[%d] -- effect WS\n", tid);
                            if (!reconnections || ret_ws_0) {
                                ret_ws_0 = invoke_effect_ws(urls[0], (vcf_record_t**) (passed_records->items + chunk_starts[j]), 
                                                            chunk_sizes[j], options_data->excludes);
                                parse_effect_response(tid, output_directory, output_directory_len, output_files, output_list, summary_count, gene_list);
                                free(effect_line[tid]);
                                effect_line[tid] = (char*) calloc (max_line_size[tid], sizeof(char));
                            }
                            
                            if (!options_data->no_phenotypes) {
                                if (!reconnections || ret_ws_1) {
                                    LOG_DEBUG_F("[%d] -- snp WS\n", omp_get_thread_num());
                                    ret_ws_1 = invoke_snp_phenotype_ws(urls[1], (vcf_record_t**) (passed_records->items + chunk_starts[j]), chunk_sizes[j]);
                                    parse_snp_phenotype_response(tid, output_list);
                                    free(snp_line[tid]);
                                    snp_line[tid] = (char*) calloc (snp_max_line_size[tid], sizeof(char));
                                }
                                 
                                if (!reconnections || ret_ws_2) {
                                    LOG_DEBUG_F("[%d] -- mutation WS\n", omp_get_thread_num());
                                    ret_ws_2 = invoke_mutation_phenotype_ws(urls[2], (vcf_record_t**) (passed_records->items + chunk_starts[j]), chunk_sizes[j]);
                                    parse_mutation_phenotype_response(tid, output_list);
                                    free(mutation_line[tid]);
                                    mutation_line[tid] = (char*) calloc (mutation_max_line_size[tid], sizeof(char));
                                }
                            }
                        }
                        
                        LOG_DEBUG_F("*** %dth web services invocation finished\n", i);
                        
                        if (ret_ws_0 || ret_ws_1 || ret_ws_2) {
                            if (ret_ws_0) {
                                LOG_ERROR_F("Effect web service error: %s\n", get_last_http_error(ret_ws_0));
                            }
                            if (ret_ws_1) {
                                LOG_ERROR_F("SNP phenotype web service error: %s\n", get_last_http_error(ret_ws_1));
                            }
                            if (ret_ws_2) {
                                LOG_ERROR_F("Mutations phenotype web service error: %s\n", get_last_http_error(ret_ws_2));
                            }
                            
                            // In presence of errors, wait 4 seconds before retrying
                            reconnections++;
                            LOG_ERROR_F("Some errors ocurred, reconnection #%d\n", reconnections);
                            sleep(4);
                        } else {
                            free(chunk_starts);
                            free(chunk_sizes);
                        }
                    } while (reconnections < max_reconnections && (ret_ws_0 || ret_ws_1 || ret_ws_2));
                }
                
                // If the maximum number of reconnections was reached still with errors, 
                // write the non-processed batch to the corresponding file
                if (reconnections == max_reconnections && (ret_ws_0 || ret_ws_1 || ret_ws_2)) {
                #pragma omp critical
                    {
                        write_vcf_batch(batch, non_processed_file);
                    }
                }
                
                // Write records that passed and failed filters to separate files, and free them
                write_filtering_output_files(passed_records, failed_records, passed_file, failed_file);
                free_filtered_records(passed_records, failed_records, batch->records);
                
                // Free batch and its contents
                vcf_batch_free(batch);
                
                i++;
            }

            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Free resources
            if (passed_file) { fclose(passed_file); }
            if (failed_file) { fclose(failed_file); }
            if (non_processed_file) { fclose(non_processed_file); }
            
            // Free filters
            for (i = 0; i < num_filters; i++) {
                filter_t *filter = filters[i];
                filter->free_func(filter);
            }
            free(filters);
            
            // Decrease list writers count
            for (i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }
        
#pragma omp section
        {
            // Thread which writes the results to all_variants, summary and one file per consequence type
            int ret = 0;
            char *line;
            list_item_t* item = NULL;
            FILE *fd = NULL;
            
            FILE *all_variants_file = cp_hashtable_get(output_files, "all_variants");
            FILE *snp_phenotype_file = cp_hashtable_get(output_files, "snp_phenotypes");
            FILE *mutation_phenotype_file = cp_hashtable_get(output_files, "mutation_phenotypes");
            
            while ((item = list_remove_item(output_list)) != NULL) {
                line = item->data_p;
                
                // Type greater than 0: consequence type identified by its SO code
                // Type equals to -1: SNP phenotype
                // Type equals to -2: mutation phenotype
                if (item->type > 0) {
                    // Write entry in the consequence type file
                    fd = cp_hashtable_get(output_files, &(item->type));
                    int ret = fprintf(fd, "%s\n", line);
                    if (ret < 0) {
                        LOG_ERROR_F("Error writing to file: '%s'\n", line);
                    }
                    
                    // Write in all_variants
                    ret = fprintf(all_variants_file, "%s\n", line);
                    if (ret < 0) {
                        LOG_ERROR_F("Error writing to all_variants: '%s'\n", line);
                    }
                    
                } else if (item->type == SNP_PHENOTYPE) {
                    ret = fprintf(snp_phenotype_file, "%s\n", line);
                    if (ret < 0) {
                        LOG_ERROR_F("Error writing to snp_phenotypes: '%s'\n", line);
                    }
                    
                } else if (item->type == MUTATION_PHENOTYPE) {
                    ret = fprintf(mutation_phenotype_file, "%s\n", line);
                    if (ret < 0) {
                        LOG_ERROR_F("Error writing to mutation_phenotypes: '%s'\n", line);
                    }
                }
                
                free(line);
                list_item_free(item);
            }
            
        }
    }

    write_summary_file(summary_count, cp_hashtable_get(output_files, "summary"));
    write_genes_with_variants_file(gene_list, output_directory);
    write_result_file(shared_options_data, options_data, summary_count, output_directory);

    free_output_data_structures(output_files, summary_count, gene_list);
    free_ws_buffers(shared_options_data->num_threads);
    free(output_list);
    vcf_close(vcf_file);
    
    update_job_status_file(100, job_status);
    close_job_status_file(job_status);
    
    return ret_code;
}



static void parse_effect_response(int tid, char *output_directory, size_t output_directory_len, cp_hashtable *output_files, 
                                  list_t *output_list, cp_hashtable *summary_count, cp_hashtable *gene_list) {
    int *SO_found = (int*) malloc (sizeof(int)); // Whether the SO code field has been found
    int *count;
    char tmp_consequence_type[128];
    
    int num_lines;
    char **split_batch = split(effect_line[tid], "\n", &num_lines);
    
    for (int i = 0; i < num_lines; i++) {
        int num_columns;
        char *copy_buf = strdup(split_batch[i]);
        char **split_result = split(copy_buf, "\t", &num_columns);
        free(copy_buf);
        
        // Find consequence type name (always after SO field)
        *SO_found = 0;
        if (num_columns == 25) {
//             LOG_DEBUG_F("gene = %s\tSO = %d\tCT = %s\n", split_result[17], atoi(split_result[18] + 3), split_result[19]);
            if (!cp_hashtable_contains(gene_list, split_result[17])) {
                cp_hashtable_put(gene_list, strdup(split_result[17]), NULL);
            }
            *SO_found = atoi(split_result[18] + 3);
           memset(tmp_consequence_type, 0, 128 * sizeof(char));
           strncat(tmp_consequence_type, split_result[19], strlen(split_result[19]));
        } else {
            if (strlen(split_batch[i]) == 0) { // Last line in batch could be only a newline
                continue;
            }
            
            LOG_INFO_F("[%d] Non-valid line found (%d fields): '%s'\n", tid, num_columns, split_batch[i]);
            
            for (int s = 0; s < num_columns; s++) {
                free(split_result[s]);
            }
            free(split_result);
            continue;
        }
        
        for (int s = 0; s < num_columns; s++) {
            free(split_result[s]);
        }
        free(split_result);
        
        if (!*SO_found) { // SO:000000 is not valid
            LOG_INFO_F("[%d] Non-valid SO found (0)\n", tid);
            continue;
        }

//         LOG_DEBUG_F("[%d] SO found = %d\n", tid, *SO_found);
        size_t consequence_type_len = strlen(tmp_consequence_type);
     
        // If file does not exist, create its descriptor and summary counter
        FILE *aux_file = cp_hashtable_get(output_files, SO_found);
        if (!aux_file) {
#pragma omp critical
            {
                // This construction avoids 2 threads trying to insert the same CT
                aux_file = cp_hashtable_get(output_files, SO_found);
                if (!aux_file) {
                    char filename[output_directory_len + consequence_type_len + 6];
                    memset(filename, 0, (output_directory_len + consequence_type_len + 6) * sizeof(char));
                    strncat(filename, output_directory, output_directory_len);
                    strncat(filename, "/", 1);
                    strncat(filename, tmp_consequence_type, consequence_type_len);
                    strncat(filename, ".txt", 4);
                    aux_file = fopen(filename, "a");
                    
                    // Add to hashtables (file descriptors and summary counters)
                    int *SO_stored = (int*) malloc (sizeof(int));
                    *SO_stored = *SO_found;
                    cp_hashtable_put(output_files, SO_stored, aux_file);

                    LOG_INFO_F("[%d] New consequence type found = %s\n", tid, tmp_consequence_type);
                }
            }
        }
        
        // Write line[tid] to file corresponding to its consequence type
        if (aux_file) { 
#pragma omp critical
            {
                // TODO move critical one level below?
                count = (int*) cp_hashtable_get(summary_count, tmp_consequence_type);
                if (count == NULL) {
                    char *consequence_type = (char*) calloc (consequence_type_len+1, sizeof(char));
                    strncat(consequence_type, tmp_consequence_type, consequence_type_len);
                    assert(!strcmp(consequence_type, tmp_consequence_type));
                    count = (int*) malloc (sizeof(int));
                    *count = 0;
                    cp_hashtable_put(summary_count, consequence_type, count);
//                     LOG_DEBUG_F("[%d] Initialized summary count for %s\n", tmp_consequence_type);
                }
                // Increment counter for summary
                (*count)++;
            }
            
//             LOG_DEBUG_F("[%d] before writing %s\n", tid, tmp_consequence_type);
            list_item_t *output_item = list_item_new(tid, *SO_found, strdup(split_batch[i]));
            list_insert_item(output_item, output_list);
//             LOG_DEBUG_F("[%d] after writing %s\n", tid, tmp_consequence_type);
        }
    }
    
    for (int i = 0; i < num_lines; i++) {
        free(split_batch[i]);
    }
    free(split_batch);
}

static void parse_snp_phenotype_response(int tid, list_t *output_list) {
    list_item_t *output_item = list_item_new(tid, SNP_PHENOTYPE, trim(strdup(snp_line[tid])));
    list_insert_item(output_item, output_list);
}

static void parse_mutation_phenotype_response(int tid, list_t *output_list) {
    list_item_t *output_item = list_item_new(tid, MUTATION_PHENOTYPE, trim(strdup(mutation_line[tid])));
    list_insert_item(output_item, output_list);
}



static int initialize_output_files(char *output_directory, size_t output_directory_len, cp_hashtable **output_files) {
    // Initialize collections of file descriptors
    *output_files = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                                  50,
                                                  cp_hash_int,
                                                  (cp_compare_fn) int_cmp,
                                                  NULL,
                                                  (cp_destructor_fn) free_file_key1,
                                                  NULL,
                                                  (cp_destructor_fn) free_file_descriptor
                                                 );
    
    char all_variants_filename[output_directory_len + 18];
    sprintf(all_variants_filename, "%s/all_variants.txt", output_directory);
    FILE *all_variants_file = fopen(all_variants_filename, "a");
    if (!all_variants_file) {   // Can't store results
        return 1;
    }
    char *key = (char*) calloc (13, sizeof(char));
    strncat(key, "all_variants", 12);
    cp_hashtable_put(*output_files, key, all_variants_file);
    
    char summary_filename[output_directory_len + 13];
    sprintf(summary_filename, "%s/summary.txt", output_directory);
    FILE *summary_file = fopen(summary_filename, "a");
    if (!summary_file) {   // Can't store results
        return 2;
    }
    key = (char*) calloc (8, sizeof(char));
    strncat(key, "summary", 7);
    cp_hashtable_put(*output_files, key, summary_file);
    
    char snp_phenotype_filename[output_directory_len + 20];
    sprintf(snp_phenotype_filename, "%s/snp_phenotypes.txt", output_directory);
    FILE *snp_phenotype_file = fopen(snp_phenotype_filename, "a");
    if (!snp_phenotype_filename) {
        return 3;
    }
    key = (char*) calloc (15, sizeof(char));
    strncat(key, "snp_phenotypes", 14);
    cp_hashtable_put(*output_files, key, snp_phenotype_file);
    
    char mutation_phenotype_filename[output_directory_len + 25];
    sprintf(mutation_phenotype_filename, "%s/mutation_phenotypes.txt", output_directory);
    FILE *mutation_phenotype_file = fopen(mutation_phenotype_filename, "a");
    if (!mutation_phenotype_filename) {
        return 3;
    }
    key = (char*) calloc (20, sizeof(char));
    strncat(key, "mutation_phenotypes", 19);
    cp_hashtable_put(*output_files, key, mutation_phenotype_file);
    
    return 0;
}

static void initialize_output_data_structures(shared_options_data_t *shared_options, list_t **output_list, cp_hashtable **summary_count, cp_hashtable **gene_list) {
    // Initialize output text list
    *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options->num_threads, shared_options->max_batches * shared_options->batch_lines, *output_list);
    
    // Initialize summary counters and genes list
    *summary_count = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                                 64,
                                                 cp_hash_istring,
                                                 (cp_compare_fn) strcasecmp,
                                                 NULL,
                                                 (cp_destructor_fn) free_file_key2,
                                                 NULL,
                                                 (cp_destructor_fn) free_summary_counter
                                                 );

    *gene_list = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                              64,
                                              cp_hash_istring,
                                              (cp_compare_fn) strcasecmp,
                                              NULL,
                                              (cp_destructor_fn) free_file_key2,
                                              NULL,
                                              NULL
                                             );
}

void free_output_data_structures(cp_hashtable *output_files, cp_hashtable *summary_count, cp_hashtable *gene_list) {
    // Free file descriptors, summary counters and gene list
    cp_hashtable_destroy(output_files);
    cp_hashtable_destroy(summary_count);
    cp_hashtable_destroy(gene_list);
}

static void free_file_key1(int *key) {
    LOG_DEBUG_F("Free file key 1: %d\n", *key);
    free(key);
}

static void free_file_descriptor(FILE *fd) {
    LOG_DEBUG("Free file descriptor\n");
    fclose(fd);
}

static void free_file_key2(char *key) {
    LOG_DEBUG_F("Free file key 2: %s\n", key);
    free(key);
}

static void free_summary_counter(int *count) {
    LOG_DEBUG_F("Free summary counter %d\n", *count);
    free(count);
}

static int int_cmp(int *a, int *b) {
    return *a - *b;
}
