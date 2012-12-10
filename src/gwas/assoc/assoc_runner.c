/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
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

#include "assoc_runner.h"

int run_association_test(shared_options_data_t* shared_options_data, assoc_options_data_t* options_data) {
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, INT_MAX, output_list);

    int ret_code = 0;
    vcf_file_t *file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
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
            LOG_DEBUG_F("Level %d: number of threads in the team - %d\n", 0, omp_get_num_threads());
            
            double start = omp_get_wtime();

            ret_code = vcf_read(file, 0,
                                (shared_options_data->batch_bytes > 0) ? shared_options_data->batch_bytes : shared_options_data->batch_lines,
                                shared_options_data->batch_bytes <= 0);

            double stop = omp_get_wtime();
            double total = stop - start;

            if (ret_code) {
                LOG_FATAL_F("Error %d while reading the file %s\n", ret_code, file->filename);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            notify_end_reading(file);
        }

#pragma omp section
        {
            LOG_DEBUG_F("Level %d: number of threads in the team - %d\n", 10, omp_get_num_threads());
            // Enable nested parallelism
            omp_set_nested(1);
            
            volatile int initialization_done = 0;
            individual_t **individuals;
            
            // Create chain of filters for the VCF file
            filter_t **filters = NULL;
            int num_filters = 0;
            if (shared_options_data->chain != NULL) {
                filters = sort_filter_chain(shared_options_data->chain, &num_filters);
            }
            FILE *passed_file = NULL, *failed_file = NULL;
            get_filtering_output_files(shared_options_data, &passed_file, &failed_file);
    
            double start = omp_get_wtime();

            double *factorial_logarithms = NULL;
            
#pragma omp parallel num_threads(shared_options_data->num_threads) shared(initialization_done, factorial_logarithms, filters, individuals)
            {
            LOG_DEBUG_F("Level %d: number of threads in the team - %d\n", 11, omp_get_num_threads()); 

            int i = 0;
            char *text_begin, *text_end;
            while(text_begin = fetch_vcf_text_batch(file)) {
                text_end = text_begin + strlen(text_begin);
                
                vcf_reader_status *status = vcf_reader_status_new(shared_options_data->batch_lines);
                if (shared_options_data->batch_bytes > 0) {
                    ret_code = run_vcf_parser(text_begin, text_end, 0, file, status);
                } else if (shared_options_data->batch_lines > 0) {
                    ret_code = run_vcf_parser(text_begin, text_end, shared_options_data->batch_lines, file, status);
                }
                
                // Initialize structures needed for association tests and write headers of output files
                if (!initialization_done) {
//                    sample_ids = associate_samples_and_positions(file);
# pragma omp critical
                {
                    // Guarantee that just one thread performs this operation
                    if (!initialization_done) {
                        // Sort individuals in PED as defined in the VCF file
                    	individuals = sort_individuals(file, ped_file);
                        
                        // Add headers associated to the defined filters
                        vcf_header_entry_t **filter_headers = get_filters_as_vcf_headers(filters, num_filters);
                        for (int j = 0; j < num_filters; j++) {
                            add_vcf_header_entry(filter_headers[j], file);
                        }
                        
                        // Write file format, header entries and delimiter
                        if (passed_file != NULL) { write_vcf_header(file, passed_file); }
                        if (failed_file != NULL) { write_vcf_header(file, failed_file); }
                        
                        LOG_DEBUG("VCF header written\n");
                        
                        if (options_data->task == FISHER) {
                            factorial_logarithms = init_logarithm_array(get_num_vcf_samples(file) * 10);
                        }
                        
                        initialization_done = 1;
                    }
                }
                }
                
                vcf_batch_t *batch = fetch_vcf_batch(file);
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

                // Launch association test over records that passed the filters
                if (passed_records->size > 0) {
//                     LOG_DEBUG_F("[%d] Test execution\n", omp_get_thread_num());
                	assert(individuals);
                    assoc_test(options_data->task, (vcf_record_t**) passed_records->items, passed_records->size, 
                                individuals, get_num_vcf_samples(file), factorial_logarithms, output_list);
                }
                
                // Write records that passed and failed to separate files
                write_filtering_output_files(passed_records, failed_records, passed_file, failed_file);
                
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
                
                i++;
            }  
            
            notify_end_parsing(file);
            }

            double stop = omp_get_wtime();
            double total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Free resources
            for (int i = 0; i < num_filters; i++) {
                filter_t *filter = filters[i];
                filter->free_func(filter);
            }
            free(filters);
            free(individuals);

            // Decrease list writers count
            for (int i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }

#pragma omp section
        {
            LOG_DEBUG_F("Level %d: number of threads in the team - %d\n", 20, omp_get_num_threads());
            
            // Thread which writes the results to the output file
            FILE *fd = NULL;
            char *path = NULL, *filename = NULL;
            size_t filename_len = 0;
            
            // Set whole path to the output file
            if (options_data->task == CHI_SQUARE) {
                filename_len = strlen("hpg-variant.chisq");
                filename = strdup("hpg-variant.chisq");
            } else if (options_data->task == FISHER) {
                filename_len = strlen("hpg-variant.fisher");
                filename = strdup("hpg-variant.fisher");
            }
            path = (char*) calloc ((output_directory_len + filename_len + 2), sizeof(char));
            sprintf(path, "%s/%s", shared_options_data->output_directory, filename);
            fd = fopen(path, "w");
            
            LOG_INFO_F("Association test output filename = %s\n", path);
            free(filename);
            
            double start = omp_get_wtime();
            
            // Write data: header + one line per variant
            list_item_t* item = NULL;
            
            if (options_data->task == CHI_SQUARE) {
                assoc_basic_result_t *result;
                fprintf(fd, "#CHR         POS       A1      C_A1    C_U1         F_A1            F_U1       A2      C_A2    C_U2         F_A2            F_U2              OR           CHISQ         P-VALUE\n");
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
                fprintf(fd, "#CHR         POS       A1      C_A1    C_U1         F_A1            F_U1       A2      C_A2    C_U2         F_A2            F_U2              OR         P-VALUE\n");
                while ((item = list_remove_item(output_list)) != NULL) {
                    result = item->data_p;
                    
                    fprintf(fd, "%s\t%8ld\t%s\t%3d\t%3d\t%6f\t%6f\t%s\t%3d\t%3d\t%6f\t%6f\t%6f\t%6f\n",
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
            
            // Sort resulting file
            char *cmd = calloc (40 + strlen(path) * 4, sizeof(char));
            sprintf(cmd, "sort -k1,1h -k2,2n %s > %s.tmp && mv %s.tmp %s", path, path, path, path);
            
            int sort_ret = system(cmd);
            if (sort_ret) {
                LOG_WARN("TDT results could not be sorted by chromosome and position, will be shown unsorted\n");
            }
            
            double stop = omp_get_wtime();
            double total = stop - start;

            LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        }
    }
   
    free(output_list);
    vcf_close(file);
    // TODO delete conflicts among frees
    ped_close(ped_file, 0);
        
    return ret_code;
}


individual_t **sort_individuals(vcf_file_t *vcf, ped_file_t *ped) {
	family_t *family;
	family_t **families = (family_t**) cp_hashtable_get_values(ped->families);
	int num_families = get_num_families(ped);

	individual_t **individuals = calloc (get_num_vcf_samples(vcf), sizeof(individual_t*));
	cp_hashtable *positions = associate_samples_and_positions(vcf);
	int *pos;

	for (int f = 0; f < num_families; f++) {
		family = families[f];
		individual_t *father = family->father;
		individual_t *mother = family->mother;
		cp_list *children = family->children;

		if (father != NULL) {
			pos = cp_hashtable_get(positions, father->id);
			individuals[*pos] = father;
		}

		if (mother != NULL) {
			pos = cp_hashtable_get(positions, mother->id);
			individuals[*pos] = mother;
		}

		cp_list_iterator *children_iterator = cp_list_create_iterator(family->children, COLLECTION_LOCK_READ);
		individual_t *child = NULL;
		while ((child = cp_list_iterator_next(children_iterator)) != NULL) {
			pos = cp_hashtable_get(positions, child->id);
			individuals[*pos] = child;
		}
        cp_list_iterator_destroy(children_iterator);
	}

	cp_hashtable_destroy(positions);

	return individuals;
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
//     char **keys = (char**) cp_hashtable_get_keys(sample_ids);
//     int num_keys = cp_hashtable_count(sample_ids);
//     for (int i = 0; i < num_keys; i++) {
//         printf("%s\t%d\n", keys[i], *((int*) cp_hashtable_get(sample_ids, keys[i])));
//     }
    
    return sample_ids;
}
