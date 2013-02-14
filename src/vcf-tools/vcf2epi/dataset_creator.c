#include "dataset_creator.h"


int create_dataset_from_vcf(shared_options_data_t* shared_options_data) {
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, INT_MAX, output_list);

    int ret_code = 0;
    double start, stop, total;
    vcf_file_t *file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
    if (!file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ped_file_t *ped_file = ped_open(shared_options_data->ped_filename);
    if (!ped_file) {
        LOG_FATAL("PED file does not exist!\n");
    }
    
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
    
    // Create dataset
    size_t num_variants = 0;
    int *destination, num_affected, num_unaffected;
    uint8_t *phenotypes;
    
    LOG_INFO("About to create epistasis dataset...\n");

#pragma omp parallel sections private(ret_code, start, stop, total) shared(destination,num_affected, num_unaffected)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Level %d: number of threads in the team - %d\n", 0, omp_get_num_threads());
            
            double start = omp_get_wtime();

            ret_code = vcf_read(file, 1,
                                (shared_options_data->batch_bytes > 0) ? shared_options_data->batch_bytes : shared_options_data->batch_lines,
                                shared_options_data->batch_bytes <= 0);

            double stop = omp_get_wtime();

            if (ret_code) {
                LOG_FATAL_F("Error %d while reading the file %s\n", ret_code, file->filename);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), stop - start);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), (stop - start) * 1000);

            notify_end_parsing(file);
        }

#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            filter_t **filters = NULL;
            int num_filters = 0;
            if (shared_options_data->chain != NULL) {
                filters = sort_filter_chain(shared_options_data->chain, &num_filters);
            }
            FILE *passed_file = NULL, *failed_file = NULL, *non_processed_file = NULL;
            get_filtering_output_files(shared_options_data, &passed_file, &failed_file);
    
            int i = 0;
            vcf_batch_t *batch = NULL;
            
            start = omp_get_wtime();

            while (batch = fetch_vcf_batch(file)) {
                if (i == 0) {
                    // Get individual phenotypes
                    phenotypes = get_individual_phenotypes(file, ped_file, &num_affected, &num_unaffected);
                    destination = group_individuals_by_phenotype(phenotypes, num_affected, num_unaffected);
                    
//                     assert(destination);
//                     
//                     printf("destination = { ");
//                     for (int k = 0; k < get_num_vcf_samples(file); k++) {
//                         printf("%d ", destination[k]);
//                     }
//                     printf("}\n");
                    
                    // Add headers associated to the defined filters
                    vcf_header_entry_t **filter_headers = get_filters_as_vcf_headers(filters, num_filters);
                    for (int j = 0; j < num_filters; j++) {
                        add_vcf_header_entry(filter_headers[j], file);
                    }
                    
                    // Write file format, header entries and delimiter
                    if (passed_file != NULL) { write_vcf_header(file, passed_file); }
                    if (failed_file != NULL) { write_vcf_header(file, failed_file); }
                }
                
                if (i % 10 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->records->size, batch->records->capacity);
                }

                // Write records that passed to a separate file, and query the WS with them as args
                array_list_t *failed_records = NULL;
                array_list_t *passed_records = filter_records(filters, num_filters, batch->records, &failed_records);
                if (passed_records->size > 0) {
                    // Divide the list of passed records in ranges of size defined in config file
//                     int num_chunks;
//                     int *chunk_sizes;
//                     int *chunk_starts = create_chunks(passed_records->size, shared_options_data->entries_per_thread, &num_chunks, &chunk_sizes);
//                     
//                     // OpenMP: Launch a thread for each range
//                     #pragma omp parallel for num_threads(shared_options_data->num_threads)
//                     for (int j = 0; j < num_chunks; j++) {
//                         // TODO task
//                     }

                    num_variants += passed_records->size;
                    uint8_t *genotypes = epistasis_dataset_process_records((vcf_record_t**) passed_records->items, passed_records->size, 
                                                                           destination, get_num_vcf_samples(file),
                                                                           shared_options_data->num_threads);
                    list_item_t *gt_item = list_item_new(1, passed_records->size, genotypes);
                    list_insert_item(gt_item, output_list);
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
            // Thread which writes the results to the output file
            LOG_DEBUG_F("Level %d: number of threads in the team - %d\n", 20, omp_get_num_threads());
            
            FILE *fp = fopen("epistasis_dataset.bin","wb");
            
            double start = omp_get_wtime();
            
            // TODO write binary file with dataset
            
            size_t num_samples;
            bool header_written = false;
            list_item_t* item = NULL;
            while (item = list_remove_item(output_list)) {
                uint8_t *genotypes = item->data_p;
                
                // First make room for the number of variants, then write the number of samples and phenotypes
                if (!header_written) {
                    num_samples = get_num_vcf_samples(file);
                    if (!fwrite(&num_variants, sizeof(size_t), 1, fp) ||
                        !fwrite(&num_affected, sizeof(uint32_t), 1, fp) ||
                        !fwrite(&num_unaffected, sizeof(uint32_t), 1, fp)) {
                        LOG_ERROR("The header of the dataset could not be written!");
                    }
                    header_written = true;
                }
                
                // Then the dataset itself
                if (!fwrite(genotypes, sizeof(uint8_t), item->type * num_samples, fp)) {
                    LOG_ERROR_F("%d variants could not be written!\n", item->type);
                }
                
                free(genotypes);
                list_item_free(item);
            }
            
            // Finally, write the real number of variants
            fseek(fp, 0, SEEK_SET);
            if (!fwrite(&num_variants, sizeof(size_t), 1, fp)) {
                LOG_ERROR("The number of variants in the dataset could not be written!");
            }
            
            fclose(fp);
    
            double stop = omp_get_wtime();

            LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), stop - start);
            LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), (stop - start) * 1000);

        }
    }
    
    free(phenotypes);
    free(output_list);
    vcf_close(file);
    // TODO delete conflicts among frees
//     ped_close(ped_file, 0);
    
    
    return ret_code;
}


/* *******************
 *  Dataset reading  *
 * *******************/

uint8_t *epistasis_dataset_process_records(vcf_record_t **variants, size_t num_variants, int *destination,
                                           int num_samples, int threads) {
    uint8_t *genotypes = malloc (num_variants * num_samples * sizeof(uint8_t));
    
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < num_variants; i++) {
        vcf_record_t *record = variants[i];
        int gt_position = get_field_position_in_format("GT", strndup(record->format, record->format_len));
        
        // For each sample, get genotype and increment counters in the dataset
        for (int k = 0; k < num_samples; k++) {
//             printf("%d\tbase = %d\tindex = %d\n", k, destination[k], i * num_samples + destination[k]);
            char *sample = strdup(array_list_get(k, record->samples));
            int allele1, allele2, gt_dataset_index;
            if (get_alleles(sample, gt_position, &allele1, &allele2)) {
                genotypes[i * num_samples + destination[k]] = 255;
                
            } else {
                if (!allele1 && !allele2) { // Homozygous in first allele
                    genotypes[i * num_samples + destination[k]] = 0;
                } else if (allele1 != allele2) { // Heterozygous
                    genotypes[i * num_samples + destination[k]] = 1;
                } else if (allele1 && allele1 == allele2) { // Homozygous in second allele
                    genotypes[i * num_samples + destination[k]] = 2;
                }
            }
            free(sample);
        }
    }
    
    return genotypes;
}


/* *******************
 *       Sorting     *
 * *******************/

uint8_t *get_individual_phenotypes(vcf_file_t *vcf, ped_file_t *ped, int *num_affected, int *num_unaffected) {
    size_t num_samples = get_num_vcf_samples(vcf);
    individual_t **individuals = sort_individuals(vcf, ped);
    uint8_t *phenotypes = malloc (num_samples * sizeof(uint8_t));
    
    *num_affected = *num_unaffected = 0;
    
    // Set phenotypes value
    for (int i = 0; i < num_samples; i++) {
        if (individuals[i]->condition == AFFECTED) {
            phenotypes[i] = 1;
            (*num_affected)++;
        } else {
            phenotypes[i] = 0;
            (*num_unaffected)++;
        }
    }
    
    free(individuals);
    
    return phenotypes;
}
 
int *group_individuals_by_phenotype(uint8_t *phenotypes, int num_affected, int num_unaffected) {
    int num_samples = num_affected + num_unaffected;
    int *destination = malloc (num_samples * sizeof(int));
    
    // Group samples depending on their phenotype (first cases, then controls)
    int affected_idx = 0, unaffected_idx = num_affected;
    for (int i = 0; i < num_samples; i++) {
        if (phenotypes[i]) {
            destination[i] = affected_idx;
            affected_idx++;
        } else {
            destination[i] = unaffected_idx;
            unaffected_idx++;
        }
        LOG_DEBUG_F("Sample #%d -> phenotype = %d, destination = %d\n", i, phenotypes[i], destination[i]);
    }
    
    return destination;
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
