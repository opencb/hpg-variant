#include "tdt_runner.h"

int verbose_tdt = 0;
int permute = 0;

int run_tdt_test(global_options_data_t* global_options_data, gwas_options_data_t* options_data) {
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("batches", 1, options_data->max_batches, read_list);
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", options_data->num_threads, MIN(10, options_data->max_batches) * options_data->batch_size, output_list);

    int ret_code = 0;
    double start, stop, total;
    vcf_file_t *file = vcf_open(global_options_data->vcf_filename);
    ped_file_t *ped_file = ped_open(global_options_data->ped_filename);
    
    LOG_INFO("About to read PED file...\n");
    // Read PED file before doing any proccessing
    ret_code = ped_read(ped_file);
    if (ret_code != 0) {
        return ret_code;
    }
    
//     create_directory(global_options_data->output_directory);
    
//     // Remove all .txt files in folder
//     ret_code = delete_files_by_extension(global_options_data->output_directory, "txt");
//     if (ret_code != 0) {
//         return ret_code;
//     }
//     char *output_directory = global_options_data->output_directory;
//     size_t output_directory_len = strlen(output_directory);
    
    LOG_INFO("About to perform TDT test...\n");

#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            ret_code = vcf_read_batches(read_list, options_data->batch_size, file, 0);

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
            
            // Create chain of filters for the VCF file
            filter_t **filters = NULL;
            int num_filters = 0;
            if (options_data->chain != NULL) {
                filters = sort_filter_chain(options_data->chain, &num_filters);
            }
    
            // TODO Create map to relate the position of individuals in the list of samples defined in the VCF file
            list_t *sample_names = file->samples_names;
            cp_hashtable *sample_ids = cp_hashtable_create(sample_names->length * 2,
                                                            cp_hash_string,
                                                            (cp_compare_fn) strcasecmp);
            list_item_t *item = sample_names->first_p;
            int *index;
            for (int i = 0; i < sample_names->length && item != NULL; i++) {
                index = (int*) malloc (sizeof(int));
                *index = i;
                cp_hashtable_put(sample_ids, ((individual_t*) item->data_p)->id, index);
                
                item = item->next_p;
            }
            item = NULL;
            
            start = omp_get_wtime();

            int i = 0;
            while ((item = list_remove_item(read_list)) != NULL) {
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

                // Write records that passed to a separate file, and query the WS with them as args
                if (passed_records->length > 0) {
                    // Divide the list of passed records in ranges of size defined in config file
                    int max_chunk_size = 1000;  // TODO define dynamically
                    int num_chunks;
                    list_item_t **chunk_starts = create_chunks(passed_records, max_chunk_size, &num_chunks);
                    
                    // OpenMP: Launch a thread for each range
                    #pragma omp parallel for
                    for (int j = 0; j < num_chunks; j++) {
                        LOG_DEBUG_F("[%d] Test execution\n", omp_get_thread_num());
                        ret_code = tdt_test(ped_file, chunk_starts[j], max_chunk_size, sample_ids);
                    }
                    free(chunk_starts);
                    
                    LOG_INFO_F("*** %dth TDT execution finished\n", i);
                    
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
                    list_free(passed_records, NULL);
                }
                if (failed_records) {
                    LOG_DEBUG_F("[Batch %d] %zu failed records\n", i, failed_records->length);
                    list_free(failed_records, NULL);
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
            if (passed_file) { fclose(passed_file); }
            if (failed_file) { fclose(failed_file); }
            
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
            // TODO write each line of the result to the output file
        }
    }
        
    return ret_code;
}


// tdt_result_t *
int
tdt_test(ped_file_t *ped_file, list_item_t *variants, int num_variants, cp_hashtable *sample_ids) {
    int ret_code = 0;
//     tdt_result_t *result = (tdt_result_t*) calloc (1, sizeof(tdt_result_t));
    cp_hashtable *families = ped_file->families;
    int num_families = get_num_families(ped_file);
    int num_samples = cp_hashtable_count(sample_ids);
    char *sample_data[num_samples];
    
    int father_allele1, father_allele2;
    int mother_allele1, mother_allele2;
    int child_allele1, child_allele2;

    ///////////////////////////////////
    // Perform analysis for each variant

    list_item_t *cur_variant = variants;
    // TODO chunks in the same way as in hpg-variant/effect
    for (int i = 0; i < num_variants && cur_variant != NULL; i++) {
        vcf_record_t *record = (vcf_record_t*) cur_variant->data_p;
//         LOG_INFO_F("Checking variant %s:%ld\n", record->chromosome, record->position);
        
    //         // Adaptive permutation, skip this SNP?
    //         if (par::adaptive_perm && (!perm.snp_test[variant])) {
    //             continue;
    //         }

        list_t *item_samples = record->samples;
        list_item_t *sample = item_samples->first_p;
        // TODO implement arraylist in order to avoid this conversion
        for (int j = 0; j < num_samples && sample != NULL; j++) {
            sample_data[j] = sample->data_p;
            sample = sample->next_p;
        }
    
        // Transmission counts
        
        double t1 = 0;
        double t2 = 0;
        
        // Count over families
        char **keys = (char**) cp_hashtable_get_keys(families);
        family_t *family;
        for (int f = 0; f < num_families; f++) {
            family = cp_hashtable_get(families, keys[f]);
            
            LOG_INFO_F("Checking suitability of family %s\n", family->id);
//             if ( !family[f]->TDT ) continue;

            int trA = 0;  // transmitted allele from first het parent
            int unA = 0;  // untransmitted allele from first het parent
            
            int trB = 0;  // transmitted allele from second het parent
            int unB = 0;  // untransmitted allele from second het parent
            
            individual_t *father = family->father;
            individual_t *mother = family->mother;
            cp_list *children = family->children;

            char *father_sample = sample_data[*((int*) cp_hashtable_get(sample_ids, father->id))];
            char *mother_sample = sample_data[*((int*) cp_hashtable_get(sample_ids, mother->id))];
            
            // If any parent's alleles can't be read or is missing, go to next family
            if (get_alleles(father_sample, &father_allele1, &father_allele2) ||
                get_alleles(mother_sample, &mother_allele1, &mother_allele2)) {
                    continue;
            }
            
            // We need two genotyped parents, with at least one het
            if (father_allele1 == father_allele2 && mother_allele1 == mother_allele2) {
                continue;
            }
            
            if ((father_allele1 && !father_allele2) || (mother_allele1 && !mother_allele2)) {
                continue;
            }


            // Consider all offspring in nuclear family
            cp_list_iterator *children_iterator = cp_list_create_iterator(family->children, COLLECTION_LOCK_READ);
            individual_t *child = NULL;
            while ((child = cp_list_iterator_next(children_iterator)) != NULL) {
                // Only consider affected children
                // TODO Accept non-default specification using 0 as unaffected and 1 as affected
                if (child->phenotype != 2.0f) { continue; }
                
                char *child_sample = sample_data[*((int*) cp_hashtable_get(sample_ids, child->id))];
                if (get_alleles(child_sample, &child_allele1, &child_allele2)) {
                    continue;
                }
                
                // Skip if offspring has missing genotype
                if (child_allele1 && !child_allele2) { continue; }
                
                // We've now established: no missing genotypes
                // and at least one heterozygous parent

                // Kid is 00

                if (!child_allele1 && !child_allele2) {
                    if ( ( (!father_allele1) && father_allele2 ) && 
                        ( (!mother_allele1) && mother_allele2 ) )
                    { trA=1; unA=2; trB=1; unB=2; }
                    else 
                    { trA=1; unA=2; } 
                }
                else if ( (!child_allele1) && child_allele2 )  // Kid is 01
                {
                    // het dad
                    if (father_allele1 != father_allele2 )
                    {
                        // het mum
                        if ( mother_allele1 != mother_allele2 )
                    { trA=1; trB=2; unA=2; unB=1; }
                        else if ( !mother_allele1 ) 
                    { trA=2; unA=1; }
                        else { trA=1; unA=2; }
                    }
                    else if ( !father_allele1 ) 
                    {
                        trA=2; unA=1; 
                    }           
                    else
                    {
                        trA=1; unA=2;
                    }
                }
                else // kid is 1/1
                {
                    
                    if ( ( (!father_allele1) && father_allele2 ) && 
                        ( (!mother_allele1) && mother_allele2 ) )
                    { trA=2; unA=1; trB=2; unB=1; }
                    else 
                    { 
                        trA=2; unA=1;
                    }
                }
                
                // We have now populated trA (first transmission) 
                // and possibly trB also 
                
                ////////////////////////////////////////
                // Permutation? 50:50 flip (precomputed)
                
                if (permute) {
//                     if (flipA[f])
//                     {
                    int t = trA;
                    trA = unA;
                    unA = t;
                    
                    t = trB;
                    trB = unB;
                    unB = t;
//                     }
                }
                
                // Increment transmission counts
                if (trA==1) { t1++; }
                if (trB==1) { t1++; }
                if (trA==2) { t2++; }
                if (trB==2) { t2++; }
                
                if (verbose_tdt) {
                    LOG_INFO_F("TDT\t%s %s : %d %d - %d %d - %f %f - %d %d - %d %d - %d %d\n", 
                           record->id, family->id, trA, unA, trB, unB, t1, t2, 
                           father_allele1, father_allele2, mother_allele1, mother_allele2, child_allele1, child_allele2);
                }
            } // next offspring in family
            cp_list_iterator_destroy(children_iterator);
        
        }  // next nuclear family

        /////////////////////////////
        // Finished counting: now compute
        // the statistics
        
        double tdt_chisq, par_chisq, com_chisq;
        tdt_chisq = par_chisq = com_chisq = -1;
        
        // Basic TDT test
        
        if (t1+t2 > 0) {
            tdt_chisq = ((t1-t2)*(t1-t2))/(t1+t2);
        }
        
//         printf("%s:%ld\t%s-%s\t%f\t%f\t%f\t%f\n",//\t%f\n", 
//                record->chromosome, record->position, record->reference, record->alternate, 
//                t1, t2, (float) t1/t2, tdt_chisq);//p_value);
        
        cur_variant = cur_variant->next_p;
    } // next variant

//     return result;
    return ret_code;
}
