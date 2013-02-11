/*
 * Copyright (c) 2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2013 Ignacio Medina (ICM-CIPF)
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

#include "epistasis_runner.h"

int run_epistasis(shared_options_data_t* shared_options_data, epistasis_options_data_t* options_data) {
    int ret_code = 0;
    
//     list_t *output_list = (list_t*) malloc (sizeof(list_t));
//     list_init("output", shared_options_data->num_threads, INT_MAX, output_list);
// 
//     int ret_code = 0;
//     vcf_file_t *file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
//     if (!file) {
//         LOG_FATAL("VCF file does not exist!\n");
//     }
//     
//     ped_file_t *ped_file = ped_open(shared_options_data->ped_filename);
//     if (!ped_file) {
//         LOG_FATAL("PED file does not exist!\n");
//     }
//     
//     LOG_INFO("About to read PED file...\n");
//     // Read PED file before doing any proccessing
//     ret_code = ped_read(ped_file);
//     if (ret_code != 0) {
//         LOG_FATAL_F("Can't read PED file: %s\n", ped_file->filename);
//     }
    
    // TODO load binary input dataset
    int num_affected, num_unaffected;
    size_t num_variants, file_len, genotypes_offset;
    
    uint8_t *input_file = epistasis_dataset_load(&num_variants, &num_affected, &num_unaffected, &file_len, &genotypes_offset, options_data->dataset_filename);
    uint8_t *genotypes = input_file + genotypes_offset;
    
    // Try to create the directory where the output files will be stored
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }
    
    /*************** Precalculate the rest of variables the algorithm needs  ***************/
    
    int num_samples = num_affected + num_unaffected;
    int num_blocks_per_dim = ceil((double) num_variants / options_data->stride);
    
    // Precalculate combinations for a given options_data->order
    int num_genotype_combinations;
    uint8_t **genotype_combinations = get_genotype_combinations(options_data->order, &num_genotype_combinations);
    
    // Ranking of best models in each repetition
    linked_list_t *best_models[options_data->num_cv_repetitions];
    
    // TODO Not used in basic MDR, so not initialized
//     array_list_t* aux_ret = array_list_new(4, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    
    /**************************** End of variables precalculus  ****************************/
    
    
    for (int r = 0; r < options_data->num_cv_repetitions; r++) {
        // Initialize folds, first block coordinates, genotype combinations and rankings for each repetition
        unsigned int *sizes, *training_sizes;
        int **folds = get_k_folds(num_affected, num_unaffected, options_data->num_folds, &sizes);
        int block_coords[options_data->order]; memset(block_coords, 0, options_data->order * sizeof(int));
        
        linked_list_t *ranking_risky[options_data->num_folds];
        for (int i = 0; i < options_data->num_folds; i++) {
            ranking_risky[i] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
        }
        
        // Calculate size of training datasets
        training_sizes = calloc(3 * options_data->num_folds, sizeof(unsigned int));
        for (int i = 0; i < options_data->num_folds; i++) {
            for (int j = 0; j < options_data->num_folds; j++) {
                if (i != j) {
                    training_sizes[3 * i] += sizes[3 * i];
                    training_sizes[3 * i + 1] += sizes[3 * i + 1];
                    training_sizes[3 * i + 2] += sizes[3 * i + 2];
                }
            }
        }
        
        enum eval_mode mode = TRAINING;
    
        do {
//         printf("BLOCK %d,%d\n", block_coords[0], block_coords[1]);
        
        uint8_t *block_starts[2];
        block_starts[0] = genotypes + block_coords[0] * options_data->stride * num_samples;
        block_starts[1] = genotypes + block_coords[1] * options_data->stride * num_samples;
        
        // Test first combination in the block
        int *comb = get_first_combination_in_block(options_data->order, block_coords, options_data->stride);
//         print_combination(comb, comb_idx, options_data->order);
        
        do  {
            // Run for each fold
            for (int i = 0; i < options_data->num_folds; i++) {
//                 printf("Combination (%d,%d) and fold %d\n", comb[0], comb[1], i);
                
                // Get genotypes of that combination
                uint8_t *training_genotypes = get_genotypes_for_combination_exclude_fold(options_data->order, comb, num_samples, sizes[3 * i], folds[i], options_data->stride, block_starts);
                
                risky_combination *risky_comb = get_model_from_combination_in_fold(options_data->order, comb, training_genotypes,
                                                                                   training_sizes[3 * i + 1], training_sizes[3 * i + 2],
                                                                                   num_genotype_combinations, genotype_combinations, NULL); // TODO aux_ret needed?
                
                if (risky_comb) {
                    // Check the model against the testing dataset
                    double accuracy = 0.0f;
                    
                    if (mode == TESTING) {
                        uint8_t *testing_genotypes = get_genotypes_for_combination_and_fold(options_data->order, risky_comb->combination, 
                                                                                            num_samples, sizes[3 * i + 1] + sizes[3 * i + 2], 
                                                                                            folds[i], options_data->stride, block_starts);
                        accuracy = test_model(options_data->order, risky_comb, testing_genotypes, sizes[3 * i + 1], sizes[3 * i + 2]);
                        free(testing_genotypes);
                    } else {
                        accuracy = test_model(options_data->order, risky_comb, training_genotypes, training_sizes[3 * i + 1], training_sizes[3 * i + 2]);
                    }
//                     printf("*  Balanced accuracy: %.3f\n", accuracy);
                    
                    int position = add_to_model_ranking(risky_comb, options_data->max_ranking_size, ranking_risky[i]);
//                     if (position >= 0) {
//                         printf("Combination inserted at position %d\n", position);
//                     } else {
//                         printf("Combination not inserted\n");
//                     }
                }
                
                free(training_genotypes);
            }
            
        } while (get_next_combination_in_block(options_data->order, comb, block_coords, options_data->stride)); // Test next combinations
        
        free(comb);
        
        } while (get_next_block(num_blocks_per_dim, options_data->order, block_coords));
        
        
        // Merge all rankings in one
        size_t repetition_ranking_size = 0;
        for (int i = 0; i < options_data->num_folds; i++) {
            repetition_ranking_size += ranking_risky[i]->size;
        }
        risky_combination *repetition_ranking[repetition_ranking_size];
        size_t current_index = 0;
        for (int i = 0; i < options_data->num_folds; i++) {
            size_t current_ranking_size = ranking_risky[i]->size;
            linked_list_iterator_t* iter = linked_list_iterator_new(ranking_risky[i]);
            
            risky_combination *element = NULL;
            while(element = linked_list_iterator_next(iter)) {
                repetition_ranking[current_index] = element;
                current_index++;
//                 printf("(%d %d - %.3f) ", element->combination[0], element->combination[1], element->accuracy);
            }
        }
        
        // qsort by coordinates
        qsort(repetition_ranking, repetition_ranking_size, sizeof(risky_combination*), compare_risky);
        
        // Sum all values of each position and get the mean of accuracies
        linked_list_t *sorted_repetition_ranking = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
        risky_combination *current = repetition_ranking[0];
        
        for (int i = 1; i < repetition_ranking_size; i++) {
            risky_combination *element = repetition_ranking[i];
            if (!compare_risky(&current, &element)) {
                current->accuracy += element->accuracy;
            } else {
                current->accuracy /= options_data->num_folds;
                int position = add_to_model_ranking(current, repetition_ranking_size, sorted_repetition_ranking);
                current = element;
            }
        }
        
        // Save the models ranking
        best_models[r] = sorted_repetition_ranking;
        
        free(folds);
        printf("\n\n");
    }
    
    
    
//     // Show the best model of each repetition
//     for (int r = 0; r < options_data->num_cv_repetitions; r++) {
//         risky_combination *element = linked_list_get(0, best_models[r]);
//         assert(element);
//         assert(element->combination);
//         printf("CV %d\t(%d %d - %.3f)\n", r, element->combination[0], element->combination[1], element->accuracy);
//     }
    
    // CVC (get the model that appears more times in the first ranking position)
    khash_t(cvc) *models_for_cvc = kh_init(cvc);
    for (int r = 0; r < options_data->num_cv_repetitions; r++) {
        risky_combination *element = linked_list_get(0, best_models[r]);
        
        int key = 0;
        for (int i = 0; i < options_data->order; i++) {
            key += element->combination[i] * pow(num_samples, options_data->order - i);
        }
        
        int ret;
        khiter_t iter = kh_get(cvc, models_for_cvc, key);
        if (iter != kh_end(models_for_cvc)) {
            (kh_value(models_for_cvc, iter))++; // Increment number of occurrences
        } else {
            iter = kh_put(cvc, models_for_cvc, key, &ret);
            if (ret) {
                kh_value(models_for_cvc, iter) = 1;
            }
        }
        linked_list_free(best_models[r], risky_combination_free);
    }
    
    int bestkey, bestvalue = 0;
    for (int k = kh_begin(models_for_cvc); k < kh_end(models_for_cvc); k++) {
        if (kh_exist(models_for_cvc, k)) {
            int key = kh_key(models_for_cvc, k);
            int value = kh_value(models_for_cvc, k);
            printf("%d -> %d\n", key, value);
            if (value > bestvalue) {
                bestkey = key;
                bestvalue = value;
            }
        }
    }
    
    int bestcomb[options_data->order];
    for (int i = 0; i < options_data->order; i++) {
        bestcomb[i] = bestkey / pow(num_samples, options_data->order - i);
        bestkey -= bestcomb[i] * pow(num_samples, options_data->order - i);
    }
    
    printf("Best model is (%d %d) with a CVC of %d/%d\n", bestcomb[0], bestcomb[1], bestvalue, options_data->num_cv_repetitions);
    kh_destroy(cvc, models_for_cvc);
    
    epistasis_dataset_close(input_file, file_len);
    
    return ret_code;
}


/* *******************
 * Output generation *
 * *******************/

FILE *get_output_file(char *output_directory, char **path) {
    char *filename = "hpg-variant.epi";
    *path = (char*) calloc ((strlen(output_directory) + strlen(filename) + 2), sizeof(char));
    sprintf(*path, "%s/%s", output_directory, filename);
    return fopen(*path, "w");
}

void write_output_header(FILE *fd) {
    assert(fd);
    fprintf(fd, "#CHR         POS       A1      A2         T       U           OR           CHISQ         P-VALUE\n");
}

void write_output_body(list_t* output_list, FILE *fd) {
    assert(fd);
    list_item_t* item = NULL;
//     while (item = list_remove_item(output_list)) {
//         tdt_result_t *result = item->data_p;
//         
//         fprintf(fd, "%s\t%8ld\t%s\t%s\t%3d\t%3d\t%6f\t%6f\t%6f\n",
//                 result->chromosome, result->position, result->reference, result->alternate, 
//                 result->t1, result->t2, result->odds_ratio, result->chi_square, result->p_value);
//         
//         tdt_result_free(result);
//         list_item_free(item);
//     }
}

/* ******************************
 *      Auxiliary functions     *
 * ******************************/

int compare_risky(const void *risky_1, const void *risky_2) {
    risky_combination *r1 = *((risky_combination**) risky_1);
    risky_combination *r2 = *((risky_combination**) risky_2);
    
    for (int i = 0; i < r1->order; i++) {
        if (r1->combination[i] < r2->combination[i]) {
            return -1;
        }
        if (r1->combination[i] > r2->combination[i]) {
            return 1;
        }
    }
    
    return 0;
}
