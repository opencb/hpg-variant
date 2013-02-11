#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <check.h>
#include <khash.h>

#include <containers/array_list.h>
#include <containers/linked_list.h>

#include "epistasis/cross_validation.h"
#include "epistasis/mdr.h"
#include "epistasis/model.h"


Suite *create_test_suite(void);

int compare_risky(const void *risky_1, const void *risky_2);
int compare_accuracies(const void *risky_1, const void *risky_2);


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/

KHASH_MAP_INIT_INT(cvc, int);

/* ******************************
 *          Unit tests          *
 * ******************************/


START_TEST (test_get_high_risk_combinations) {
    int num_counts = 6, num_risky = 0;
    int counts[] = { 8, 40, 4, 75, 9, 20, 8, 63 } ;
    int num_affected = 10, num_unaffected = 80;
    
    int *combinations = get_high_risk_combinations(counts, num_counts, num_affected, num_unaffected, 
                                                   &num_risky, NULL, mdr_high_risk_combinations);
    
    fail_if(num_risky != 2, "There should be 2 risky combination");
    fail_if(combinations[0] != 0, "Combination 0 (zero) should be risky");
    fail_if(combinations[1] != 2, "Combination 2 should be risky");
}
END_TEST


START_TEST (test_mdr_steps_2_6) {
    const int max_ranking_size = 10;
    const int order = 2, stride = 3, num_blocks = 3;
    const int num_affected = 6, num_unaffected = 6;
    
    const int num_samples = 12, num_folds = 3, num_samples_in_fold = num_samples / num_folds;
    const int aff_per_fold = num_affected / num_folds, aff_per_training_fold = 2 * aff_per_fold;
    const int unaff_per_fold = num_unaffected / num_folds, unaff_per_training_fold = 2 * unaff_per_fold;
    
    int folds[3][4] = { { 0, 1, 6, 7 }, { 4, 5, 10, 11 }, { 2, 3, 8, 9 } };
    uint8_t genotypes[] = { 0, 1, 1, 0, 1, 2, 2, 1, 2, 1, 1, 2,   // First block
                            1, 2, 2, 1, 2, 0, 0, 2, 0, 2, 2, 0,
                            2, 0, 0, 2, 0, 1, 1, 0, 1, 0, 0, 1,
                            0, 1, 0, 1, 1, 2, 1, 2, 2, 0, 2, 0,   // Second block
                            1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1,
                            1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 1, 2,
                            1, 0, 1, 2, 0, 2, 0, 2, 0, 2, 2, 1,   // Third block
                            1, 2, 2, 1, 1, 0, 1, 0, 1, 2, 0, 2,
                            2, 0, 0, 0, 1, 2, 0, 2, 2, 1, 1, 0 };
    
    int block_2d[] = { 0, 0 };
    
    // Precalculate combinations for a given order
    int num_genotype_combinations;
    uint8_t **genotype_combinations = get_genotype_combinations(order, &num_genotype_combinations);
    
    linked_list_t *ranking_risky[num_folds];
    for (int i = 0; i < num_folds; i++) {
        ranking_risky[i] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    }
    
    do {
//         printf("BLOCK %d,%d\n", block_2d[0], block_2d[1]);
        
        uint8_t *block_starts[2];
        block_starts[0] = genotypes + block_2d[0] * stride * num_samples;
        block_starts[1] = genotypes + block_2d[1] * stride * num_samples;
        
        // Test first combination in the block
        int *comb = get_first_combination_in_block(order, block_2d, stride);
//         print_combination(comb, comb_idx, order);
        
        do  {
            // Run for each fold
            for (int i = 0; i < num_folds; i++) {
//                 printf("Combination (%d,%d) and fold %d\n", comb[0], comb[1], i);
                
                uint8_t *training_genotypes = get_genotypes_for_combination_exclude_fold(order, comb, num_samples, num_samples_in_fold, folds[i], stride, block_starts);
                
                risky_combination *risky_comb = get_model_from_combination_in_fold(order, comb, training_genotypes, aff_per_training_fold, unaff_per_training_fold, 
                                                                                   num_genotype_combinations, genotype_combinations);
                
                if (risky_comb) {
                    // Check the model against the testing dataset
                    double accuracy = test_model(order, risky_comb, training_genotypes, aff_per_training_fold, unaff_per_training_fold);
//                     printf("*  Balanced accuracy: %.3f\n", accuracy);
                    
                    int position = add_to_model_ranking(risky_comb, max_ranking_size, ranking_risky[i]);
//                     if (position >= 0) {
//                         printf("Combination inserted at position %d\n", position);
//                     } else {
//                         printf("Combination not inserted\n");
//                     }
                }
            }
            
//             printf("-------------\n");
            
        } while (get_next_combination_in_block(order, comb, block_2d, stride)); // Test next combinations
        
        free(comb);
        
        
//         printf("\n==============\n");
    } while (get_next_block(num_blocks, order, block_2d));
}
END_TEST


START_TEST(test_mdr_5_repetitions_5_fold) {
    const int num_repetitions = 5, num_variants = 10, num_folds = 5, order = 2, stride = 4, num_blocks = ceil((double) num_variants / stride);
    const int max_ranking_size = 20;
    const unsigned int num_samples = 40, num_affected = 20, num_unaffected = 20;
    
    const uint8_t genotypes[10 * 40] = { 
        1, 1, 0, 1, 0, 1, 2, 1, 0, 0, 2, 1, 1, 1, 2, 1, 0, 2, 0, 1, 1, 0, 2, 2, 1, 1, 1, 1, 0, 0, 1, 2, 0, 2, 1, 0, 0, 0, 1, 1,
        1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 2, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 2, 0, 1, 0, 1,
        0, 0, 2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 2, 0, 0, 0, 1,
        0, 0, 2, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 0, 0, 0, 1,
        0, 0, 2, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 0, 0, 0, 1,
        0, 0, 2, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 0, 0, 0, 1 };
    
    // Precalculate combinations for a given order
    int num_genotype_combinations;
    uint8_t **genotype_combinations = get_genotype_combinations(order, &num_genotype_combinations);
    
    // Ranking of best models in each repetition
    linked_list_t *best_models[num_repetitions];
    
    for (int r = 0; r < num_repetitions; r++) {
        // Initialize folds, first block coordinates, genotype combinations and rankings for each repetition
        unsigned int *sizes, *training_sizes;
        int **folds = get_k_folds(num_affected, num_unaffected, num_folds, &sizes);
        int block_coords[order]; memset(block_coords, 0, order * sizeof(int));
        
        linked_list_t *ranking_risky[num_folds];
        for (int i = 0; i < num_folds; i++) {
            ranking_risky[i] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
        }
        
        // Calculate size of training datasets
        training_sizes = calloc(3 * num_folds, sizeof(unsigned int));
        for (int i = 0; i < num_folds; i++) {
            for (int j = 0; j < num_folds; j++) {
                if (i != j) {
                    training_sizes[3 * i] += sizes[3 * i];
                    training_sizes[3 * i + 1] += sizes[3 * i + 1];
                    training_sizes[3 * i + 2] += sizes[3 * i + 2];
                }
            }
        }
        
        enum eval_mode mode = TRAINING;
    
        // TODO Run MDR steps 2-6
        do {
//         printf("BLOCK %d,%d\n", block_coords[0], block_coords[1]);
        
        uint8_t *block_starts[2];
        block_starts[0] = genotypes + block_coords[0] * stride * num_samples;
        block_starts[1] = genotypes + block_coords[1] * stride * num_samples;
        
        // Test first combination in the block
        int *comb = get_first_combination_in_block(order, block_coords, stride);
//         print_combination(comb, comb_idx, order);
        
        do  {
            // Run for each fold
            for (int i = 0; i < num_folds; i++) {
//                 printf("Combination (%d,%d) and fold %d\n", comb[0], comb[1], i);
                
                // Get genotypes of that combination
                uint8_t *training_genotypes = get_genotypes_for_combination_exclude_fold(order, comb, num_samples, sizes[3 * i], folds[i], stride, block_starts);
                
                risky_combination *risky_comb = get_model_from_combination_in_fold(order, comb, training_genotypes,
                                                                                   training_sizes[3 * i + 1], training_sizes[3 * i + 2],
                                                                                   num_genotype_combinations, genotype_combinations);
                
                if (risky_comb) {
                    // Check the model against the testing dataset
                    double accuracy = 0.0f;
                    
                    if (mode == TESTING) {
                        uint8_t *testing_genotypes = get_genotypes_for_combination_and_fold(order, risky_comb->combination, 
                                                                                            num_samples, sizes[3 * i + 1] + sizes[3 * i + 2], 
                                                                                            folds[i], stride, block_starts);
//                         accuracy = test_model(order, risky_comb, folds[i], num_samples, sizes[3 * i + 1], sizes[3 * i + 2], stride, block_starts);
                        accuracy = test_model(order, risky_comb, testing_genotypes, sizes[3 * i + 1], sizes[3 * i + 2]);
                        free(testing_genotypes);
                    } else {
                        accuracy = test_model(order, risky_comb, training_genotypes, training_sizes[3 * i + 1], training_sizes[3 * i + 2]);
                    }
//                     printf("*  Balanced accuracy: %.3f\n", accuracy);
                    
                    int position = add_to_model_ranking(risky_comb, max_ranking_size, ranking_risky[i]);
//                     if (position >= 0) {
//                         printf("Combination inserted at position %d\n", position);
//                     } else {
//                         printf("Combination not inserted\n");
//                     }
                }
                
                free(training_genotypes);
            }
            
//             printf("-------------\n");
            
        } while (get_next_combination_in_block(order, comb, block_coords, stride)); // Test next combinations
        
        free(comb);
        
        
//         printf("\n==============\n");
        } while (get_next_block(num_blocks, order, block_coords));
            
        // Print rankings
        size_t repetition_ranking_size = 0;
        for (int i = 0; i < num_folds; i++) {
//             size_t current_ranking_size = ranking_risky[i]->size;
//             linked_list_iterator_t* iter = linked_list_iterator_new(ranking_risky[i]);
//             risky_combination *last_element = current_ranking_size > 0 ? linked_list_get_last(ranking_risky[i]) : NULL;
//             risky_combination *element = NULL;
//             
//             printf("Ranking %d (size %zu) = { ", i, current_ranking_size);
//             while(element = linked_list_iterator_next(iter)) {
//                 printf("(%d %d - %.3f) ", element->combination[0], element->combination[1], element->accuracy);
//             }
//             printf("}\n");
            
            repetition_ranking_size += ranking_risky[i]->size;
        }
    //     printf("\n");
        
        // Merge all rankings in one
        risky_combination *repetition_ranking[repetition_ranking_size];
        size_t current_index = 0;
        for (int i = 0; i < num_folds; i++) {
            size_t current_ranking_size = ranking_risky[i]->size;
            linked_list_iterator_t* iter = linked_list_iterator_new(ranking_risky[i]);
            
            risky_combination *element = NULL;
            while(element = linked_list_iterator_next(iter)) {
                repetition_ranking[current_index] = element;
                current_index++;
    //             printf("(%d %d - %.3f) ", element->combination[0], element->combination[1], element->accuracy);
            }
        }
    //     printf("\n\n");
        
        // qsort by coordinates
        qsort(repetition_ranking, repetition_ranking_size, sizeof(risky_combination*), compare_risky);
    //     for (int i = 0; i < repetition_ranking_size; i++) {
    //         risky_combination *element = repetition_ranking[i];
    //         printf("(%d %d - %.3f) ", element->combination[0], element->combination[1], element->accuracy);
    //     }
    //     printf("\n\n");
        
        linked_list_t *sorted_repetition_ranking = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
        
        // Sum all values of each position and get the mean of accuracies
        risky_combination *current = repetition_ranking[0];
        
        for (int i = 1; i < repetition_ranking_size; i++) {
            risky_combination *element = repetition_ranking[i];
            if (!compare_risky(&current, &element)) {
                current->accuracy += element->accuracy;
            } else {
                current->accuracy /= num_folds;
                int position = add_to_model_ranking(current, repetition_ranking_size, sorted_repetition_ranking);
                current = element;
            }
        }
        
//         linked_list_iterator_t* iter = linked_list_iterator_new(sorted_repetition_ranking);
//         risky_combination *element = NULL;
//         while(element = linked_list_iterator_next(iter)) {
//             printf("(%d %d - %.3f) ", element->combination[0], element->combination[1], element->accuracy);
//         }
//         linked_list_iterator_free(iter);
        
        best_models[r] = sorted_repetition_ranking;
        
        free(folds);
        printf("\n\n");
    }
    
//     // Show the best model of each repetition
//     for (int r = 0; r < num_repetitions; r++) {
//         risky_combination *element = linked_list_get(0, best_models[r]);
//         assert(element);
//         assert(element->combination);
//         printf("CV %d\t(%d %d - %.3f)\n", r, element->combination[0], element->combination[1], element->accuracy);
//     }
    
    // CVC (get the model that appears more times in the first position)
    khash_t(cvc) *models_for_cvc = kh_init(cvc);
    for (int r = 0; r < num_repetitions; r++) {
        risky_combination *element = linked_list_get(0, best_models[r]);
        
        int key = 0;
        for (int i = 0; i < order; i++) {
            key += element->combination[i] * pow(num_samples, order - i);
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
    
    int bestcomb[order];
    for (int i = 0; i < order; i++) {
        bestcomb[i] = bestkey / pow(num_samples, order - i);
        bestkey -= bestcomb[i] * pow(num_samples, order - i);
    }
    
    printf("Best model is (%d %d) with a CVC of %d/%d\n", bestcomb[0], bestcomb[1], bestvalue, num_repetitions);
    kh_destroy(cvc, models_for_cvc);
}
END_TEST


/* ******************************
 *      Main entry point        *
 * ******************************/

int main (int argc, char *argv) {
    Suite *fs = create_test_suite();
    SRunner *fs_runner = srunner_create(fs);
    srunner_run_all(fs_runner, CK_NORMAL);
    int number_failed = srunner_ntests_failed (fs_runner);
    srunner_free (fs_runner);
    
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


Suite *create_test_suite(void) {
    TCase *tc_counts = tcase_create("Risk classification");
    tcase_add_test(tc_counts, test_get_high_risk_combinations);
    
    TCase *tc_cv = tcase_create("Cross validation process");
    tcase_add_test(tc_cv, test_mdr_steps_2_6);
    tcase_add_test(tc_cv, test_mdr_5_repetitions_5_fold);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("MDR");
    suite_add_tcase(fs, tc_counts);
    suite_add_tcase(fs, tc_cv);
    
    return fs;
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

int compare_accuracies(const void *risky_1, const void *risky_2) {
    risky_combination *r1 = *((risky_combination**) risky_1);
    risky_combination *r2 = *((risky_combination**) risky_2);
    
    return (r1->accuracy - r2->accuracy) * 100;
}
