#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <check.h>

#include <containers/array_list.h>

#include "epistasis/cross_validation.h"
#include "epistasis/mdr.h"
#include "epistasis/model.h"


Suite *create_test_suite(void);


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/



/* ******************************
 *          Unit tests          *
 * ******************************/


START_TEST (test_get_high_risk_combinations) {
    int num_counts = 6, num_risky = 0;
    int counts[] = { 8, 40, 4, 75, 9, 20, 8, 63 } ;
    int num_affected = 10, num_unaffected = 80;
    array_list_t* aux_ret = array_list_new(4, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    
    int *combinations = get_high_risk_combinations(counts, num_counts, num_affected, num_unaffected, 
                                                   &num_risky, aux_ret, mdr_high_risk_combinations);
    
    fail_if(num_risky != 2, "There should be 2 risky combination");
    fail_if(combinations[0] != 0, "Combination 0 (zero) should be risky");
    fail_if(combinations[1] != 2, "Combination 2 should be risky");
}
END_TEST


START_TEST (test_mdr_steps_2_5) {
    const int max_ranking_size = 10;
    const int order = 2;
    const int num_variants = 9, stride = 3, num_blocks = 3;
    const int num_affected = 6, num_unaffected = 6;
    
    const int num_samples = 12, num_folds = 3, num_samples_in_fold = num_samples / num_folds;
    const int aff_per_fold = num_affected / num_folds, aff_per_training_fold = 2 * aff_per_fold;
    const int unaff_per_fold = num_unaffected / num_folds, unaff_per_training_fold = 2 * unaff_per_fold;
    
    int folds[3][4] = { { 0, 1, 6, 7 }, { 4, 5, 10, 11 }, { 2, 3, 8, 9 } };
    unsigned int *sizes;
    
    uint8_t genotypes[] = { 0, 1, 1, 0, 1, 2, 2, 1, 2, 1, 1, 2,   // First block
                            1, 2, 2, 1, 2, 0, 0, 2, 0, 2, 2, 0,
                            2, 0, 0, 2, 0, 1, 1, 0, 1, 0, 0, 1,
                            0, 1, 0, 1, 1, 2, 1, 2, 2, 0, 2, 0,   // Second block
                            1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1,
                            1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 1, 2,
                            1, 0, 1, 2, 0, 2, 0, 2, 0, 2, 2, 1,   // Third block
                            1, 2, 2, 1, 1, 0, 1, 0, 1, 2, 0, 2,
                            2, 0, 0, 0, 1, 2, 0, 2, 2, 1, 1, 0 };
    
    array_list_t* aux_ret = array_list_new(4, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    
    int block_2d[] = { 0, 0 };
    int block_3d[] = { 0, 0, 0 };
    
    // Precalculate combinations for a given order
    int num_genotype_combinations;
    uint8_t **genotype_combinations = get_genotype_combinations(order, &num_genotype_combinations);
    
    array_list_t *ranking_risky[num_folds];
    for (int i = 0; i < num_folds; i++) {
        ranking_risky[i] = array_list_new(100, 1.5, COLLECTION_MODE_ASYNCHRONIZED);
    }
    
    do {
        printf("BLOCK %d,%d\n", block_2d[0], block_2d[1]);
//         printf("BLOCK %d,%d,%d\n", block_3d[0], block_3d[1], block_3d[2]);
        
        uint8_t *block_starts[2];
        block_starts[0] = genotypes + block_2d[0] * stride * num_samples;
        block_starts[1] = genotypes + block_2d[1] * stride * num_samples;
        
        // Test first combination in the block
        int *comb = get_first_combination_in_block(order, block_2d, stride);
//         print_combination(comb, comb_idx, order);
        
        do  {
            // Run for each fold
            for (int i = 0; i < num_folds; i++) {
                printf("Combination (%d,%d) and fold %d\n", comb[0], comb[1], i);
                
                risky_combination *risky_comb = get_model_from_combination_in_fold(order, comb, folds[i], num_samples_in_fold, 
                                                                                   num_samples, aff_per_training_fold, unaff_per_training_fold,
                                                                                   stride, block_starts, num_genotype_combinations, genotype_combinations,
                                                                                   aux_ret);
                
                if (risky_comb) {
                    // Check the model against the testing dataset
                    double accuracy = test_model(order, risky_comb, folds[i], num_samples, aff_per_fold, unaff_per_fold, stride, block_starts);
                    printf("*  Balanced accuracy: %.3f\n", accuracy);
                    
                    int position = add_to_model_ranking(risky_comb, max_ranking_size, ranking_risky[i]);
                    if (position >= 0) {
                        printf("Combination inserted at position %d\n", position);
                    } else {
                        printf("Combination not inserted\n");
                    }
                }
            }
            
            printf("\n-------------\n");
            
        } while (get_next_combination_in_block(order, comb, block_2d, stride)); // Test next combinations
        
        free(comb);
        
        
        printf("\n==============\n");
    } while (get_next_block(num_blocks, order, block_2d));
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
    tcase_add_test(tc_cv, test_mdr_steps_2_5);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("MDR");
    suite_add_tcase(fs, tc_counts);
    suite_add_tcase(fs, tc_cv);
    
    return fs;
}

