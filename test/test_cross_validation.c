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

void check_copied_genotypes_fold_0(uint8_t *genotypes, uint8_t *block_starts[2], int snp_offset[2], int num_samples);

void check_copied_genotypes_fold_1(uint8_t *genotypes, uint8_t *block_starts[2], int snp_offset[2], int num_samples);

void check_copied_genotypes_fold_2(uint8_t *genotypes, uint8_t *block_starts[2], int snp_offset[2], int num_samples);


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/



/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST (test_get_k_folds) {
    int **folds;
    unsigned int *sizes;
    unsigned int k;
    int i;
    
    // 400 { 200 aff, 200 unaff} in 10-fold
    k = 10;
    folds = get_k_folds(200, 200, k, &sizes);
    for (i = 0; i < k; i++) {
        fail_unless (sizes[3 * i] == 40 && sizes[3 * i + 1] == 20 && sizes[3 * i + 2] == 20, "Fold size must be (40,20,20)");
        
        for (int j = 0; j < sizes[3 * i + 1]; j++) {
            fail_if(folds[i][j] >= 200, "(200,200,10) Affected samples can't have the index of unaffected ones");
        }
        for (int j = 0; j < sizes[3 * i + 2]; j++) {
            fail_if(folds[i][sizes[3 * i + 1] + j] < 200, "(200,200,10) Unaffected samples can't have the index of affected ones");
        }
    }
    
    // 400 { 150 aff, 250 unaff} in 4-fold
    k = 4;
    folds = get_k_folds(150, 250, 4, &sizes);
    i = 0; fail_unless (sizes[3 * i] == 101 && sizes[3 * i + 1] == 38 && sizes[3 * i + 2] == 63, "Fold #0 size must be (101,38,63)");
    i = 1; fail_unless (sizes[3 * i] == 101 && sizes[3 * i + 1] == 38 && sizes[3 * i + 2] == 63, "Fold #1 size must be (101,38,63)");
    i = 2; fail_unless (sizes[3 * i] == 99 && sizes[3 * i + 1] == 37 && sizes[3 * i + 2] == 62, "Fold #2 size must be (99,37,62)");
    i = 3; fail_unless (sizes[3 * i] == 99 && sizes[3 * i + 1] == 37 && sizes[3 * i + 2] == 62, "Fold #3 size must be (99,37,62)");
    
    for (i = 0; i < k; i++) {
        for (int j = 0; j < sizes[3 * i + 1]; j++) {
            fail_if(folds[i][j] >= 150, "(150,250,4) Affected samples can't have the index of unaffected ones");
        }
        for (int j = 0; j < sizes[3 * i + 2]; j++) {
            fail_if(folds[i][sizes[3 * i + 1] + j] < 150, "(150,250,4) Unaffected samples can't have the index of affected ones");
        }
    }
    
    // 400 { 150 aff, 250 unaff} in 10-fold
    k = 10;
    folds = get_k_folds(150, 250, k, &sizes);
    for (i = 0; i < k; i++) {
        fail_unless (sizes[3 * i] == 40 && sizes[3 * i + 1] == 15 && sizes[3 * i + 2] == 25, "Fold size must be (40,15,25)");
        
        for (int j = 0; j < sizes[3 * i + 1]; j++) {
            fail_if(folds[i][j] >= 150, "(150,250,10) Affected samples can't have the index of unaffected ones");
        }
        for (int j = 0; j < sizes[3 * i + 2]; j++) {
            fail_if(folds[i][sizes[3 * i + 1] + j] < 150, "(150,250,10) Unaffected samples can't have the index of affected ones");
        }
    }
    
    
    // 125 { 50 aff, 75 unaff} in 5-fold
    k = 5;
    folds = get_k_folds(50, 75, k, &sizes);
    for (i = 0; i < k; i++) {
        fail_unless (sizes[3 * i] == 25 && sizes[3 * i + 1] == 10 && sizes[3 * i + 2] == 15, "Fold size must be (25,10,15)");
        
        for (int j = 0; j < sizes[3 * i + 1]; j++) {
            fail_if(folds[i][j] >= 50, "(50,75,5) Affected samples can't have the index of unaffected ones");
        }
        for (int j = 0; j < sizes[3 * i + 2]; j++) {
            fail_if(folds[i][sizes[3 * i + 1] + j] < 50, "(50,75,5) Unaffected samples can't have the index of affected ones");
        }
    }
    
    // 125 { 50 aff, 75 unaff} in 7-fold
    k = 7;
    folds = get_k_folds(50, 75, k, &sizes);
    i = 0; fail_unless (sizes[3 * i] == 19 && sizes[3 * i + 1] == 8 && sizes[3 * i + 2] == 11, "Fold #0 size must be (19,8,11)");
    for (i = 1; i < k-2; i++) {
        fail_unless (sizes[3 * i] == 18 && sizes[3 * i + 1] == 7 && sizes[3 * i + 2] == 11, "Fold size must be (18,7,11)");
    }
    i = k-2; fail_unless (sizes[3 * i] == 17 && sizes[3 * i + 1] == 7 && sizes[3 * i + 2] == 10, "Fold k-2 size must be (17,7,10)");
    i = k-1; fail_unless (sizes[3 * i] == 17 && sizes[3 * i + 1] == 7 && sizes[3 * i + 2] == 10, "Fold k-1 size must be (17,7,10)");
    
    for (i = 0; i < k; i++) {
        for (int j = 0; j < sizes[3 * i + 1]; j++) {
            fail_if(folds[i][j] >= 50, "(50,75,7) Affected samples can't have the index of unaffected ones");
        }
        for (int j = 0; j < sizes[3 * i + 2]; j++) {
            fail_if(folds[i][sizes[3 * i + 1] + j] < 50, "(50,75,7) Unaffected samples can't have the index of affected ones");
        }
    }
    
    // 125 { 50 aff, 75 unaff} in 10-fold
    k = 10;
    folds = get_k_folds(50, 75, k, &sizes);
    for (i = 0; i < 5; i++) {
        fail_unless (sizes[3 * i] == 13 && sizes[3 * i + 1] == 5 && sizes[3 * i + 2] == 8, "Fold size must be (13,5,8)");
    }
    for (i = 5; i < k; i++) {
        fail_unless (sizes[3 * i] == 12 && sizes[3 * i + 1] == 5 && sizes[3 * i + 2] == 7, "Fold size must be (12,5,7)");
    }
    
    for (i = 0; i < k; i++) {
        for (int j = 0; j < sizes[3 * i + 1]; j++) {
            fail_if(folds[i][j] >= 50, "(50,75,10) Affected samples can't have the index of unaffected ones");
        }
        for (int j = 0; j < sizes[3 * i + 2]; j++) {
            fail_if(folds[i][sizes[3 * i + 1] + j] < 50, "(50,75,10) Unaffected samples can't have the index of affected ones");
        }
    }
    
    // 20 { 8 aff, 12 unaff) in 8-fold
    k = 8;
    folds = get_k_folds(8, 12, k, &sizes);
    for (i = 0; i < 4; i++) {
        fail_unless (sizes[3 * i] == 3 && sizes[3 * i + 1] == 1 && sizes[3 * i + 2] == 2, "Fold size must be (3,1,2)");
    }
    for (i = 4; i < k; i++) {
        fail_unless (sizes[3 * i] == 2 && sizes[3 * i + 1] == 1 && sizes[3 * i + 2] == 1, "Fold size must be (2,1,1)");
    }
    
    for (i = 0; i < k; i++) {
        for (int j = 0; j < sizes[3 * i + 1]; j++) {
            fail_if(folds[i][j] >= 8, "(8,12,8) Affected samples can't have the index of unaffected ones");
        }
        for (int j = 0; j < sizes[3 * i + 2]; j++) {
            fail_if(folds[i][sizes[3 * i + 1] + j] < 8, "(8,12,8) Unaffected samples can't have the index of affected ones");
        }
    }
    
    // 20 { 8 aff, 12 unaff) in 10-fold
    k = 10;
    folds = get_k_folds(8, 12, k, &sizes);
    for (i = 0; i < 2; i++) {
        fail_unless (sizes[3 * i] == 3 && sizes[3 * i + 1] == 1 && sizes[3 * i + 2] == 2, "Fold size must be (3,1,2)");
    }
    for (i = 2; i < 8; i++) {
        fail_unless (sizes[3 * i] == 2 && sizes[3 * i + 1] == 1 && sizes[3 * i + 2] == 1, "Fold size must be (2,1,1)");
    }
    for (i = 9; i < k; i++) {
        fail_unless (sizes[3 * i] == 1 && sizes[3 * i + 1] == 0 && sizes[3 * i + 2] == 1, "Fold size must be (1,0,1)");
    }
    
    for (i = 0; i < k; i++) {
        for (int j = 0; j < sizes[3 * i + 1]; j++) {
            fail_if(folds[i][j] >= 8, "(8,12,10) Affected samples can't have the index of unaffected ones");
        }
        for (int j = 0; j < sizes[3 * i + 2]; j++) {
            fail_if(folds[i][sizes[3 * i + 1] + j] < 8, "(8,12,10) Unaffected samples can't have the index of affected ones");
        }
    }
    
    
    // 20 { 16 aff, 4 unaff) in 5-fold
    k = 5;
    folds = get_k_folds(16, 4, k, &sizes);
    i = 0; fail_unless (sizes[3 * i] == 5 && sizes[3 * i + 1] == 4 && sizes[3 * i + 2] == 1, "Fold #0 size must be (5,4,1)");
    i = 1; fail_unless (sizes[3 * i] == 4 && sizes[3 * i + 1] == 3 && sizes[3 * i + 2] == 1, "Fold #1 size must be (4,3,1)");
    i = 2; fail_unless (sizes[3 * i] == 4 && sizes[3 * i + 1] == 3 && sizes[3 * i + 2] == 1, "Fold #2 size must be (4,3,1)");
    i = 3; fail_unless (sizes[3 * i] == 4 && sizes[3 * i + 1] == 3 && sizes[3 * i + 2] == 1, "Fold #3 size must be (4,3,1)");
    i = 4; fail_unless (sizes[3 * i] == 3 && sizes[3 * i + 1] == 3 && sizes[3 * i + 2] == 0, "Fold #4 size must be (3,3,0)");
    
    for (i = 0; i < k; i++) {
        for (int j = 0; j < sizes[3 * i + 1]; j++) {
            fail_if(folds[i][j] >= 16, "(16,4,5) Affected samples can't have the index of unaffected ones");
        }
        for (int j = 0; j < sizes[3 * i + 2]; j++) {
            fail_if(folds[i][sizes[3 * i + 1] + j] < 16, "(16,4,5) Unaffected samples can't have the index of affected ones");
        }
    }
    
    // 20 { 16 aff, 4 unaff) in 10-fold
    k = 10;
    folds = get_k_folds(16, 4, k, &sizes);
    for (i = 0; i < 4; i++) {
        fail_unless (sizes[3 * i] == 3 && sizes[3 * i + 1] == 2 && sizes[3 * i + 2] == 1, "Fold size must be (3,2,1)");
    }
    for (i = 4; i < 6; i++) {
        fail_unless (sizes[3 * i] == 2 && sizes[3 * i + 1] == 2 && sizes[3 * i + 2] == 0, "Fold size must be (2,2,0)");
    }
    for (i = 6; i < k; i++) {
        fail_unless (sizes[3 * i] == 1 && sizes[3 * i + 1] == 1 && sizes[3 * i + 2] == 0, "Fold size must be (1,1,0)");
    }
    
    for (i = 0; i < k; i++) {
        for (int j = 0; j < sizes[3 * i + 1]; j++) {
            fail_if(folds[i][j] >= 16, "(16,4,10) Affected samples can't have the index of unaffected ones");
        }
        for (int j = 0; j < sizes[3 * i + 2]; j++) {
            fail_if(folds[i][sizes[3 * i + 1] + j] < 16, "(16,4,10) Unaffected samples can't have the index of affected ones");
        }
    }
    
}
END_TEST


START_TEST (test_get_genotypes_for_combination_and_fold) {
    int order = 2, num_folds = 3, stride = 3, num_blocks = 3;
    unsigned int *sizes;
    
    int num_samples = 12, num_samples_in_fold = 4;
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
        
//     printf("Genotypes = {\n");
//     for (int i = 0; i < 9; i++) {
//         printf("\t{ ");
//         for (int j = 0; j < 12; j++) {
//             printf("%d ", genotypes[i * 12 + j]);
//         }
//         printf(" }\n");
//     }
//     printf(" }\n\n");
    
    int block_2d[] = { 0, 0 };
    do {
//         printf("BLOCK %d,%d\n", block_2d[0], block_2d[1]);
        uint8_t *block_starts[2];
        block_starts[0] = genotypes + block_2d[0] * stride * num_samples;
        block_starts[1] = genotypes + block_2d[1] * stride * num_samples;
        
        // Test first combination in the block
        int *comb = get_first_combination_in_block(order, block_2d, stride);
        int snp_offset[] = { comb[0] % stride, comb[1] % stride };
        
        // Run for each fold
        for (int f = 0; f < num_folds; f++) {
//             printf("Combination (%d,%d) and fold %d\n", comb[0], comb[1], f);
            
            int next_fold = (f < 2) ? f + 1 : 0 ;
            uint8_t *val = get_genotypes_for_combination_exclude_fold(order, comb, num_samples, num_samples_in_fold, folds[f], stride, block_starts);
            
//             printf("* val = { ");
//             for (int i = 0; i < order; i++) {
//                 for (int j = 0; j < (num_samples - num_samples_in_fold); j++) {
//                     printf("%u ", val[i * (num_samples - num_samples_in_fold) + j]);
//                 }
//                 printf("\t");
//             }
//             printf("}\n\n");
            
            if (f == 0) {
                check_copied_genotypes_fold_0(val, block_starts, snp_offset, num_samples);
            } else if (f == 1) {
                check_copied_genotypes_fold_1(val, block_starts, snp_offset, num_samples);
            } else if (f == 2) {
                check_copied_genotypes_fold_2(val, block_starts, snp_offset, num_samples);
            }
            
            free(val);
        }
        
        // Test next combinations
        while (get_next_combination_in_block(order, comb, block_2d, stride)) {
            int snp_offset[] = { comb[0] % stride, comb[1] % stride };
            
            for (int f = 0; f < num_folds; f++) {
//                 printf("Combination (%d,%d) and fold %d\n", comb[0], comb[1], f);
            
                int next_fold = (f < 2) ? f + 1 : 0 ;
                uint8_t *val = get_genotypes_for_combination_exclude_fold(order, comb, num_samples, num_samples_in_fold, folds[f], stride, block_starts);
                
//                 printf("* val = { ");
//                 for (int j = 0; j < order; j++) {
//                     for (int k = 0; k < (num_samples - num_samples_in_fold); k++) {
//                         printf("%u ", val[j * (num_samples - num_samples_in_fold) + k]);
//                     }
//                     printf("\t");
//                 }
//                 printf("}\n\n");
            
                if (f == 0) {
                    check_copied_genotypes_fold_0(val, block_starts, snp_offset, num_samples);
                } else if (f == 1) {
                    check_copied_genotypes_fold_1(val, block_starts, snp_offset, num_samples);
                } else if (f == 2) {
                    check_copied_genotypes_fold_2(val, block_starts, snp_offset, num_samples);
                }
            
                free(val);
            }
        }
        
        free(comb);
        
//         printf("\n-----------\n");
    } while (get_next_block(num_blocks, order, block_2d));
}
END_TEST


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
    
//     array_list_t **risky_combinations = malloc(num_folds * sizeof(risky_combination*));
    array_list_t *risky_combinations[num_folds];
    for (int i = 0; i < num_folds; i++) {
        risky_combinations[i] = array_list_new(100, 1.5, COLLECTION_MODE_ASYNCHRONIZED);
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
        
        // Run for each fold
        for (int i = 0; i < num_folds; i++) {
            printf("Combination (%d,%d) and fold %d\n", comb[0], comb[1], i);
            
            risky_combination *risky_comb = check_risk_of_combination_in_fold(order, comb, i, folds[i], num_samples_in_fold, 
                                                                              num_samples, aff_per_training_fold, unaff_per_training_fold,
                                                                              stride, block_starts,
                                                                              num_genotype_combinations, genotype_combinations,
                                                                              aux_ret);
            
            if (risky_comb) {
                // TODO step 5 -> check against the testing dataset
                // TODO step 6 -> ellaborate a ranking of best N combinations
                // TODO would be better to filter before inserting into the list of risky combinations?
                array_list_insert(risky_comb, risky_combinations[i]);
            }
        }
        
        // Test next combinations
        while (get_next_combination_in_block(order, comb, block_2d, stride)) {
            for (int i = 0; i < num_folds; i++) {
                printf("Combination (%d,%d) and fold %d\n", comb[0], comb[1], i);
                
                risky_combination *risky_comb = check_risk_of_combination_in_fold(order, comb, i, folds[i], num_samples_in_fold, 
                                                                                num_samples, aff_per_training_fold, unaff_per_training_fold,
                                                                                stride, block_starts,
                                                                                num_genotype_combinations, genotype_combinations,
                                                                                aux_ret);
                
                if (risky_comb) {
                    // Insert into the list of risky combination of that fold
                    array_list_insert(risky_comb, risky_combinations[i]);
                }
            }
        }
        
        free(comb);
        
        
        // TODO step 5 -> check against the testing dataset
        // TODO step 6 -> ellaborate a ranking of best N combinations
        // TODO would be better to filter before inserting into the list of risky combinations?
        for (int i = 0; i < num_folds; i++) {
            printf("Checking against fold %d\n", i);
            
            for (int j = 0; j < array_list_size(risky_combinations[i]); j++) {
                risky_combination *risky_comb = array_list_get(j, risky_combinations[i]);
                uint8_t *val = get_genotypes_for_combination_and_fold(order, risky_comb->combination, num_samples, num_samples_in_fold, folds[i], stride, block_starts);
            
                // Function for getting the matrix containing {FP,FN,TP,TN}
                unsigned int *confusion_matrix = get_confusion_matrix(order, risky_comb, aff_per_fold, unaff_per_fold, val);
                
                printf("confusion matrix = { ");
                for (int k = 0; k < 4; k++) {
                    printf("%u ", confusion_matrix[k]);
                }
                printf("}\n");
                
                // Function for evaluating the model, based on the confusion matrix
                double eval = evaluate_model(confusion_matrix, BA);
                
                printf("risky combination = {\n  SNP: ");
                print_combination(risky_comb->combination, 0, order);
                printf("  GT: ");
                for (int j = 0; j < risky_comb->num_risky * 2; j++) {
                    if (j % 2) {
                        printf("%d), ", risky_comb->genotypes[j]);
                    } else {
                        printf("(%d ", risky_comb->genotypes[j]);
                    }
                }
                printf("\n  Balanced accuracy: %.3f\n}\n", eval);

            }
            printf("-------------\n");
        }
        
        
        
        
        printf("\n-----------\n");
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
    TCase *tc_k_fold = tcase_create("k-fold creation");
    tcase_add_test(tc_k_fold, test_get_k_folds);
    
    TCase *tc_genotypes = tcase_create("Genotype and fold association");
    tcase_add_test(tc_genotypes, test_get_genotypes_for_combination_and_fold);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("k-fold cross validation");
    suite_add_tcase(fs, tc_k_fold);
    suite_add_tcase(fs, tc_genotypes);
    
    return fs;
}


/* ******************************
 *          Auxiliary           *
 * *****************************/

void check_copied_genotypes_fold_0(uint8_t *genotypes, uint8_t *block_starts[2], int snp_offset[2], int num_samples) {
    // Test first SNP
//     printf("%u = %u\n", genotypes[0], *(block_starts[0] + snp_offset[0] * num_samples + 2));
    fail_if(genotypes[0] != *(block_starts[0] + snp_offset[0] * num_samples + 2), "F0: Value 0 must be aligned with offset 2");
    fail_if(genotypes[1] != *(block_starts[0] + snp_offset[0] * num_samples + 3), "F0: Value 1 must be aligned with offset 3");
    fail_if(genotypes[2] != *(block_starts[0] + snp_offset[0] * num_samples + 4), "F0: Value 2 must be aligned with offset 4");
    fail_if(genotypes[3] != *(block_starts[0] + snp_offset[0] * num_samples + 5), "F0: Value 3 must be aligned with offset 5");
    fail_if(genotypes[4] != *(block_starts[0] + snp_offset[0] * num_samples + 8), "F0: Value 4 must be aligned with offset 8");
    fail_if(genotypes[5] != *(block_starts[0] + snp_offset[0] * num_samples + 9), "F0: Value 5 must be aligned with offset 9");
    fail_if(genotypes[6] != *(block_starts[0] + snp_offset[0] * num_samples + 10), "F0: Value 6 must be aligned with offset 10");
    fail_if(genotypes[7] != *(block_starts[0] + snp_offset[0] * num_samples + 11), "F0: Value 7 must be aligned with offset 11");

    // Test second SNP
//     printf("%u = %d\n", genotypes[8], *(block_starts[1] + snp_offset[1] * num_samples + 2));
    fail_if(genotypes[8] != *(block_starts[1] + snp_offset[1] * num_samples + 2), "F0: Value 8 must be aligned with offset 2");
    fail_if(genotypes[9] != *(block_starts[1] + snp_offset[1] * num_samples + 3), "F0: Value 9 must be aligned with offset 3");
    fail_if(genotypes[10] != *(block_starts[1] + snp_offset[1] * num_samples + 4), "F0: Value 10 must be aligned with offset 4");
    fail_if(genotypes[11] != *(block_starts[1] + snp_offset[1] * num_samples + 5), "F0: Value 11 must be aligned with offset 5");
    fail_if(genotypes[12] != *(block_starts[1] + snp_offset[1] * num_samples + 8), "F0: Value 12 must be aligned with offset 8");
    fail_if(genotypes[13] != *(block_starts[1] + snp_offset[1] * num_samples + 9), "F0: Value 13 must be aligned with offset 9");
    fail_if(genotypes[14] != *(block_starts[1] + snp_offset[1] * num_samples + 10), "F0: Value 14 must be aligned with offset 10");
    fail_if(genotypes[15] != *(block_starts[1] + snp_offset[1] * num_samples + 11), "F0: Value 15 must be aligned with offset 11");
}

void check_copied_genotypes_fold_1(uint8_t *genotypes, uint8_t *block_starts[2], int snp_offset[2], int num_samples) {
    // Test first SNP
//     printf("%u = %u\n", genotypes[0], *(block_starts[0] + snp_offset[0] * num_samples + 2));
    fail_if(genotypes[0] != *(block_starts[0] + snp_offset[0] * num_samples + 0), "F1: Value 0 must be aligned with offset 0");
    fail_if(genotypes[1] != *(block_starts[0] + snp_offset[0] * num_samples + 1), "F1: Value 1 must be aligned with offset 1");
    fail_if(genotypes[2] != *(block_starts[0] + snp_offset[0] * num_samples + 2), "F1: Value 2 must be aligned with offset 2");
    fail_if(genotypes[3] != *(block_starts[0] + snp_offset[0] * num_samples + 3), "F1: Value 3 must be aligned with offset 3");
    fail_if(genotypes[4] != *(block_starts[0] + snp_offset[0] * num_samples + 6), "F1: Value 4 must be aligned with offset 6");
    fail_if(genotypes[5] != *(block_starts[0] + snp_offset[0] * num_samples + 7), "F1: Value 5 must be aligned with offset 7");
    fail_if(genotypes[6] != *(block_starts[0] + snp_offset[0] * num_samples + 8), "F1: Value 6 must be aligned with offset 8");
    fail_if(genotypes[7] != *(block_starts[0] + snp_offset[0] * num_samples + 9), "F1: Value 7 must be aligned with offset 9");
    
    // Test second SNP
//     printf("%u = %d\n", genotypes[8], *(block_starts[1] + snp_offset[1] * num_samples + 2));
    fail_if(genotypes[8] != *(block_starts[1] + snp_offset[1] * num_samples + 0), "F1: Value 8 must be aligned with offset 0");
    fail_if(genotypes[9] != *(block_starts[1] + snp_offset[1] * num_samples + 1), "F1: Value 9 must be aligned with offset 1");
    fail_if(genotypes[10] != *(block_starts[1] + snp_offset[1] * num_samples + 2), "F1: Value 10 must be aligned with offset 2");
    fail_if(genotypes[11] != *(block_starts[1] + snp_offset[1] * num_samples + 3), "F1: Value 11 must be aligned with offset 3");
    fail_if(genotypes[12] != *(block_starts[1] + snp_offset[1] * num_samples + 6), "F1: Value 12 must be aligned with offset 6");
    fail_if(genotypes[13] != *(block_starts[1] + snp_offset[1] * num_samples + 7), "F1: Value 13 must be aligned with offset 7");
    fail_if(genotypes[14] != *(block_starts[1] + snp_offset[1] * num_samples + 8), "F1: Value 14 must be aligned with offset 8");
    fail_if(genotypes[15] != *(block_starts[1] + snp_offset[1] * num_samples + 9), "F1: Value 15 must be aligned with offset 9");
}

void check_copied_genotypes_fold_2(uint8_t *genotypes, uint8_t *block_starts[2], int snp_offset[2], int num_samples) {
    // Test first SNP
//     printf("%u = %u\n", genotypes[0], *(block_starts[0] + snp_offset[0] * num_samples + 2));
    fail_if(genotypes[0] != *(block_starts[0] + snp_offset[0] * num_samples + 0), "F2: Value 0 must be aligned with offset 0");
    fail_if(genotypes[1] != *(block_starts[0] + snp_offset[0] * num_samples + 1), "F2: Value 1 must be aligned with offset 1");
    fail_if(genotypes[2] != *(block_starts[0] + snp_offset[0] * num_samples + 4), "F2: Value 2 must be aligned with offset 4");
    fail_if(genotypes[3] != *(block_starts[0] + snp_offset[0] * num_samples + 5), "F2: Value 3 must be aligned with offset 5");
    fail_if(genotypes[4] != *(block_starts[0] + snp_offset[0] * num_samples + 6), "F2: Value 4 must be aligned with offset 6");
    fail_if(genotypes[5] != *(block_starts[0] + snp_offset[0] * num_samples + 7), "F2: Value 5 must be aligned with offset 7");
    fail_if(genotypes[6] != *(block_starts[0] + snp_offset[0] * num_samples + 10), "F2: Value 6 must be aligned with offset 10");
    fail_if(genotypes[7] != *(block_starts[0] + snp_offset[0] * num_samples + 11), "F2: Value 7 must be aligned with offset 11");
    
    // Test second SNP
//     printf("%u = %d\n", genotypes[8], *(block_starts[1] + snp_offset[1] * num_samples + 2));
    fail_if(genotypes[8] != *(block_starts[1] + snp_offset[1] * num_samples + 0), "F2: Value 8 must be aligned with offset 0");
    fail_if(genotypes[9] != *(block_starts[1] + snp_offset[1] * num_samples + 1), "F2: Value 9 must be aligned with offset 1");
    fail_if(genotypes[10] != *(block_starts[1] + snp_offset[1] * num_samples + 4), "F2: Value 10 must be aligned with offset 4");
    fail_if(genotypes[11] != *(block_starts[1] + snp_offset[1] * num_samples + 5), "F2: Value 11 must be aligned with offset 5");
    fail_if(genotypes[12] != *(block_starts[1] + snp_offset[1] * num_samples + 6), "F2: Value 12 must be aligned with offset 6");
    fail_if(genotypes[13] != *(block_starts[1] + snp_offset[1] * num_samples + 7), "F2: Value 13 must be aligned with offset 7");
    fail_if(genotypes[14] != *(block_starts[1] + snp_offset[1] * num_samples + 10), "F2: Value 14 must be aligned with offset 10");
    fail_if(genotypes[15] != *(block_starts[1] + snp_offset[1] * num_samples + 11), "F2: Value 15 must be aligned with offset 11");
}
