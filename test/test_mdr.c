#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <check.h>

#include <containers/array_list.h>

#include "epistasis/mdr.h"
#include "epistasis/model.h"


Suite *create_test_suite(void);


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
    int folds[3][4] = { { 0, 1, 6, 7 }, { 4, 5, 10, 11 }, { 8, 9, 2, 3 } };
    
    uint8_t genotypes[] = { 0, 1, 1, 0, 1, 2, 2, 1, 2, 1, 1, 2,   // First block
                            1, 2, 2, 1, 2, 0, 0, 2, 0, 2, 2, 0,
                            2, 0, 0, 2, 0, 1, 1, 0, 1, 0, 0, 1,
                            0, 1, 0, 1, 1, 2, 1, 2, 2, 0, 2, 0,   // Second block
                            1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1,
                            1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 1, 2,
                            1, 0, 1, 2, 0, 2, 0, 2, 0, 2, 2, 1,   // Third block
                            1, 2, 2, 1, 1, 0, 1, 0, 1, 2, 0, 2,
                            2, 0, 0, 0, 1, 2, 0, 2, 2, 1, 1, 0 };
        
    printf("Genotypes = {\n");
    for (int i = 0; i < 9; i++) {
        printf("\t{ ");
        for (int j = 0; j < 12; j++) {
            printf("%d ", genotypes[i * 12 + j]);
        }
        printf(" }\n");
    }
    printf(" }\n\n");
    
    int block_2d[] = { 0, 0 };
    do {
        printf("BLOCK %d,%d\n", block_2d[0], block_2d[1]);
        uint8_t *block_starts[2];
        block_starts[0] = genotypes + block_2d[0] * stride * num_samples;
        block_starts[1] = genotypes + block_2d[1] * stride * num_samples;
        
        // Test first combination in the block
        int *comb = get_first_combination_in_block(order, block_2d, stride);
        int snp_offset[] = { comb[0] % stride, comb[1] % stride };
        
        // Run for each fold
        for (int i = 0; i < num_folds; i++) {
            int next_fold = (i < 2) ? i + 1 : 0 ;
            int *val = get_genotypes_for_combination_and_fold(order, comb, num_samples, num_samples_in_fold, folds[i], block_2d, stride, block_starts);
           
//             printf("Combination (%d,%d) and fold %d\n", comb[0], comb[1], i);
            
            // Test first SNP
//             printf("%d = %d\n", val[0], *(block_starts[0] + snp_offset[0] * num_samples + i * num_samples_in_fold));
            fail_if(val[0] != *(block_starts[0] + snp_offset[0] * num_samples + i * num_samples_in_fold), "Value 0 must be aligned with offset: fold * 4");
            fail_if(val[1] != *(block_starts[0] + snp_offset[0] * num_samples + i * num_samples_in_fold + 1), "Value 1 must be aligned with offset: fold * 4 + 1");
            fail_if(val[2] != *(block_starts[0] + snp_offset[0] * num_samples + next_fold * num_samples_in_fold + 2), "Value 2 must be aligned with offset: (fold+1) * 4 + 2");
            fail_if(val[3] != *(block_starts[0] + snp_offset[0] * num_samples + next_fold * num_samples_in_fold + 3), "Value 3 must be aligned with offset: (fold+1) * 4 + 3");
            
            // Test second SNP
//             printf("%d = %d\n", val[4], *(block_starts[1] + snp_offset[1] * num_samples + i * num_samples_in_fold));
            fail_if(val[4] != *(block_starts[1] + snp_offset[1] * num_samples + i * num_samples_in_fold), "Value 4 must be aligned with offset: fold * 4");
            fail_if(val[5] != *(block_starts[1] + snp_offset[1] * num_samples + i * num_samples_in_fold + 1), "Value 5 must be aligned with offset: fold * 4 + 1");
            fail_if(val[6] != *(block_starts[1] + snp_offset[1] * num_samples + next_fold * num_samples_in_fold + 2), "Value 6 must be aligned with offset: (fold+1) * 4 + 2");
            fail_if(val[7] != *(block_starts[1] + snp_offset[1] * num_samples + next_fold * num_samples_in_fold + 3), "Value 7 must be aligned with offset: (fold+1) * 4 + 3");
            
            free(val);
        }
        
        // Test next combinations
        while (get_next_combination_in_block(order, comb, block_2d, stride)) {
            int snp_offset[] = { comb[0] % stride, comb[1] % stride };
            
            for (int i = 0; i < num_folds; i++) {
                int next_fold = (i < 2) ? i + 1 : 0 ;
                int *val = get_genotypes_for_combination_and_fold(order, comb, num_samples, num_samples_in_fold, folds[i], block_2d, stride, block_starts);
            
                printf("Combination (%d,%d) and fold %d\n", comb[0], comb[1], i);
                
                // Test first SNP
//                 printf("%d = %d\n", val[0], *(block_starts[0] + snp_offset[0] * num_samples + i * num_samples_in_fold));
                fail_if(val[0] != *(block_starts[0] + snp_offset[0] * num_samples + i * num_samples_in_fold), "Value 0 must be aligned with offset: fold * 4");
                fail_if(val[1] != *(block_starts[0] + snp_offset[0] * num_samples + i * num_samples_in_fold + 1), "Value 1 must be aligned with offset: fold * 4 + 1");
                fail_if(val[2] != *(block_starts[0] + snp_offset[0] * num_samples + next_fold * num_samples_in_fold + 2), "Value 2 must be aligned with offset: (fold+1) * 4 + 2");
                fail_if(val[3] != *(block_starts[0] + snp_offset[0] * num_samples + next_fold * num_samples_in_fold + 3), "Value 3 must be aligned with offset: (fold+1) * 4 + 3");
                
                // Test second SNP
//                 printf("[%d] %d = %d\n", val[4], *(block_starts[1] + snp_offset[1] * num_samples + i * num_samples_in_fold));
                fail_if(val[4] != *(block_starts[1] + snp_offset[1] * num_samples + i * num_samples_in_fold), "Value 4 must be aligned with offset: fold * 4");
                fail_if(val[5] != *(block_starts[1] + snp_offset[1] * num_samples + i * num_samples_in_fold + 1), "Value 5 must be aligned with offset: fold * 4 + 1");
                fail_if(val[6] != *(block_starts[1] + snp_offset[1] * num_samples + next_fold * num_samples_in_fold + 2), "Value 6 must be aligned with offset: (fold+1) * 4 + 2");
                fail_if(val[7] != *(block_starts[1] + snp_offset[1] * num_samples + next_fold * num_samples_in_fold + 3), "Value 7 must be aligned with offset: (fold+1) * 4 + 3");
                
                free(val);
            }
        }
        
        free(comb);
        
        printf("\n-----------\n");
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
    // 100 { 30 aff, 70 unaff} in 10-fold
    const int num_variants = 100, block_size = 10, num_blocks = 10;
    const int num_affected = 30, num_unaffected = 70;
    
    const int order = 2;
    const int num_folds = 10;
    unsigned int *sizes;
    
    int num_counts;
    int* counts;
    
    int **folds = get_k_folds(num_affected, num_unaffected, num_folds, &sizes);
    
    
    int comb_idx = 0;
    int block_2d[] = { 0, 0 };
    int block_3d[] = { 0, 0, 0 };
    do {
        printf("BLOCK %d,%d\n", block_2d[0], block_2d[1]);
//         printf("BLOCK %d,%d,%d\n", block_3d[0], block_3d[1], block_3d[2]);
        
        // TODO test first combination in the block
        int *comb = get_first_combination_in_block(order, block_2d, 10);
        print_combination(comb, comb_idx, order);
        
        // TODO run for each fold
        for (int i = 0; i < num_folds; i++) {
            // TODO get genotypes of that combination
            /*
             * arguments needed:
             * 1) order
             * 2) combination
             * 3) block coordinates
             * 4) stride
             * 5) block start
             * 6) fold samples
             * 7) NUMBER of fold samples
             */
           
            
            // TODO get counts for that combination
//             counts = get_counts(order, genotypes, num_affected, num_unaffected, &num_counts);
            
            // TODO get high risk pairs for these counts
            
            
//             free(counts);
        }
        
        // TODO get next combination
        while (get_next_combination_in_block(order, comb, block_2d, 10)) {
            print_combination(comb, comb_idx * 10, order);
        
            // TODO get counts for that combination
//             counts = get_counts(order, genotypes, num_affected, num_unaffected, &num_counts);
            
            // TODO get high risk pairs for these counts
            
//             free(counts);
        }
        
        free(comb);
        
        comb_idx++;
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
    TCase *tc_k_fold = tcase_create("k-fold and cross-validation");
    tcase_add_test(tc_k_fold, test_get_k_folds);
    tcase_add_test(tc_k_fold, test_get_genotypes_for_combination_and_fold);
//     tcase_add_test(tc_k_fold, test_mdr_steps_2_5);
    
    TCase *tc_counts = tcase_create("Risk classification");
    tcase_add_test(tc_counts, test_get_high_risk_combinations);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("MDR");
    suite_add_tcase(fs, tc_k_fold);
    suite_add_tcase(fs, tc_counts);
    
    return fs;
}
