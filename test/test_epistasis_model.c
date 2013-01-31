#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <check.h>

#include <containers/array_list.h>

#include "epistasis/model.h"


Suite *create_test_suite(void);

int num_samples = 8;
int num_affected = 4;
int num_unaffected = 4;

uint8_t genotypes[] = { 0, 0, 1, 0, 2, 1, 0, 2, 
                        0, 1, 1, 0, 0, 0, 1, 1,
                        1, 2, 0, 1, 0, 2, 0, 0,
                        0, 0, 0, 2, 1, 1, 0, 2 };


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/



/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST(test_get_masks) {
    int order = 4;
    
    int num_masks;
    uint8_t *masks = get_masks(order, genotypes, num_samples, &num_masks);
    
    fail_unless(num_masks == 3 * order * num_samples, "There should be 3 GTs * 4 SNPs * 8 samples = 96 elements in the masks array");
    
    uint8_t mask_0[] = { 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1 };
    uint8_t mask_1[] = { 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t mask_2[] = { 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0 };
    uint8_t mask_3[] = { 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1 };
    
//     for (int i = 0; i < order; i++) {
//         print_gt_combination(masks + i * 3 * num_samples, i, 3 * num_samples);
//     }
//     printf("\n");
//     print_gt_combination(mask_0, 0, 3 * num_samples);
//     print_gt_combination(mask_1, 1, 3 * num_samples);
//     print_gt_combination(mask_2, 2, 3 * num_samples);
//     print_gt_combination(mask_3, 3, 3 * num_samples);
    
    for (int j = 0; j < num_samples; j++) {
        char msg[26];
        sprintf(msg, "Mask SNP0: %d should be %d\n", j, mask_0[j]);
        fail_if(masks[j] != mask_0[j], msg);
    }
    
    for (int j = 0; j < num_samples; j++) {
        char msg[26];
        sprintf(msg, "Mask SNP1: %d should be %d\n", j, mask_1[j]);
        fail_if(masks[num_samples * 3 + j] != mask_1[j], msg);
    }
    
    for (int j = 0; j < num_samples; j++) {
        char msg[26];
        sprintf(msg, "Mask SNP2: %d should be %d\n", j, mask_2[j]);
        fail_if(masks[num_samples * 6 + j] != mask_2[j], msg);
    }
    
    for (int j = 0; j < num_samples; j++) {
        char msg[26];
        sprintf(msg, "Mask SNP3: %d should be %d\n", j, mask_3[j]);
        fail_if(masks[num_samples * 9 + j] != mask_3[j], msg);
    }
}
END_TEST

START_TEST(test_get_counts) {
    int order = 2;
    
    int num_counts;
    int num_combinations;
    uint8_t **combinations = get_genotype_combinations(order, &num_combinations);
    
    // Order 2
    int* counts = get_counts(order, genotypes, combinations, num_combinations, num_affected, num_unaffected, &num_counts);
    fail_if(num_counts != 18, "There must be 9 order 2 combinations * 2 phenotypes");
    
    fail_if(counts[0] != 2 && counts[1] != 0, "A(0,0) -> 2,0");
    fail_if(counts[2] != 1 && counts[3] != 1, "A(0,1) -> 1,1");
    fail_if(counts[4] != 0 && counts[5] != 0, "A(0,2) -> 0,0");
    fail_if(counts[6] != 0 && counts[7] != 1, "A(1,0) -> 0,1");
    fail_if(counts[8] != 1 && counts[9] != 0, "A(1,1) -> 1,0");
    fail_if(counts[10] != 0 && counts[11] != 0, "A(1,2) -> 0,0");
    fail_if(counts[12] != 0 && counts[13] != 1, "A(2,0) -> 0,1");
    fail_if(counts[14] != 0 && counts[15] != 1, "A(2,1) -> 0,1");
    fail_if(counts[16] != 0 && counts[17] != 0, "A(2,2) -> 0,0");
    
    free(combinations);
    
    // Order 3
    order = 3;
    combinations = get_genotype_combinations(order, &num_combinations);
    counts = get_counts(order, genotypes, combinations, num_combinations, num_affected, num_unaffected, &num_counts);
    fail_if(num_counts != 54, "There must be 27 order 3 combinations * 2 phenotypes");
    
    int idx = 0;
    fail_if(counts[idx] != 0 && counts[idx + 1] != 0, "A(0,0,0) -> 0,0");
    fail_if(counts[idx + 2] != 2 && counts[idx + 3] != 0, "A(0,0,1) -> 2,0");
    fail_if(counts[idx + 4] != 0 && counts[idx + 5] != 0, "A(0,0,2) -> 0,0");
    
    fail_if(counts[idx + 6] != 0 && counts[idx + 7] != 1, "A(0,1,0) -> 0,1");
    fail_if(counts[idx + 8] != 0 && counts[idx + 9] != 0, "A(0,1,1) -> 0,0");
    fail_if(counts[idx + 10] != 1 && counts[idx + 11] != 0, "A(0,1,2) -> 1,0");
    
    fail_if(counts[idx + 12] != 0 && counts[idx + 13] != 0, "A(0,2,0) -> 0,0");
    fail_if(counts[idx + 14] != 0 && counts[idx + 15] != 0, "A(0,2,1) -> 0,0");
    fail_if(counts[idx + 16] != 0 && counts[idx + 17] != 0, "A(0,2,2) -> 0,0");
    
    idx = 18;
    fail_if(counts[idx] != 0 && counts[idx + 1] != 0, "A(1,0,0) -> 0,0");
    fail_if(counts[idx + 2] != 0 && counts[idx + 3] != 0, "A(1,0,1) -> 0,0");
    fail_if(counts[idx + 4] != 0 && counts[idx + 5] != 1, "A(1,0,2) -> 0,1");
    
    fail_if(counts[idx + 6] != 1 && counts[idx + 7] != 0, "A(1,1,0) -> 1,0");
    fail_if(counts[idx + 8] != 0 && counts[idx + 9] != 0, "A(1,1,1) -> 0,0");
    fail_if(counts[idx + 10] != 0 && counts[idx + 11] != 0, "A(1,1,2) -> 0,0");
    
    fail_if(counts[idx + 12] != 0 && counts[idx + 13] != 0, "A(1,2,0) -> 0,0");
    fail_if(counts[idx + 14] != 0 && counts[idx + 15] != 0, "A(1,2,1) -> 0,0");
    fail_if(counts[idx + 16] != 0 && counts[idx + 17] != 0, "A(1,2,2) -> 0,0");
    
    idx = 36;
    fail_if(counts[idx] != 0 && counts[idx + 1] != 1, "A(2,0,0) -> 0,1");
    fail_if(counts[idx + 2] != 0 && counts[idx + 3] != 0, "A(2,0,1) -> 0,0");
    fail_if(counts[idx + 4] != 0 && counts[idx + 5] != 0, "A(2,0,2) -> 0,0");
    
    fail_if(counts[idx + 6] != 0 && counts[idx + 7] != 1, "A(2,1,0) -> 0,1");
    fail_if(counts[idx + 8] != 0 && counts[idx + 9] != 0, "A(2,1,1) -> 0,0");
    fail_if(counts[idx + 10] != 0 && counts[idx + 11] != 0, "A(2,1,2) -> 0,0");
    
    fail_if(counts[idx + 12] != 0 && counts[idx + 13] != 0, "A(2,2,0) -> 0,0");
    fail_if(counts[idx + 14] != 0 && counts[idx + 15] != 0, "A(2,2,1) -> 0,0");
    fail_if(counts[idx + 16] != 0 && counts[idx + 17] != 0, "A(2,2,2) -> 0,0");
}
END_TEST


START_TEST(test_get_confusion_matrix) {
    int order = 2;
    int num_combinations;
    
    // ---------- order 2 ------------
    
    uint8_t genotypes_2d[] = { 1, 1, 0, 2, 2, 2, 1, 0, 0, 0, 1, 2,
                               0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2 };
    uint8_t **possible_2d = get_genotype_combinations(order, &num_combinations);
    // Risky combinations: (0,1), (2,1), (2,2)
    risky_combination *combination_2d = risky_combination_new(order, (int[2]){ 0, 1 }, possible_2d, 3, (int[3]){ 3, 7, 8 } );
    
    // 7 affected, 5 unaffected
    int *confusion_matrix = get_confusion_matrix(order, combination_2d, 7, 5, genotypes_2d);
    
    fail_if(confusion_matrix[0] != 6, "(7 aff,5 unaff) TP = 6");
    fail_if(confusion_matrix[1] != 1, "(7 aff,5 unaff) FN = 1");
    fail_if(confusion_matrix[2] != 1, "(7 aff,5 unaff) FP = 1");
    fail_if(confusion_matrix[3] != 4, "(7 aff,5 unaff) TN = 4");
    
    // 4 affected, 8 unaffected
    confusion_matrix = get_confusion_matrix(order, combination_2d, 4, 8, genotypes_2d);
    
    fail_if(confusion_matrix[0] != 3, "(4 aff,8 unaff) TP = 3");
    fail_if(confusion_matrix[1] != 1, "(4 aff,8 unaff) FN = 1");
    fail_if(confusion_matrix[2] != 4, "(4 aff,8 unaff) FP = 4");
    fail_if(confusion_matrix[3] != 4, "(4 aff,8 unaff) TN = 4");
    
    free(confusion_matrix);
    risky_combination_free(combination_2d);
    
    // ---------- order 3 ------------
    order = 3;
    uint8_t genotypes_3d[] = { 1, 1, 0, 2, 2, 2, 1, 0, 0, 0, 1, 2,
                               0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2,
                               1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0 };
    uint8_t **possible_3d = get_genotype_combinations(order, &num_combinations);
    // Risky combinations: (1,0,1), (2,1,0), (2,2,1)
    risky_combination *combination_3d = risky_combination_new(order, (int[3]){ 0, 1, 2 }, possible_3d, 4, (int[4]){ 4, 10, 21, 25 } );
    
    // 6 affected, 6 unaffected
    confusion_matrix = get_confusion_matrix(order, combination_3d, 6, 6, genotypes_3d);
    
    fail_if(confusion_matrix[0] != 6, "(6 aff,6 unaff) TP = 6");
    fail_if(confusion_matrix[1] != 0, "(6 aff,6 unaff) FN = 0");
    fail_if(confusion_matrix[2] != 3, "(6 aff,6 unaff) FP = 3");
    fail_if(confusion_matrix[3] != 3, "(6 aff,6 unaff) TN = 3");
    
    free(confusion_matrix);
    risky_combination_free(combination_3d);
                               
}
END_TEST


START_TEST(test_model_evaluation_formulas) {
    const unsigned int TP_1 = 40, FP_1 = 4, FN_1 = 2, TN_1 = 10;
    const unsigned int TP_2 = 20, FP_2 = 10, FN_2 = 10, TN_2 = 20;
    
    // Accuracy
    fail_if(evaluate_model(TP_1, TN_1, FP_1, FN_1, CA) - 0.89285714 > 1e-6, "CA(1) = 0.89285714");
    fail_if(evaluate_model(TP_2, TN_2, FP_2, FN_2, CA) - 0.66666666 > 1e-6, "CA(2) = 0.6666...");
    
    // Balanced accuracy
    fail_if(evaluate_model(TP_1, TN_1, FP_1, FN_1, BA) - 0.83333333 > 1e-6, "BA(1) = 0.83333...");
    fail_if(evaluate_model(TP_2, TN_2, FP_2, FN_2, BA) - 0.66666666 > 1e-6, "BA(2) = 0.6666...");
    
    // Gamma
    fail_if(evaluate_model(TP_1, TN_1, FP_1, FN_1, GAMMA) - 0.96078431 > 1e-6, "GAMMA(1) = 0.96078431");
    fail_if(evaluate_model(TP_2, TN_2, FP_2, FN_2, GAMMA) - 0.6 > 1e-6, "GAMMA(2) = 0.6");
    
    // Kendall's Tau-B
    fail_if(evaluate_model(TP_1, TN_1, FP_1, FN_1, TAU_B) - 0.70352647 > 1e-6, "TAU_B(1) = 0.70352647");
    fail_if(evaluate_model(TP_2, TN_2, FP_2, FN_2, TAU_B) - 0.33333333 > 1e-6, "TAU_B(2) = 0.3333...");
    
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
    TCase *tc_counts = tcase_create("Genotype counts");
    tcase_add_test(tc_counts, test_get_masks);
    tcase_add_test(tc_counts, test_get_counts);
    
    TCase *tc_ranking = tcase_create("Evaluation and ranking");
    tcase_add_test(tc_ranking, test_get_confusion_matrix);
    tcase_add_test(tc_ranking, test_model_evaluation_formulas);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Epistasis model");
    suite_add_tcase(fs, tc_counts);
    suite_add_tcase(fs, tc_ranking);
    
    return fs;
}
