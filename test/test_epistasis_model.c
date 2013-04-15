#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <check.h>

#include <containers/array_list.h>

#include "gwas/epistasis/model.h"


Suite *create_test_suite(void);

int num_samples = 8;
int num_affected = 4;
int num_unaffected = 4;

uint8_t genotypes[4][8] = { { 0, 0, 1, 0, 2, 1, 0, 2 }, 
                            { 0, 1, 1, 0, 0, 0, 1, 1 },
                            { 1, 2, 0, 1, 0, 2, 0, 0 },
                            { 0, 0, 0, 2, 1, 1, 0, 2 } };

/* ******************************
 *      Unchecked fixtures      *
 * ******************************/



/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST(test_get_masks) {
    int order = 4;
    int num_samples_with_padding = 32;
    
    // Genotypes of all SNPs
    uint8_t masks_gt0[] = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt1[] = { 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt2[] = { 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt3[] = { 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *masks_genotypes[] = { masks_gt0, masks_gt1, masks_gt2, masks_gt3 };
    
    masks_info info; masks_info_new(order, num_affected, num_unaffected, &info);
    
    uint8_t *masks = set_genotypes_masks(order, masks_genotypes, info);
    
    fail_unless(info.num_masks == 3 * order * num_samples_with_padding, "There should be 3 GTs * 4 SNPs * 4+4 samples (padded) = 384 elements in the masks array");
    
    uint8_t mask_0[] = { 255, 255, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    // Genotype 0
                         0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    // Genotype 1
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };  // Genotype 2
    uint8_t mask_1[] = { 255, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t mask_2[] = { 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         255, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t mask_3[] = { 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
/*
     for (int i = 0; i < order; i++) {
         print_gt_combination(masks + i * 3 * num_samples_with_padding, i, 3 * num_samples_with_padding);
     }
     printf("\n");
     print_gt_combination(mask_0, 0, 3 * num_samples_with_padding);
     print_gt_combination(mask_1, 1, 3 * num_samples_with_padding);
     print_gt_combination(mask_2, 2, 3 * num_samples_with_padding);
     print_gt_combination(mask_3, 3, 3 * num_samples_with_padding);
*/
    
    for (int j = 0; j < num_samples; j++) {
        char msg[26];
        sprintf(msg, "Mask SNP0: %d should be %d\n", j, mask_0[j]);
        fail_if(masks[j] != mask_0[j], msg);
    }
    
    for (int j = 0; j < num_samples; j++) {
        char msg[26];
        sprintf(msg, "Mask SNP1: %d should be %d\n", j, mask_1[j]);
        fail_if(masks[num_samples_with_padding * 3 + j] != mask_1[j], msg);
    }
    
    for (int j = 0; j < num_samples; j++) {
        char msg[26];
        sprintf(msg, "Mask SNP2: %d should be %d\n", j, mask_2[j]);
        fail_if(masks[num_samples_with_padding * 6 + j] != mask_2[j], msg);
    }
    
    for (int j = 0; j < num_samples; j++) {
        char msg[26];
        sprintf(msg, "Mask SNP3: %d should be %d\n", j, mask_3[j]);
        fail_if(masks[num_samples_with_padding * 9 + j] != mask_3[j], msg);
    }
}
END_TEST

START_TEST(test_get_counts) {
    int order = 2;
    
    // Genotypes of all SNPs
    uint8_t masks_gt0[] = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt1[] = { 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt2[] = { 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *masks_genotypes[] = { masks_gt0, masks_gt1, masks_gt2 };
    
    masks_info info; masks_info_new(order, num_affected, num_unaffected, &info);
    
    uint8_t *masks = set_genotypes_masks(order, masks_genotypes, info);
    
    int num_combinations;
    uint8_t **combinations = get_genotype_combinations(order, &num_combinations);
    
    // Order 2
    int num_counts = 2 * pow(NUM_GENOTYPES, order);
    int counts_o2[num_counts];
    combination_counts(order, masks, combinations, num_combinations, counts_o2, info);
    fail_if(num_counts != 18, "There must be 9 order 2 combinations * 2 phenotypes");
    
    fail_if(counts_o2[0] != 2 || counts_o2[1] != 0, "A(0,0) -> 2,0");
    fail_if(counts_o2[2] != 1 || counts_o2[3] != 1, "A(0,1) -> 1,1");
    fail_if(counts_o2[4] != 0 || counts_o2[5] != 0, "A(0,2) -> 0,0");
    fail_if(counts_o2[6] != 0 || counts_o2[7] != 1, "A(1,0) -> 0,1");
    fail_if(counts_o2[8] != 1 || counts_o2[9] != 0, "A(1,1) -> 1,0");
    fail_if(counts_o2[10] != 0 || counts_o2[11] != 0, "A(1,2) -> 0,0");
    fail_if(counts_o2[12] != 0 || counts_o2[13] != 1, "A(2,0) -> 0,1");
    fail_if(counts_o2[14] != 0 || counts_o2[15] != 1, "A(2,1) -> 0,1");
    fail_if(counts_o2[16] != 0 || counts_o2[17] != 0, "A(2,2) -> 0,0");
    
    free(combinations);
    
    // Order 3
    order = 3;
    combinations = get_genotype_combinations(order, &num_combinations);
    masks = set_genotypes_masks(order, masks_genotypes, info);
    num_counts = 2 * pow(NUM_GENOTYPES, order);
    int counts_o3[num_counts];
    combination_counts(order, masks, combinations, num_combinations, counts_o3, info);
    fail_if(num_counts != 54, "There must be 27 order 3 combinations * 2 phenotypes");
    
//     for (int i = 0; i < order; i++) {
//         print_gt_combination(masks + i * 3 * info.num_samples_per_mask, i, 3 * info.num_samples_per_mask);
//     }
    
    int idx = 0;
    fail_if(counts_o3[idx] != 0 || counts_o3[idx + 1] != 0, "A(0,0,0) -> 0,0");
    fail_if(counts_o3[idx + 2] != 2 || counts_o3[idx + 3] != 0, "A(0,0,1) -> 2,0");
    fail_if(counts_o3[idx + 4] != 0 || counts_o3[idx + 5] != 0, "A(0,0,2) -> 0,0");
    
    fail_if(counts_o3[idx + 6] != 0 || counts_o3[idx + 7] != 1, "A(0,1,0) -> 0,1");
    fail_if(counts_o3[idx + 8] != 0 || counts_o3[idx + 9] != 0, "A(0,1,1) -> 0,0");
    fail_if(counts_o3[idx + 10] != 1 || counts_o3[idx + 11] != 0, "A(0,1,2) -> 1,0");
    
    fail_if(counts_o3[idx + 12] != 0 || counts_o3[idx + 13] != 0, "A(0,2,0) -> 0,0");
    fail_if(counts_o3[idx + 14] != 0 || counts_o3[idx + 15] != 0, "A(0,2,1) -> 0,0");
    fail_if(counts_o3[idx + 16] != 0 || counts_o3[idx + 17] != 0, "A(0,2,2) -> 0,0");
    
    idx = 18;
    fail_if(counts_o3[idx] != 0 || counts_o3[idx + 1] != 0, "A(1,0,0) -> 0,0");
    fail_if(counts_o3[idx + 2] != 0 || counts_o3[idx + 3] != 0, "A(1,0,1) -> 0,0");
    fail_if(counts_o3[idx + 4] != 0 || counts_o3[idx + 5] != 1, "A(1,0,2) -> 0,1");
    
    fail_if(counts_o3[idx + 6] != 1 || counts_o3[idx + 7] != 0, "A(1,1,0) -> 1,0");
    fail_if(counts_o3[idx + 8] != 0 || counts_o3[idx + 9] != 0, "A(1,1,1) -> 0,0");
    fail_if(counts_o3[idx + 10] != 0 || counts_o3[idx + 11] != 0, "A(1,1,2) -> 0,0");
    
    fail_if(counts_o3[idx + 12] != 0 || counts_o3[idx + 13] != 0, "A(1,2,0) -> 0,0");
    fail_if(counts_o3[idx + 14] != 0 || counts_o3[idx + 15] != 0, "A(1,2,1) -> 0,0");
    fail_if(counts_o3[idx + 16] != 0 || counts_o3[idx + 17] != 0, "A(1,2,2) -> 0,0");
    
    idx = 36;
    fail_if(counts_o3[idx] != 0 || counts_o3[idx + 1] != 1, "A(2,0,0) -> 0,1");
    fail_if(counts_o3[idx + 2] != 0 || counts_o3[idx + 3] != 0, "A(2,0,1) -> 0,0");
    fail_if(counts_o3[idx + 4] != 0 || counts_o3[idx + 5] != 0, "A(2,0,2) -> 0,0");
    
    fail_if(counts_o3[idx + 6] != 0 || counts_o3[idx + 7] != 1, "A(2,1,0) -> 0,1");
    fail_if(counts_o3[idx + 8] != 0 || counts_o3[idx + 9] != 0, "A(2,1,1) -> 0,0");
    fail_if(counts_o3[idx + 10] != 0 || counts_o3[idx + 11] != 0, "A(2,1,2) -> 0,0");
    
    fail_if(counts_o3[idx + 12] != 0 || counts_o3[idx + 13] != 0, "A(2,2,0) -> 0,0");
    fail_if(counts_o3[idx + 14] != 0 || counts_o3[idx + 15] != 0, "A(2,2,1) -> 0,0");
    fail_if(counts_o3[idx + 16] != 0 || counts_o3[idx + 17] != 0, "A(2,2,2) -> 0,0");
    
    free(combinations);
}
END_TEST


START_TEST(test_get_confusion_matrix) {
    int order = 2;
    int num_combinations;
    masks_info info;
    
    // ---------- order 2 ------------
    
    uint8_t **possible_2d = get_genotype_combinations(order, &num_combinations);
    // Risky combinations: (1,0), (2,1), (2,2)
    risky_combination *combination_2d = risky_combination_new(order, (int[2]){ 0, 1 }, possible_2d, 3, (int[3]){ 3, 7, 8 }, NULL);
    
    // 7 affected, 5 unaffected
    uint8_t gt2_0a[] = { 1, 1, 0, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t gt2_1a[] = { 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *genotypes_2da[2] = { gt2_0a, gt2_1a };
    
    masks_info_new(order, 7, 5, &info);
    int *confusion_matrix = calculate_confusion_matrix(order, combination_2d, info, genotypes_2da);
    
    //printf("matrix = { %d, %d, %d, %d }\n", confusion_matrix[0], confusion_matrix[1], confusion_matrix[2], confusion_matrix[3]);
    fail_if(confusion_matrix[0] != 6, "(7 aff,5 unaff) TP = 6");
    fail_if(confusion_matrix[1] != 1, "(7 aff,5 unaff) FN = 1");
    fail_if(confusion_matrix[2] != 1, "(7 aff,5 unaff) FP = 1");
    fail_if(confusion_matrix[3] != 4, "(7 aff,5 unaff) TN = 4");
    
    // 4 affected, 8 unaffected
    uint8_t gt2_0b[] = { 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t gt2_1b[] = { 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *genotypes_2db[2] = { gt2_0b, gt2_1b };
    
    masks_info_new(order, 4, 8, &info);
    confusion_matrix = calculate_confusion_matrix(order, combination_2d, info, genotypes_2db);
    
    //printf("matrix = { %d, %d, %d, %d }\n", confusion_matrix[0], confusion_matrix[1], confusion_matrix[2], confusion_matrix[3]);
    fail_if(confusion_matrix[0] != 3, "(4 aff,8 unaff) TP = 3");
    fail_if(confusion_matrix[1] != 1, "(4 aff,8 unaff) FN = 1");
    fail_if(confusion_matrix[2] != 4, "(4 aff,8 unaff) FP = 4");
    fail_if(confusion_matrix[3] != 4, "(4 aff,8 unaff) TN = 4");
    
    free(confusion_matrix);
    risky_combination_free(combination_2d);
    
    // ---------- order 3 ------------
    order = 3;
    uint8_t gt3_0[] = { 1, 1, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t gt3_1[] = { 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t gt3_2[] = { 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *genotypes_3d[3] = { gt3_0, gt3_1, gt3_2 };
    
    uint8_t **possible_3d = get_genotype_combinations(order, &num_combinations);
    // Risky combinations: (1,0,1), (2,1,0), (2,2,1)
    risky_combination *combination_3d = risky_combination_new(order, (int[3]){ 0, 1, 2 }, possible_3d, 4, (int[4]){ 4, 10, 21, 25 }, NULL);
    
    // 6 affected, 6 unaffected
    masks_info_new(order, 6, 6, &info);
    confusion_matrix = calculate_confusion_matrix(order, combination_3d, info, genotypes_3d);
    
    //printf("matrix = { %d, %d, %d, %d }\n", confusion_matrix[0], confusion_matrix[1], confusion_matrix[2], confusion_matrix[3]);
    fail_if(confusion_matrix[0] != 6, "(6 aff,6 unaff) TP = 6");
    fail_if(confusion_matrix[1] != 0, "(6 aff,6 unaff) FN = 0");
    fail_if(confusion_matrix[2] != 3, "(6 aff,6 unaff) FP = 3");
    fail_if(confusion_matrix[3] != 3, "(6 aff,6 unaff) TN = 3");
    
    free(confusion_matrix);
    risky_combination_free(combination_3d);
                               
}
END_TEST


START_TEST(test_model_evaluation_formulas) {
    unsigned int confusion_matrix_1[] = { 40, 2, 4, 10 };
    unsigned int confusion_matrix_2[] = { 20, 10, 10, 20 };
    
    // Accuracy
    fail_if(evaluate_model(confusion_matrix_1, CA) - 0.89285714 > 1e-6, "CA(1) = 0.89285714");
    fail_if(evaluate_model(confusion_matrix_2, CA) - 0.66666666 > 1e-6, "CA(2) = 0.6666...");
    
    // Balanced accuracy
    fail_if(evaluate_model(confusion_matrix_1, BA) - 0.83333333 > 1e-6, "BA(1) = 0.83333...");
    fail_if(evaluate_model(confusion_matrix_2, BA) - 0.66666666 > 1e-6, "BA(2) = 0.6666...");
    
    // Gamma
    fail_if(evaluate_model(confusion_matrix_1, GAMMA) - 0.96078431 > 1e-6, "GAMMA(1) = 0.96078431");
    fail_if(evaluate_model(confusion_matrix_2, GAMMA) - 0.6 > 1e-6, "GAMMA(2) = 0.6");
    
    // Kendall's Tau-B
    fail_if(evaluate_model(confusion_matrix_1, TAU_B) - 0.70352647 > 1e-6, "TAU_B(1) = 0.70352647");
    fail_if(evaluate_model(confusion_matrix_2, TAU_B) - 0.33333333 > 1e-6, "TAU_B(2) = 0.3333...");
    
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
