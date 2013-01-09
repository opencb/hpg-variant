#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <check.h>

#include "epistasis/model.h"


Suite *create_test_suite(void);

int order = 4;
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
    int num_counts;
    
    // Order 2
    int* counts = get_counts(2, genotypes, num_affected, num_unaffected, &num_counts);
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
    
    // Order 3
    counts = get_counts(3, genotypes, num_affected, num_unaffected, &num_counts);
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
    TCase *tc_masks = tcase_create("Genotype masks");
    tcase_add_test(tc_masks, test_get_masks);
    
    TCase *tc_counts = tcase_create("Genotype counts");
    tcase_add_test(tc_masks, test_get_counts);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Epistasis model");
    suite_add_tcase(fs, tc_masks);
    suite_add_tcase(fs, tc_counts);
    
    return fs;
}
