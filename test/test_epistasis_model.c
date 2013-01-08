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

#define samples_per_order   (order * num_samples)

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

void print_combination(uint8_t comb[], unsigned long idx, int k) {
    printf("%lu -> {", idx);
    int i;
    for (i = 0; i < k; ++i)
        printf("%d, ", comb[i]);
    printf("\b\b}\n");
}

START_TEST(test_get_masks) {
    int num_masks;
    uint8_t *masks = get_masks(order, genotypes, num_samples, &num_masks);
    
    fail_unless(num_masks == 3 * samples_per_order, "There should be 3 GTs * 4 SNPs * 8 samples = 96 elements in the masks array");
    
    uint8_t mask_0[] = { 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0 };
    uint8_t mask_1[] = { 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0 };
    uint8_t mask_2[] = { 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1 };
    
    print_combination(masks, 0, samples_per_order);
    print_combination(mask_0, 0, samples_per_order);
    print_combination(masks + samples_per_order, 1, samples_per_order);
    print_combination(mask_1, 1, num_samples * 4);
    print_combination(masks + samples_per_order * 2, 2, samples_per_order);
    print_combination(mask_2, 2, num_samples * 4);
    
    for (int j = 0; j < samples_per_order; j++) {
        char msg[23];
        sprintf(msg, "Mask 0: %d should be %d\n", j, mask_0[j]);
        fail_if(masks[j] != mask_0[j], msg);
    }
    
    for (int j = 0; j < samples_per_order; j++) {
        char msg[23];
        sprintf(msg, "Mask 1: %d should be %d\n", j, mask_1[j]);
        fail_if(masks[num_samples * 4 + j] != mask_1[j], msg);
    }
    
    for (int j = 0; j < samples_per_order; j++) {
        char msg[23];
        sprintf(msg, "Mask 2: %d should be %d\n", j, mask_2[j]);
        fail_if(masks[num_samples * 8 + j] != mask_2[j], msg);
    }
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
    
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Epistasis model");
    suite_add_tcase(fs, tc_masks);
    
    return fs;
}
