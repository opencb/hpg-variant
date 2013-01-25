#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <check.h>

#include <containers/array_list.h>

#include "epistasis/mdr.h"
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
    }
    
    // 400 { 150 aff, 250 unaff} in 4-fold
    k = 4;
    folds = get_k_folds(150, 250, 4, &sizes);
    i = 0; fail_unless (sizes[3 * i] == 101 && sizes[3 * i + 1] == 38 && sizes[3 * i + 2] == 63, "Fold #0 size must be (101,38,63)");
    i = 1; fail_unless (sizes[3 * i] == 101 && sizes[3 * i + 1] == 38 && sizes[3 * i + 2] == 63, "Fold #1 size must be (101,38,63)");
    i = 2; fail_unless (sizes[3 * i] == 100 && sizes[3 * i + 1] == 37 && sizes[3 * i + 2] == 63, "Fold #2 size must be (100,37,63)");
    i = 3; fail_unless (sizes[3 * i] == 98 && sizes[3 * i + 1] == 37 && sizes[3 * i + 2] == 61, "Fold #3 size must be (98,37,61)");
    
    // 400 { 150 aff, 250 unaff} in 10-fold
    k = 10;
    folds = get_k_folds(150, 250, k, &sizes);
    for (i = 0; i < k; i++) {
        fail_unless (sizes[3 * i] == 40 && sizes[3 * i + 1] == 15 && sizes[3 * i + 2] == 25, "Fold size must be (40,15,25)");
    }
    
    
    // 125 { 50 aff, 75 unaff} in 5-fold
    k = 5;
    folds = get_k_folds(50, 75, k, &sizes);
    for (i = 0; i < k; i++) {
        fail_unless (sizes[3 * i] == 25 && sizes[3 * i + 1] == 10 && sizes[3 * i + 2] == 15, "Fold size must be (25,10,15)");
    }
    
    // 125 { 50 aff, 75 unaff} in 7-fold
    k = 7;
    folds = get_k_folds(50, 75, k, &sizes);
    i = 0; fail_unless (sizes[3 * i] == 19 && sizes[3 * i + 1] == 8 && sizes[3 * i + 2] == 11, "Fold #0 size must be (19,8,11)");
    for (i = 1; i < k-1; i++) {
        fail_unless (sizes[3 * i] == 18 && sizes[3 * i + 1] == 7 && sizes[3 * i + 2] == 11, "Fold size must be (19,7,11)");
    }
    i = k-1; fail_unless (sizes[3 * i] == 16 && sizes[3 * i + 1] == 7 && sizes[3 * i + 2] == 9, "Fold #0 size must be (16,7,9)");
    
    // 125 { 50 aff, 75 unaff} in 10-fold
    k = 10;
    folds = get_k_folds(50, 75, k, &sizes);
    for (i = 0; i < 5; i++) {
        fail_unless (sizes[3 * i] == 13 && sizes[3 * i + 1] == 5 && sizes[3 * i + 2] == 8, "Fold size must be (13,5,8)");
    }
    for (i = 5; i < k; i++) {
        fail_unless (sizes[3 * i] == 12 && sizes[3 * i + 1] == 5 && sizes[3 * i + 2] == 7, "Fold size must be (12,5,7)");
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
    
    
    // 20 { 16 aff, 4 unaff) in 5-fold
    k = 5;
    folds = get_k_folds(16, 4, k, &sizes);
    i = 0; fail_unless (sizes[3 * i] == 5 && sizes[3 * i + 1] == 4 && sizes[3 * i + 2] == 1, "Fold #0 size must be (5,4,1)");
    i = 1; fail_unless (sizes[3 * i] == 4 && sizes[3 * i + 1] == 3 && sizes[3 * i + 2] == 1, "Fold #1 size must be (4,3,1)");
    i = 2; fail_unless (sizes[3 * i] == 4 && sizes[3 * i + 1] == 3 && sizes[3 * i + 2] == 1, "Fold #2 size must be (4,3,1)");
    i = 3; fail_unless (sizes[3 * i] == 4 && sizes[3 * i + 1] == 3 && sizes[3 * i + 2] == 1, "Fold #3 size must be (4,3,1)");
    i = 4; fail_unless (sizes[3 * i] == 3 && sizes[3 * i + 1] == 3 && sizes[3 * i + 2] == 0, "Fold #4 size must be (3,3,0)");
    
    // 20 { 16 aff, 4 unaff) in 10-fold
    k = 10;
    folds = get_k_folds(16, 4, k, &sizes);
    for (i = 0; i < 4; i++) {
        fail_unless (sizes[3 * i] == 3 && sizes[3 * i + 1] == 2 && sizes[3 * i + 2] == 1, "Fold size must be (3,2,1)");
    }
    for (i = 4; i < 8; i++) {
        fail_unless (sizes[3 * i] == 2 && sizes[3 * i + 1] == 2 && sizes[3 * i + 2] == 0, "Fold size must be (2,2,0)");
    }
    for (i = 9; i < k; i++) {
        fail_unless (sizes[3 * i] == 0 && sizes[3 * i + 1] == 0 && sizes[3 * i + 2] == 0, "Fold size must be (0,0,0)");
    }
    
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
    
    TCase *tc_counts = tcase_create("Risk classification");
    tcase_add_test(tc_counts, test_get_high_risk_combinations);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("MDR");
    suite_add_tcase(fs, tc_k_fold);
    suite_add_tcase(fs, tc_counts);
    
    return fs;
}
