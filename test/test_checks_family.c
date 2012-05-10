#include <stdio.h>
#include <stdlib.h>

#include <check.h>

#include "../gwas/checks_family.h"


Suite *create_test_suite(void);



/* ******************************
 *      Unchecked fixtures      *
 * ******************************/



/* ******************************
 *          Unit tests         *
 * ******************************/

START_TEST (mendel_valid_families) {
    fail_if(check_mendel("20", 0, 0, 0, 0, 0, 0, MALE) != 0,    "00x00 -> 00 is a valid family");
    
    fail_if(check_mendel("20", 0, 0, 0, 1, 0, 0, FEMALE) != 0,  "00x01 -> 00 is a valid family");
    fail_if(check_mendel("20", 0, 0, 0, 1, 0, 1, FEMALE) != 0,  "00x01 -> 01 is a valid family");
    fail_if(check_mendel("20", 0, 0, 0, 1, 1, 0, FEMALE) != 0,  "00x01 -> 10 is a valid family");
    
    fail_if(check_mendel("20", 0, 0, 1, 0, 0, 0, MALE) != 0,    "00x10 -> 00 is a valid family");
    fail_if(check_mendel("20", 0, 0, 1, 0, 0, 1, MALE) != 0,    "00x10 -> 01 is a valid family");
    fail_if(check_mendel("20", 0, 0, 1, 0, 1, 0, MALE) != 0,    "00x10 -> 10 is a valid family");
    
    fail_if(check_mendel("20", 0, 1, 0, 0, 0, 0, MALE) != 0,    "01x00 -> 00 is a valid family");
    fail_if(check_mendel("20", 0, 1, 0, 0, 0, 1, FEMALE) != 0,  "01x00 -> 01 is a valid family");
    fail_if(check_mendel("20", 0, 1, 0, 0, 1, 0, MALE) != 0,    "01x00 -> 10 is a valid family");
    
    fail_if(check_mendel("20", 0, 1, 0, 1, 0, 0, FEMALE) != 0,  "01x01 -> 00 is a valid family");
    fail_if(check_mendel("20", 0, 1, 0, 1, 0, 1, FEMALE) != 0,  "01x01 -> 01 is a valid family");
    fail_if(check_mendel("20", 0, 1, 0, 1, 1, 0, FEMALE) != 0,  "01x01 -> 10 is a valid family");
    fail_if(check_mendel("20", 0, 1, 0, 1, 1, 1, FEMALE) != 0,  "01x01 -> 11 is a valid family");
    
    fail_if(check_mendel("20", 0, 1, 1, 0, 0, 0, MALE) != 0,    "01x10 -> 00 is a valid family");
    fail_if(check_mendel("20", 0, 1, 1, 0, 0, 1, MALE) != 0,    "01x10 -> 01 is a valid family");
    fail_if(check_mendel("20", 0, 1, 1, 0, 1, 0, MALE) != 0,    "01x10 -> 10 is a valid family");
    fail_if(check_mendel("20", 0, 1, 1, 0, 1, 1, MALE) != 0,    "01x10 -> 11 is a valid family");
    
    fail_if(check_mendel("20", 0, 1, 1, 0, 1, 1, MALE) != 0,    "01x11 -> 01 is a valid family");
    fail_if(check_mendel("20", 0, 1, 1, 0, 1, 1, MALE) != 0,    "01x11 -> 11 is a valid family");
}
END_TEST

START_TEST (mendel_type_1) {
    fail_if(check_mendel("20", 0, 0, 0, 0, 0, 1, MALE) != 1,    "00x00 -> 01 is a mendel error 1");
    fail_if(check_mendel("20", 0, 0, 0, 0, 1, 0, MALE) != 1,    "00x00 -> 10 is a mendel error 1");
}
END_TEST

START_TEST (mendel_type_2) {
    fail_if(check_mendel("20", 1, 1, 1, 1, 0, 1, MALE) != 2,    "11x11 -> 01 is a mendel error 2");
    fail_if(check_mendel("20", 1, 1, 1, 1, 1, 0, MALE) != 2,    "11x11 -> 10 is a mendel error 2");
}
END_TEST

START_TEST (mendel_type_3) {
    
}
END_TEST

START_TEST (mendel_type_4) {
    
}
END_TEST

START_TEST (mendel_type_5) {
    
}
END_TEST

START_TEST (mendel_type_6) {
    
}
END_TEST

START_TEST (mendel_type_7) {
    
}
END_TEST

START_TEST (mendel_type_8) {
    
}
END_TEST

START_TEST (mendel_type_9) {
    
}
END_TEST

START_TEST (mendel_type_10) {
    
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


Suite *create_test_suite(void)
{
    TCase *tc_mendel = tcase_create("Mendelian errors");
    tcase_add_test(tc_mendel, mendel_valid_families);
    tcase_add_test(tc_mendel, mendel_type_1);
    tcase_add_test(tc_mendel, mendel_type_2);
    tcase_add_test(tc_mendel, mendel_type_3);
    tcase_add_test(tc_mendel, mendel_type_4);
    tcase_add_test(tc_mendel, mendel_type_5);
    tcase_add_test(tc_mendel, mendel_type_6);
    tcase_add_test(tc_mendel, mendel_type_7);
    tcase_add_test(tc_mendel, mendel_type_8);
    tcase_add_test(tc_mendel, mendel_type_9);
    tcase_add_test(tc_mendel, mendel_type_10);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Checks for families correctnes");
    suite_add_tcase(fs, tc_mendel);
    
    return fs;
}
