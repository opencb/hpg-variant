#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <check.h>

#include <bioformats/vcf/vcf_file_structure.h>

#include "epistasis/dataset.h"


Suite *create_test_suite(void);

int num_samples = 20;
int phenotypes[20];

int num_records = 3;
char *format = "GT";
char *possible_gts[3] = { "0/0", "0/1", "1/1" };
vcf_record_t **records;

epistasis_dataset *dataset;


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/

void setup_dataset(void) {
    // Create phenotypes
    for (int i = 0; i < num_samples; i++) {
        phenotypes[i] = i % 2;
    }
    
    // Create records
    records = malloc (num_records * sizeof(vcf_record_t*));
    
    // Record 0 gets samples { 0/0, 0/1, 1/1, 0/0, ... }
    records[0] = vcf_record_new();
    set_vcf_record_chromosome("1", 1, records[0]);
    set_vcf_record_position(1, records[0]);
    set_vcf_record_format(format, 2, records[0]);
    for (int i = 0; i < num_samples - 1; i++) {
        add_vcf_record_sample(possible_gts[i % 3], 3, records[0]);
    }
    add_vcf_record_sample("./.", 3, records[0]);
    
    // Record 1 gets samples { 0/1, 1/1, 0/0, 0/1, ... }
    records[1] = vcf_record_new();
    set_vcf_record_chromosome("1", 1, records[1]);
    set_vcf_record_position(2, records[1]);
    set_vcf_record_format(format, 2, records[1]);
    for (int i = 0; i < num_samples; i++) {
        add_vcf_record_sample(possible_gts[(i + 1) % 3], 3, records[1]);
    }
    // Record 2 gets samples { 1/1, 0/0, 0/1, 1/1, ... }
    records[2] = vcf_record_new();
    set_vcf_record_chromosome("1", 1, records[2]);
    set_vcf_record_position(3, records[2]);
    set_vcf_record_format(format, 2, records[2]);
    for (int i = 0; i < num_samples; i++) {
        add_vcf_record_sample(possible_gts[(i + 2) % 3], 3, records[2]);
    }
    
    // Create dataset
    dataset = epistasis_dataset_new();
}

void teardown_dataset(void) {
    epistasis_dataset_free(dataset);
    free(records);
}


/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST (test_process_records) {
    epistasis_dataset_process_records(records, num_records, num_samples, phenotypes, dataset);
    
    fail_unless(epistasis_dataset_get_num_variants(dataset) == 3, "The number of variants in the dataset must be 3");
    
    uint8_t *gt0 = (uint8_t*) array_list_get(0, dataset->genotype_counts);
    uint8_t *gt1 = (uint8_t*) array_list_get(1, dataset->genotype_counts);
    uint8_t *gt2 = (uint8_t*) array_list_get(2, dataset->genotype_counts);
    
    fail_unless(gt0[0] == 0, "Variant 1:1: Sample 0 must have genotype 0/0");
    fail_unless(gt0[1] == 1, "Variant 1:1, Sample 1 must have genotype 0/1");
    fail_unless(gt0[2] == 2, "Variant 1:1, Sample 2 must have genotype 1/1");
    fail_unless(gt0[3] == 0, "Variant 1:1, Sample 3 must have genotype 0/0");
    fail_unless(gt0[4] == 1, "Variant 1:1, Sample 4 must have genotype 0/1");
    fail_unless(gt0[5] == 2, "Variant 1:1, Sample 5 must have genotype 1/1");
    fail_unless(gt0[6] == 0, "Variant 1:1: Sample 6 must have genotype 0/0");
    fail_unless(gt0[7] == 1, "Variant 1:1, Sample 7 must have genotype 0/1");
    fail_unless(gt0[8] == 2, "Variant 1:1, Sample 8 must have genotype 1/1");
    fail_unless(gt0[9] == 0, "Variant 1:1, Sample 9 must have genotype 0/0");
    fail_unless(gt0[10] == 1, "Variant 1:1, Sample 10 must have genotype 0/1");
    fail_unless(gt0[11] == 2, "Variant 1:1, Sample 11 must have genotype 1/1");
    fail_unless(gt0[18] == 0, "Variant 1:1, Sample 18 must have genotype 0/0");
    fail_unless(gt0[19] == 255, "Variant 1:1, Sample 19 must have genotype ./.");
    
    for (int i = 0; i < num_samples - 2; i++) {
        fail_unless(gt1[i] == gt0[i+1], "Variant 1:2: Sample must have the next genotype to the one in first variant");
    }
    fail_unless(gt1[19] == 2, "Variant 1:2: Sample 19 must have genotype 1/1");
    
    for (int i = 0; i < num_samples - 3; i++) {
        fail_unless(gt2[i] == gt0[i+2], "Variant 1:3: Sample must have the next genotype to the one in first variant");
    }
    fail_unless(gt2[18] == 2, "Variant 1:3: Sample 18 must have genotype 1/1");
    fail_unless(gt2[19] == 0, "Variant 1:3: Sample 19 must have genotype 0/0");
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
    TCase *tc_creation = tcase_create("Dataset creation");
    tcase_add_unchecked_fixture(tc_creation, setup_dataset, teardown_dataset);
    tcase_add_test(tc_creation, test_process_records);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Epistasis dataset");
    suite_add_tcase(fs, tc_creation);
    
    return fs;
}
