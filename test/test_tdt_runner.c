#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <check.h>

#include <family.h>
#include <vcf_batch.h>
#include <vcf_file_structure.h>
#include <vcf_read.h>

#include "../gwas/gwas.h"
#include "../gwas/tdt_runner.h"

Suite *create_test_suite(void);


static ped_file_t *ped;

static list_item_t *variant;
static vcf_record_t *record;

static cp_hashtable *sample_ids;
static list_t *output_list;

static family_t *family, *familyB;
static individual_t *father, *mother, *child;
static individual_t *fatherB, *motherB, *childB;
static char *father_sample, *mother_sample, *child_sample;
static char *father_sampleB, *mother_sampleB, *child_sampleB;

static int *pos0, *pos1, *pos2;
static int *pos3, *pos4, *pos5;


/* ******************************
 *       Unchecked fixtures     *
 * ******************************/

void setup_positions(void) {
    pos0 = (int*) malloc (sizeof(int)); *pos0 = 0;
    pos1 = (int*) malloc (sizeof(int)); *pos1 = 1;
    pos2 = (int*) malloc (sizeof(int)); *pos2 = 2;
    pos3 = (int*) malloc (sizeof(int)); *pos3 = 3;
    pos4 = (int*) malloc (sizeof(int)); *pos4 = 4;
    pos5 = (int*) malloc (sizeof(int)); *pos5 = 5;
}

void teardown_positions(void) {
    free(pos0);
    free(pos1);
    free(pos2);
    free(pos3);
    free(pos4);
    free(pos5);
}

/* ******************************
 *        Checked fixtures      *
 * ******************************/



void setup_tdt_function(void) {
    ped = ped_open("tdt_files/newjob.ped");
    family = family_new("TESTFAM");
    add_family(family, ped);
    
    father_sample = (char*) calloc (3, sizeof(char));
    mother_sample = (char*) calloc (3, sizeof(char));
    child_sample = (char*) calloc (3, sizeof(char));
    
    record = create_record();
    set_record_chromosome("1", record);
    set_record_position(111111, record);
    set_record_reference("C", record);
    set_record_alternate("T", record);
    variant = list_item_new(1, 0, record);
    
    output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", 1, 1, output_list);
}

void teardown_tdt_function(void) {
    free(father_sample);
    free(mother_sample);
    free(child_sample);
}

void setup_pipeline(void) {
    
}

void teardown_pipeline(void) {
    
}

/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST (family_01_01_00) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, NULL, NULL, family);
    child = individual_new("CHILD00", 2.0, MALE, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/0");
    
    list_item_t *item = list_item_new(1, 0, father_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(2, 0, mother_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(3, 0, child_sample);
    list_insert_item(item, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD00", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(ped, variant, 1, sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 2, "In 01-01->00, b=2");
    fail_unless(result->t2 == 0, "In 01-01->00, c=0");
}
END_TEST

START_TEST (family_01_00_00) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, NULL, NULL, family);
    mother = individual_new("MOT00", 2.0, FEMALE, NULL, NULL, family);
    child = individual_new("CHILD00", 2.0, MALE, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/0");
    strcat(child_sample, "0/0");
    
    list_item_t *item = list_item_new(1, 0, father_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(2, 0, mother_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(3, 0, child_sample);
    list_insert_item(item, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT00", pos1);
    cp_hashtable_put(sample_ids, "CHILD00", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(ped, variant, 1, sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 1, "In 01-00->00, b=1");
    fail_unless(result->t2 == 0, "In 01-00->00, c=0");
}
END_TEST

START_TEST (family_01_01_01) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, NULL, NULL, family);
    child = individual_new("CHILD01", 2.0, MALE, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/1");
    
    list_item_t *item = list_item_new(1, 0, father_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(2, 0, mother_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(3, 0, child_sample);
    list_insert_item(item, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD01", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(ped, variant, 1, sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 1, "In 01-01->01, b=1");
    fail_unless(result->t2 == 1, "In 01-01->01, c=1");
}
END_TEST

START_TEST (family_01_00_01) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, NULL, NULL, family);
    mother = individual_new("MOT00", 2.0, FEMALE, NULL, NULL, family);
    child = individual_new("CHILD01", 2.0, MALE, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/0");
    strcat(child_sample, "0/1");
    
    list_item_t *item = list_item_new(1, 0, father_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(2, 0, mother_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(3, 0, child_sample);
    list_insert_item(item, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT00", pos1);
    cp_hashtable_put(sample_ids, "CHILD01", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(ped, variant, 1, sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 0, "In 01-00->01, b=0");
    fail_unless(result->t2 == 1, "In 01-00->01, c=1");
}
END_TEST

START_TEST (family_01_11_01) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, NULL, NULL, family);
    mother = individual_new("MOT11", 2.0, FEMALE, NULL, NULL, family);
    child = individual_new("CHILD01", 2.0, MALE, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "1/1");
    strcat(child_sample, "0/1");
    
    list_item_t *item = list_item_new(1, 0, father_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(2, 0, mother_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(3, 0, child_sample);
    list_insert_item(item, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT11", pos1);
    cp_hashtable_put(sample_ids, "CHILD01", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(ped, variant, 1, sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 1, "In 01-11->01, b=1");
    fail_unless(result->t2 == 0, "In 01-11->01, c=0");
}
END_TEST

START_TEST (family_00_01_01) {
    // Create family
    father = individual_new("FAT00", 2.0, MALE, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, NULL, NULL, family);
    child = individual_new("CHILD01", 2.0, MALE, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/0");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/1");
    
    list_item_t *item = list_item_new(1, 0, father_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(2, 0, mother_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(3, 0, child_sample);
    list_insert_item(item, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT00", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD01", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(ped, variant, 1, sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 0, "In 00-01->01, b=0");
    fail_unless(result->t2 == 1, "In 00-01->01, c=1");
}
END_TEST

START_TEST (family_11_01_01) {
    // Create family
    father = individual_new("FAT11", 2.0, MALE, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, NULL, NULL, family);
    child = individual_new("CHILD01", 2.0, MALE, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "1/1");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/1");
    
    list_item_t *item = list_item_new(1, 0, father_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(2, 0, mother_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(3, 0, child_sample);
    list_insert_item(item, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT11", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD01", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(ped, variant, 1, sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 1, "In 11-01->01, b=1");
    fail_unless(result->t2 == 0, "In 11-01->01, c=0");
}
END_TEST


START_TEST (combined_families) {
    // TODO combine various situations from the previous test case
    
    // Create family #1
    father = individual_new("FAT01", 2.0, MALE, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, NULL, NULL, family);
    child = individual_new("CHILD00", 2.0, MALE, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    // Create family #2
    familyB = family_new("TESTFAM-B");
    add_family(familyB, ped);
    fatherB = individual_new("FAT01B", 2.0, MALE, NULL, NULL, familyB);
    motherB = individual_new("MOT00B", 2.0, FEMALE, NULL, NULL, familyB);
    childB = individual_new("CHILD00B", 2.0, MALE, fatherB, motherB, familyB);
    family_set_parent(fatherB, familyB);
    family_set_parent(motherB, familyB);
    family_add_child(childB, familyB);
    
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/0");
    
    father_sampleB = (char*) calloc (3, sizeof(char));
    mother_sampleB = (char*) calloc (3, sizeof(char));
    child_sampleB = (char*) calloc (3, sizeof(char));
    strcat(father_sampleB, "0/1");
    strcat(mother_sampleB, "0/0");
    strcat(child_sampleB, "0/0");
    
    list_item_t *item = list_item_new(1, 0, father_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(2, 0, mother_sample);
    list_insert_item(item, record->samples);
    item = list_item_new(3, 0, child_sample);
    list_insert_item(item, record->samples);
    
    item = list_item_new(4, 0, father_sampleB);
    list_insert_item(item, record->samples);
    item = list_item_new(5, 0, mother_sampleB);
    list_insert_item(item, record->samples);
    item = list_item_new(6, 0, child_sampleB);
    list_insert_item(item, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(12, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD00", pos2);
    cp_hashtable_put(sample_ids, "FAT01B", pos3);
    cp_hashtable_put(sample_ids, "MOT00B", pos4);
    cp_hashtable_put(sample_ids, "CHILD00B", pos5);
    
    // Launch and verify execution
    fail_unless(tdt_test(ped, variant, 1, sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 3, "In 01-01->00, b=3");
    fail_unless(result->t2 == 0, "In 01-01->00, c=0");
}
END_TEST


START_TEST (whole_test) {
    // TODO invoke run_tdt_test and check results in base to output file: 
    // get each line, grep plink.tdt and compare fields
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
    TCase *tc_tdt_test_function = tcase_create("TDT test function");
    tcase_add_unchecked_fixture(tc_tdt_test_function, setup_positions, teardown_positions);
    tcase_add_checked_fixture(tc_tdt_test_function, setup_tdt_function, teardown_tdt_function);
    tcase_add_test(tc_tdt_test_function, family_00_01_01);
    tcase_add_test(tc_tdt_test_function, family_01_00_00);
    tcase_add_test(tc_tdt_test_function, family_01_00_01);
    tcase_add_test(tc_tdt_test_function, family_01_01_00);
    tcase_add_test(tc_tdt_test_function, family_01_01_01);
    tcase_add_test(tc_tdt_test_function, family_01_11_01);
    tcase_add_test(tc_tdt_test_function, family_11_01_01);
    tcase_add_test(tc_tdt_test_function, combined_families);
    
//     TCase *tc_pipeline = tcase_create("Pipeline integration");
//     tcase_add_checked_fixture(tc_pipeline, setup_pipeline, teardown_pipeline);
//     tcase_add_test(tc_pipeline, whole_test);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("TDT test");
    suite_add_tcase(fs, tc_tdt_test_function);
//     suite_add_tcase(fs, tc_pipeline);
    
    return fs;
}

