#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <check.h>

#include <bioformats/family/family.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <commons/string_utils.h>
#include <containers/array_list.h>

#include "gwas/gwas.h"
#include "gwas/tdt/tdt_runner.h"

Suite *create_test_suite(void);


static ped_file_t *ped;
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
    ped = ped_open("tdt_files/4K_variants_147_samples.ped");
    family = family_new("TESTFAM");
    add_family(family, ped);
    
    father_sample = (char*) calloc (3, sizeof(char));
    mother_sample = (char*) calloc (3, sizeof(char));
    child_sample = (char*) calloc (3, sizeof(char));
    
    record = vcf_record_new();
    set_vcf_record_chromosome("1", 1, record);
    set_vcf_record_position(111111, record);
    set_vcf_record_reference("C", 1, record);
    set_vcf_record_alternate("T", 1, record);
    set_vcf_record_format("GT", 2, record);
    
    output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", 1, 1, output_list);
}

void teardown_tdt_function(void) {
    free(father_sample);
    free(mother_sample);
    free(child_sample);
}


/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST (family_unaffected_child) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, AFFECTED, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, AFFECTED, NULL, NULL, family);
    child = individual_new("CHILD00", 2.0, MALE, UNAFFECTED, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/0");
    
    array_list_insert(father_sample, record->samples);
    array_list_insert(mother_sample, record->samples);
    array_list_insert(child_sample, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD00", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(&record, 1, (family_t **) cp_hashtable_get_values(ped->families), get_num_families(ped), 
                         sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 0, "With unaffected child, b=0");
    fail_unless(result->t2 == 0, "With unaffected child, c=0");
}
END_TEST

START_TEST (family_01_01_00) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, AFFECTED, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, AFFECTED, NULL, NULL, family);
    child = individual_new("CHILD00", 2.0, MALE, AFFECTED, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/0");
    
    array_list_insert(father_sample, record->samples);
    array_list_insert(mother_sample, record->samples);
    array_list_insert(child_sample, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD00", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(&record, 1, (family_t **) cp_hashtable_get_values(ped->families), get_num_families(ped), 
                         sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 2, "In 01-01->00, b=2");
    fail_unless(result->t2 == 0, "In 01-01->00, c=0");
}
END_TEST

START_TEST (family_01_00_00) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, AFFECTED, NULL, NULL, family);
    mother = individual_new("MOT00", 2.0, FEMALE, AFFECTED, NULL, NULL, family);
    child = individual_new("CHILD00", 2.0, MALE, AFFECTED, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/0");
    strcat(child_sample, "0/0");
    
    array_list_insert(father_sample, record->samples);
    array_list_insert(mother_sample, record->samples);
    array_list_insert(child_sample, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT00", pos1);
    cp_hashtable_put(sample_ids, "CHILD00", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(&record, 1, (family_t **) cp_hashtable_get_values(ped->families), get_num_families(ped), 
                         sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 1, "In 01-00->00, b=1");
    fail_unless(result->t2 == 0, "In 01-00->00, c=0");
}
END_TEST

START_TEST (family_01_01_01) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, AFFECTED, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, AFFECTED, NULL, NULL, family);
    child = individual_new("CHILD01", 2.0, MALE, AFFECTED, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/1");
    
    array_list_insert(father_sample, record->samples);
    array_list_insert(mother_sample, record->samples);
    array_list_insert(child_sample, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD01", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(&record, 1, (family_t **) cp_hashtable_get_values(ped->families), get_num_families(ped), 
                         sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 1, "In 01-01->01, b=1");
    fail_unless(result->t2 == 1, "In 01-01->01, c=1");
}
END_TEST

START_TEST (family_01_00_01) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, AFFECTED, NULL, NULL, family);
    mother = individual_new("MOT00", 2.0, FEMALE, AFFECTED, NULL, NULL, family);
    child = individual_new("CHILD01", 2.0, MALE, AFFECTED, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/0");
    strcat(child_sample, "0/1");
    
    array_list_insert(father_sample, record->samples);
    array_list_insert(mother_sample, record->samples);
    array_list_insert(child_sample, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT00", pos1);
    cp_hashtable_put(sample_ids, "CHILD01", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(&record, 1, (family_t **) cp_hashtable_get_values(ped->families), get_num_families(ped), 
                         sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 0, "In 01-00->01, b=0");
    fail_unless(result->t2 == 1, "In 01-00->01, c=1");
}
END_TEST

START_TEST (family_01_11_01) {
    // Create family
    father = individual_new("FAT01", 2.0, MALE, AFFECTED, NULL, NULL, family);
    mother = individual_new("MOT11", 2.0, FEMALE, AFFECTED, NULL, NULL, family);
    child = individual_new("CHILD01", 2.0, MALE, AFFECTED, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "1/1");
    strcat(child_sample, "0/1");
    
    array_list_insert(father_sample, record->samples);
    array_list_insert(mother_sample, record->samples);
    array_list_insert(child_sample, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT11", pos1);
    cp_hashtable_put(sample_ids, "CHILD01", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(&record, 1, (family_t **) cp_hashtable_get_values(ped->families), get_num_families(ped), 
                         sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 1, "In 01-11->01, b=1");
    fail_unless(result->t2 == 0, "In 01-11->01, c=0");
}
END_TEST

START_TEST (family_00_01_01) {
    // Create family
    father = individual_new("FAT00", 2.0, MALE, AFFECTED, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, AFFECTED, NULL, NULL, family);
    child = individual_new("CHILD01", 2.0, MALE, AFFECTED, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/0");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/1");
    
    array_list_insert(father_sample, record->samples);
    array_list_insert(mother_sample, record->samples);
    array_list_insert(child_sample, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT00", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD01", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(&record, 1, (family_t **) cp_hashtable_get_values(ped->families), get_num_families(ped), 
                         sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 0, "In 00-01->01, b=0");
    fail_unless(result->t2 == 1, "In 00-01->01, c=1");
}
END_TEST

START_TEST (family_11_01_01) {
    // Create family
    father = individual_new("FAT11", 2.0, MALE, AFFECTED, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, AFFECTED, NULL, NULL, family);
    child = individual_new("CHILD01", 2.0, MALE, AFFECTED, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    
    // Create VCF record to insert into variants
    strcat(father_sample, "1/1");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/1");
    
    array_list_insert(father_sample, record->samples);
    array_list_insert(mother_sample, record->samples);
    array_list_insert(child_sample, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(6, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT11", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD01", pos2);
    
    // Launch and verify execution
    fail_unless(tdt_test(&record, 1, (family_t **) cp_hashtable_get_values(ped->families), get_num_families(ped), 
                         sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 1, "In 11-01->01, b=1");
    fail_unless(result->t2 == 0, "In 11-01->01, c=0");
}
END_TEST


START_TEST (combined_families) {
    // TODO combine various situations from the previous test case
    
    // Create family #1
    father = individual_new("FAT01", 2.0, MALE, AFFECTED, NULL, NULL, family);
    mother = individual_new("MOT01", 2.0, FEMALE, AFFECTED, NULL, NULL, family);
    child = individual_new("CHILD00", 2.0, MALE, AFFECTED, father, mother, family);
    family_set_parent(father, family);
    family_set_parent(mother, family);
    family_add_child(child, family);
    // Create family #2
    familyB = family_new("TESTFAM-B");
    add_family(familyB, ped);
    fatherB = individual_new("FAT01B", 2.0, MALE, AFFECTED, NULL, NULL, familyB);
    motherB = individual_new("MOT00B", 2.0, FEMALE, AFFECTED, NULL, NULL, familyB);
    childB = individual_new("CHILD00B", 2.0, MALE, AFFECTED, fatherB, motherB, familyB);
    family_set_parent(fatherB, familyB);
    family_set_parent(motherB, familyB);
    family_add_child(childB, familyB);
    
    
    // Create VCF record to insert into variants
    strcat(father_sample, "0/1");
    strcat(mother_sample, "0/1");
    strcat(child_sample, "0/0");
    
    father_sampleB = (char*) calloc (4, sizeof(char));
    mother_sampleB = (char*) calloc (4, sizeof(char));
    child_sampleB = (char*) calloc (4, sizeof(char));
    strcat(father_sampleB, "0/1");
    strcat(mother_sampleB, "0/0");
    strcat(child_sampleB, "0/0");
    
    array_list_insert(father_sample, record->samples);
    array_list_insert(mother_sample, record->samples);
    array_list_insert(child_sample, record->samples);
    
    array_list_insert(father_sampleB, record->samples);
    array_list_insert(mother_sampleB, record->samples);
    array_list_insert(child_sampleB, record->samples);
    
    // Create ordering structure
    sample_ids = cp_hashtable_create(12, cp_hash_string, (cp_compare_fn) strcasecmp);
    cp_hashtable_put(sample_ids, "FAT01", pos0);
    cp_hashtable_put(sample_ids, "MOT01", pos1);
    cp_hashtable_put(sample_ids, "CHILD00", pos2);
    cp_hashtable_put(sample_ids, "FAT01B", pos3);
    cp_hashtable_put(sample_ids, "MOT00B", pos4);
    cp_hashtable_put(sample_ids, "CHILD00B", pos5);
    
    // Launch and verify execution
    fail_unless(tdt_test(&record, 1, (family_t **) cp_hashtable_get_values(ped->families), get_num_families(ped), 
                         sample_ids, output_list) == 0, "TDT test terminated with errors");
    fail_if(output_list->length == 0, "There must be one result inserted");
    
    tdt_result_t *result = output_list->first_p->data_p;
    fail_unless(result->t1 == 3, "In 01-01->00, b=3");
    fail_unless(result->t2 == 0, "In 01-01->00, c=0");
}
END_TEST


START_TEST (whole_test) {
    // Invoke hpg-variant/genome-analysis --tdt
    int tdt_ret = system("../bin/hpg-var-gwas tdt --vcf-file tdt_files/4K_variants_147_samples.vcf \
                                                  --ped-file tdt_files/4K_variants_147_samples.ped \
                                                  --outdir ./ ");
    fail_unless(tdt_ret == 0, "hpg-var-gwas exited with errors");
    
    // Check results in base to output file: get each line, grep plink.tdt and compare fields
    FILE *cut_proc = popen("cut -f2,5,6 ./hpg-variant.tdt", "r");
    FILE *grep_proc;
    int line_len = 256;
    char line[line_len];
    char grepcmd[line_len];
    char grep_line[line_len];
    char *aux;
    int t1, t2;
    int t1_plink, t2_plink;
    
    fail_if(fgets(&line, line_len, cut_proc) == NULL, "The TDT results file cannot be read");
    
    int i = 0;
    int entries_tested = 2000;
    
    while (fgets(&line, line_len, cut_proc) && i < entries_tested) {
        if (i % 50 == 0) {
            printf("checked record #%d\n", i);
        }
        
//         printf("Line: '%s'\n", line);
        aux = trim(strtok(line, "\t"));
        t1 = atoi(strtok(NULL, "\t"));
        t2 = atoi(strtok(NULL, "\t"));
//         printf("t1 = %d\tt2 = %d\n", t1, t2);
        
        memset(grepcmd, 0, line_len * sizeof(char));
        strncat(grepcmd, "grep -w ", 8);
        strncat(grepcmd, aux, strlen(aux));
//         strncat(grepcmd, " tdt_files/plink_dirty.tdt", 26);
        strncat(grepcmd, " tdt_files/plink.tdt", 20);
        strncat(grepcmd, " | awk '{print $6 \" \" $7}'", 28);
        
//         printf("grepcmd = %s", grepcmd);
        
        grep_proc = popen(grepcmd, "r");
        fgets(&grep_line, line_len, grep_proc);
        
        t1_plink = atoi(strtok(grep_line, " "));
        t2_plink = atoi(strtok(NULL, " "));
        
        if (t1 != t1_plink) { 
            printf("[%d, %s] t1 = %d\tt1_plink = %d\n", i, aux, t1, t1_plink);
        }
        
        fail_if(t1 != t1_plink, "t1 must be the same in both files");
        fail_if(t2 != t2_plink, "t2 must be the same in both files");
        
        pclose(grep_proc);
        i++;
    }
    
    pclose(cut_proc);
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
    tcase_add_test(tc_tdt_test_function, family_unaffected_child);
    tcase_add_test(tc_tdt_test_function, family_00_01_01);
    tcase_add_test(tc_tdt_test_function, family_01_00_00);
    tcase_add_test(tc_tdt_test_function, family_01_00_01);
    tcase_add_test(tc_tdt_test_function, family_01_01_00);
    tcase_add_test(tc_tdt_test_function, family_01_01_01);
    tcase_add_test(tc_tdt_test_function, family_01_11_01);
    tcase_add_test(tc_tdt_test_function, family_11_01_01);
    tcase_add_test(tc_tdt_test_function, combined_families);
    
    TCase *tc_pipeline = tcase_create("Pipeline integration");
    tcase_add_test(tc_pipeline, whole_test);
    tcase_set_timeout(tc_pipeline, 0);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("TDT test");
    suite_add_tcase(fs, tc_tdt_test_function);
    suite_add_tcase(fs, tc_pipeline);
    
    return fs;
}

