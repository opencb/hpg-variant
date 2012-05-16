#include <stdio.h>
#include <stdlib.h>

#include <check.h>

#include <vcf_file_structure.h>
#include <vcf_file.h>
#include <vcf_read.h>

#include "../stats/stats.h"


static vcf_record_t *record;
static list_t *output_list;
static list_item_t *record_item;


Suite *create_test_suite(void);


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/

void setup_stats(void) {
    record = create_record();
    record->chromosome = (char*) calloc (2, sizeof(char));
    strcat(record->chromosome, "1");
    
    record->id = (char*) calloc (5, sizeof(char));
    strcat(record->id, "rs12");
    
    record->position = 1234567;
    
    record->reference = (char*) calloc (2, sizeof(char));
    strcat(record->reference, "G");
    
    record->filter = (char*) calloc (5, sizeof(char));
    strcat(record->filter, "PASS");
    
    record->info = (char*) calloc (5, sizeof(char));
    strcat(record->info, "NS=3");
    
    output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", 1, 8, output_list);
    record_item = list_item_new(0, 0, record);
}

void teardown_stats(void) {
//     vcf_record_free(record);
}

/* ******************************
 *          Unit tests         *
 * ******************************/

START_TEST (biallelic) {
    record->alternate = (char*) calloc (2, sizeof(char));
    strcat(record->alternate, "T");
    record->format = (char*) calloc (6, sizeof(char));
    strcat(record->format, "GC:GT");
    
    size_t sample_idx = 0;
    char *sample = (char*) calloc (6, sizeof(char*));
    strcat(sample, "1:0/0");
    add_record_sample(sample, record, &sample_idx);
    sample = (char*) calloc (6, sizeof(char*));
    strcat(sample, "2:1/0");
    add_record_sample(sample, record, &sample_idx);
    sample = (char*) calloc (6, sizeof(char*));
    strcat(sample, "1:0/1");
    add_record_sample(sample, record, &sample_idx);
    sample = (char*) calloc (6, sizeof(char*));
    strcat(sample, "3:0/0");
    add_record_sample(sample, record, &sample_idx);
    sample = (char*) calloc (6, sizeof(char*));
    strcat(sample, "1:1/1");
    add_record_sample(sample, record, &sample_idx);
    sample = (char*) calloc (6, sizeof(char*));
    strcat(sample, "1:./1");
    add_record_sample(sample, record, &sample_idx);
    
    get_variants_stats(record_item, 1, output_list);
    fail_if(output_list->length == 0, "There must be an element processed");
    
    /*
     * typedef struct {
    char *ref_allele;
    char **alternates;
    
    int num_alleles;
    int *alleles_count;
    int *genotypes_count;
    
    int missing_alleles;
    int missing_genotypes;
} variant_stats_t;*/
    
    variant_stats_t *result = (variant_stats_t*) output_list->first_p->data_p;
    
    fail_unless(strcmp(result->ref_allele, "G") == 0, "The reference allele should be G");
    fail_unless(result->num_alleles == 2, "There should be 2 alleles");
    fail_unless(strcmp(result->alternates[0], "T") == 0, "The alternate allele should be T");
    fail_unless(result->alleles_count[0] == 6, "There should be 6 reference alleles read");
    fail_unless(result->alleles_count[1] == 5, "There should be 5 alternate alleles read");
    fail_unless(result->missing_alleles == 1, "There should be 1 missing allele");
    fail_unless(result->missing_genotypes == 1, "There should be 1 missing genotype");
}
END_TEST

START_TEST (multiallelic) {
    record->alternate = (char*) calloc (5, sizeof(char));
    strcat(record->alternate, "T,GT");
    record->format = (char*) calloc (6, sizeof(char));
    strcat(record->format, "GT:GC");
}
END_TEST

START_TEST (homozygous) {
    record->alternate = (char*) calloc (2, sizeof(char));
    strcat(record->alternate, ".");
    record->format = (char*) calloc (12, sizeof(char));
    strcat(record->format, "GT:GQ:DP:HQ");
}
END_TEST

START_TEST (whole_test) {
    
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
    TCase *tc_stats_function = tcase_create("Stats function");
    tcase_add_unchecked_fixture(tc_stats_function, setup_stats, teardown_stats);
    tcase_add_test(tc_stats_function, biallelic);
    tcase_add_test(tc_stats_function, multiallelic);
    tcase_add_test(tc_stats_function, homozygous);
    
    
    TCase *tc_system = tcase_create("System test");
    tcase_add_test(tc_system, whole_test);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("hpg-vcf-tools stats");
    suite_add_tcase(fs, tc_stats_function);
    suite_add_tcase(fs, tc_system);
    
    return fs;
}
