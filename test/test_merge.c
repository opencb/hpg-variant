#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <check.h>

#include <bioformats/vcf/vcf_file.h>

#include "merge/merge.h"


Suite *create_test_suite(void);

vcf_record_t *create_example_record_0();
vcf_record_t *create_example_record_1();
vcf_record_t *create_example_record_2();
vcf_record_t *create_example_record_3();

vcf_file_t *files[3];
merge_options_data_t *options;


/* ******************************
 *       Checked fixtures       *
 * ******************************/

void setup_merge_process(void) {
    options = calloc(1, sizeof(merge_options_data_t));
    options->missing_mode = MISSING;
    
    files[0] = calloc (1, sizeof(vcf_file_t));
    files[0]->filename = "input0.vcf";
    files[0]->samples_names = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);
    add_vcf_sample_name("S01", 3, files[0]);
    add_vcf_sample_name("S02", 3, files[0]);
    
    files[1] = calloc (1, sizeof(vcf_file_t));
    files[1]->filename = "input1.vcf";
    files[1]->samples_names = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);
    add_vcf_sample_name("S11", 3, files[1]);
    add_vcf_sample_name("S12", 3, files[1]);
    add_vcf_sample_name("S13", 3, files[1]);
    
    files[2] = calloc (1, sizeof(vcf_file_t));
    files[2]->filename = "input2.vcf";
    files[2]->samples_names = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);
    add_vcf_sample_name("S21", 3, files[2]);
}

void teardown_merge_process(void) {
    
}


/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST (merge_position_in_one_file) {
    vcf_record_t *input;
    vcf_record_file_link position;
    position.file = files[1];
    position.record = input = create_example_record_0();
    
    vcf_record_t *result = merge_unique_position(&position, files, 3, options);
    
    fail_if(result == NULL, "A VCF record must be returned after merging");
    fail_if(strcmp(result->chromosome, input->chromosome), "After merging a position present only in a file, the chromosome must stay the same");
    fail_if(result->position != input->position, "After merging a position present only in a file, the position must stay the same");
    fail_if(strcmp(result->id, input->id), "After merging a position present only in a file, the ID must stay the same");
    fail_if(strcmp(result->reference, input->reference), "After merging a position present only in a file, the reference must stay the same");
    fail_if(strcmp(result->alternate, input->alternate), "After merging a position present only in a file, the alternate must stay the same");
    fail_if(result->quality != input->quality, "After merging a position present only in a file, the quality must stay the same");
    fail_if(strcmp(result->filter, input->filter), "After merging a position present only in a file, the filter must stay the same");
    fail_if(strcmp(result->format, input->format), "After merging a position present only in a file, the format must stay the same");
    
    fail_if(result->samples->size != 6, "There must be 6 samples in the resulting record");
    fail_if(strcmp(result->samples->items[0], "./.:.:.:."), "Sample 0 must be empty");
    fail_if(strcmp(result->samples->items[1], "./.:.:.:."), "Sample 1 must be empty");
    fail_if(strcmp(result->samples->items[2], "1/1:20:40:30"), "Sample 2 must be 1/1:20:40:30");
    fail_if(strcmp(result->samples->items[3], "0/1:10:60:50"), "Sample 3 must be 0/1:10:60:50");
    fail_if(strcmp(result->samples->items[4], "0/0:30:50:70"), "Sample 4 must be 0/0:30:50:70");
    fail_if(strcmp(result->samples->items[5], "./.:.:.:."), "Sample 5 must be empty");
    
    write_vcf_record(input, stdout);
}
END_TEST

START_TEST (merge_id_test) {
    vcf_record_t *input[4];
    input[0] = create_example_record_0();
    input[1] = create_example_record_1();
    input[2] = create_example_record_2();
    input[3] = create_example_record_3();
    
    // Merge (rs123456, ".", rs654321) -> rs123456
    vcf_record_file_link **links = calloc (3, sizeof(vcf_record_file_link*));
    for (int i = 0; i < 3; i++) {
        links[i] = malloc(sizeof(vcf_record_file_link));
        links[i]->file = files[i];
        links[i]->record = input[i];
        assert(links[i]->record->id);
    }
    fail_if(strcmp(merge_id_field(links, 3), "rs123456"), "After merging (0,1,2), the ID must be rs123456");
    
    // Merge (".", rs654321, ".") -> rs654321
    for (int i = 0; i < 3; i++) {
        links[i]->record = input[i+1];
        assert(links[i]->record->id);
    }
    fail_if(strcmp(merge_id_field(links, 3), "rs654321"), "After merging (1,2,3), the ID must be rs654321");
    
    // Merge (".", ".") -> "."
    links[0]->record = input[1];
    links[1]->record = input[3];
    
    char *id = merge_id_field(links, 2);
    fail_if(strcmp(merge_id_field(links, 2), "."), "After merging (1,3), the ID must be '.'");
}
END_TEST

START_TEST (merge_alternate_test) {
    
}
END_TEST

START_TEST (merge_quality_test) {
    
}
END_TEST

START_TEST (merge_filter_test) {
    
}
END_TEST

START_TEST (merge_info_test) {
    
}
END_TEST

START_TEST (merge_format_test) {
    
}
END_TEST

START_TEST (merge_samples_test) {
    
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
    TCase *tc_merge = tcase_create("Merge position in one file");
    tcase_add_checked_fixture(tc_merge, setup_merge_process, teardown_merge_process);
    tcase_add_test(tc_merge, merge_position_in_one_file);
    
    TCase *tc_repeated = tcase_create("Merge position in several files");
    tcase_add_checked_fixture(tc_repeated, setup_merge_process, teardown_merge_process);
    tcase_add_test(tc_merge, merge_id_test);
    tcase_add_test(tc_merge, merge_alternate_test);
    tcase_add_test(tc_merge, merge_quality_test);
    tcase_add_test(tc_merge, merge_filter_test);
    tcase_add_test(tc_merge, merge_info_test);
    tcase_add_test(tc_merge, merge_format_test);
    tcase_add_test(tc_merge, merge_samples_test);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Check for hpg-vcf/merge");
    suite_add_tcase(fs, tc_merge);
    suite_add_tcase(fs, tc_repeated);
    
    return fs;
}


/* ******************************
 *      Auxiliary functions     *
 * ******************************/

vcf_record_t *create_example_record_0() {
    vcf_record_t *input = vcf_record_new();
    input->chromosome = "1";
    input->chromosome_len = strlen(input->chromosome);
    input->position = 21111111111;
    input->id = "rs123456";
    input->id_len = strlen(input->id);
    input->reference = "A";
    input->reference_len = strlen(input->reference);
    input->alternate = "T";
    input->alternate_len = strlen(input->alternate);
    input->quality = 20;
    input->filter = "PASS";
    input->filter_len = strlen(input->filter);
    input->info = "NS=3;DP=14;H2";
    input->info_len = strlen(input->info);
    input->format = "GT:GQ:DP:HQ";
    input->format_len = strlen(input->format);
    add_vcf_record_sample("1/1:20:40:30", 12, input);
    add_vcf_record_sample("0/1:10:60:50", 12, input);
    add_vcf_record_sample("0/0:30:50:70", 12, input);
    
    return input;
}

vcf_record_t *create_example_record_1() {
    vcf_record_t *input = vcf_record_new();
    input->chromosome = "1";
    input->chromosome_len = strlen(input->chromosome);
    input->position = 21111111111;
    input->id = ".";
    input->id_len = strlen(input->id);
    input->reference = "A";
    input->reference_len = strlen(input->reference);
    input->alternate = "G";
    input->alternate_len = strlen(input->alternate);
    input->quality = 30;
    input->filter = ".";
    input->filter_len = strlen(input->filter);
    input->info = "DP=10;NS=4;AF=0.5;H2";
    input->info_len = strlen(input->info);
    input->format = "GT:DP";
    input->format_len = strlen(input->format);
    add_vcf_record_sample("1/1:40", 6, input);
    add_vcf_record_sample("0/1:60", 6, input);
    add_vcf_record_sample("0/0:50", 6, input);
    
    return input;
}

vcf_record_t *create_example_record_2() {
    vcf_record_t *input = vcf_record_new();
    input->chromosome = "1";
    input->chromosome_len = strlen(input->chromosome);
    input->position = 21111111111;
    input->id = "rs654321";
    input->id_len = strlen(input->id);
    input->reference = "A";
    input->reference_len = strlen(input->reference);
    input->alternate = "T";
    input->alternate_len = strlen(input->alternate);
    input->quality = 20;
    input->filter = "PASS";
    input->filter_len = strlen(input->filter);
    input->info = "AF=0.5;NS=3;DP=14;DB;H2";
    input->info_len = strlen(input->info);
    input->format = "DP:HQ:GT:GQ";
    input->format_len = strlen(input->format);
    add_vcf_record_sample("1/1:20:40:30", 12, input);
    add_vcf_record_sample("0/1:10:60:50", 12, input);
    add_vcf_record_sample("0/0:30:50:70", 12, input);
    
    return input;
}

vcf_record_t *create_example_record_3() {
    vcf_record_t *input = vcf_record_new();
    input->chromosome = "1";
    input->chromosome_len = strlen(input->chromosome);
    input->position = 21111111111;
    input->id = ".";
    input->id_len = strlen(input->id);
    input->reference = "A";
    input->reference_len = strlen(input->reference);
    input->alternate = "T";
    input->alternate_len = strlen(input->alternate);
    input->quality = 20;
    input->filter = "PASS";
    input->filter_len = strlen(input->filter);
    input->info = "DB;H2";
    input->info_len = strlen(input->info);
    input->format = "GT";
    input->format_len = strlen(input->format);
    add_vcf_record_sample("1/1", 3, input);
    add_vcf_record_sample("0/1", 3, input);
    add_vcf_record_sample("0/0", 3, input);
    
    return input;
}
