#include <stdio.h>
#include <stdlib.h>

#include <check.h>

#include <bioformats/vcf/vcf_file.h>

#include "../merge/merge.h"


Suite *create_test_suite(void);


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
    add_sample_name("S01", files[0]);
    add_sample_name("S02", files[0]);
    
    files[1] = calloc (1, sizeof(vcf_file_t));
    files[1]->filename = "input1.vcf";
    files[1]->samples_names = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);
    add_sample_name("S11", files[1]);
    add_sample_name("S12", files[1]);
    add_sample_name("S13", files[1]);
    
    files[2] = calloc (1, sizeof(vcf_file_t));
    files[2]->filename = "input2.vcf";
    files[2]->samples_names = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);
    add_sample_name("S21", files[2]);
}

void teardown_merge_process(void) {
    
}


/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST (merge_position_in_one_file) {
    vcf_record_file_link position;
    position.file = files[1];
    
    vcf_record_t *input = create_record();
    input->chromosome = "1";
    input->position = 11111111111;
    input->id = "rs123456";
    input->reference = "A";
    input->alternate = "T";
    input->quality = 20;
    input->filter = "PASS";
    input->info = "NS=3;DP=14;AF=0.5;DB;H2";
    input->format = "GT:GQ:DP:HQ";
    add_record_sample("1/1:20:40:30", input);
    add_record_sample("0/1:10:60:50", input);
    add_record_sample("0/0:30:50:70", input);
    position.record = input;
    
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
}
END_TEST

START_TEST (merge_position_in_several_files) {
    
}
END_TEST

START_TEST (merge_batch) {
    
}
END_TEST

START_TEST (merge_whole_process) {
    
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
    TCase *tc_merge = tcase_create("Merge tool test");
    tcase_add_checked_fixture(tc_merge, setup_merge_process, teardown_merge_process);
    tcase_add_test(tc_merge, merge_position_in_one_file);
//     tcase_add_test(tc_merge, merge_position_in_several_files);
//     tcase_add_test(tc_merge, merge_batch);
//     tcase_add_test(tc_merge, merge_whole_process);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Check for hpg-vcf/merge");
    suite_add_tcase(fs, tc_merge);
    
    return fs;
}
