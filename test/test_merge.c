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

vcf_file_t *files[4];
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
    add_vcf_sample_name("S03", 3, files[0]);
    
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
    
    files[3] = calloc (1, sizeof(vcf_file_t));
    files[3]->filename = "input3.vcf";
    files[3]->samples_names = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);
    add_vcf_sample_name("S31", 3, files[3]);
    add_vcf_sample_name("S32", 3, files[3]);
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
    fail_if(strncmp(input->chromosome, result->chromosome, input->chromosome_len), "After merging a position present only in a file, the chromosome must stay the same");
    fail_if(result->position != input->position, "After merging a position present only in a file, the position must stay the same");
    fail_if(strncmp(input->id, result->id, input->id_len), "After merging a position present only in a file, the ID must stay the same");
    fail_if(strncmp(input->reference, result->reference, input->reference_len), "After merging a position present only in a file, the reference must stay the same");
    fail_if(strncmp(input->alternate, result->alternate, input->alternate_len), "After merging a position present only in a file, the alternate must stay the same");
    fail_if(result->quality != input->quality, "After merging a position present only in a file, the quality must stay the same");
    fail_if(strncmp(input->filter, result->filter, input->filter_len), "After merging a position present only in a file, the filter must stay the same");
    fail_if(strncmp(input->format, result->format, input->format_len), "After merging a position present only in a file, the format must stay the same");
    
    fail_if(result->samples->size != 7, "There must be 7 samples in the resulting record");
    fail_if(strcmp(result->samples->items[0], "./.:.:.:."), "Sample 0 must be empty");
    fail_if(strcmp(result->samples->items[1], "./.:.:.:."), "Sample 1 must be empty");
    fail_if(strcmp(result->samples->items[2], "./.:.:.:."), "Sample 2 must be empty");
    fail_if(strcmp(result->samples->items[3], "1/1:20:40:30"), "Sample 3 must be 1/1:20:40:30");
    fail_if(strcmp(result->samples->items[4], "0/1:10:60:50"), "Sample 4 must be 0/1:10:60:50");
    fail_if(strcmp(result->samples->items[5], "0/0:30:50:70"), "Sample 5 must be 0/0:30:50:70");
    fail_if(strcmp(result->samples->items[6], "./.:.:.:."), "Sample 6 must be empty");
    
    write_vcf_record(input, stdout);
    write_vcf_record(result, stdout);
}
END_TEST

START_TEST (merge_id_test) {
    vcf_record_t *input[4];
    input[0] = create_example_record_0();
    input[1] = create_example_record_1();
    input[2] = create_example_record_2();
    input[3] = create_example_record_3();
    
    vcf_record_file_link **links = calloc (3, sizeof(vcf_record_file_link*));
    for (int i = 0; i < 3; i++) {
        links[i] = malloc(sizeof(vcf_record_file_link));
        links[i]->file = files[i];
    }
    
    // Merge (rs123456, ".", rs654321) -> rs123456
    for (int i = 0; i < 3; i++) {
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
    vcf_record_t *input[4];
    input[0] = create_example_record_0();
    input[1] = create_example_record_1();
    input[2] = create_example_record_2();
    input[3] = create_example_record_3();
    
    vcf_record_file_link **links = calloc (4, sizeof(vcf_record_file_link*));
    for (int i = 0; i < 4; i++) {
        links[i] = malloc(sizeof(vcf_record_file_link));
        links[i]->file = files[i];
        links[i]->record = input[i];
    }
    
    cp_hashtable *alleles_table = cp_hashtable_create(8, cp_hash_istring, (cp_compare_fn) strcasecmp);
    
    // Merge (0,1,2,3) = (T,G,CT,T) -> "T,G,CT"
    fail_if(strcmp(merge_alternate_field(links, 4, alleles_table), "T,G,CT"), "After merging (0,1,2,3), the alternate must be 'T,G,CT'");
    cp_hashtable_remove_all(alleles_table);
    
    // Merge (0,1,2) = (T,G,CT) -> "T,G,CT"
    fail_if(strcmp(merge_alternate_field(links, 3, alleles_table), "T,G,CT"), "After merging (0,1,2), the alternate must be 'T,G,CT'");
    cp_hashtable_remove_all(alleles_table);
    
    // Merge (0,1,3) = (T,G,T) -> "T,G"
    links[0]->record = input[0];
    links[1]->record = input[1];
    links[2]->record = input[3];
    fail_if(strcmp(merge_alternate_field(links, 3, alleles_table), "T,G"), "After merging (0,1,3), the alternate must be 'T,G'");
    cp_hashtable_remove_all(alleles_table);
    
    // Merge (1,2,3) = (G,CT,T) -> "G,CT,T"
    links[0]->record = input[1];
    links[1]->record = input[2];
    links[2]->record = input[3];
    fail_if(strcmp(merge_alternate_field(links, 3, alleles_table), "G,CT,T"), "After merging (1,2,3), the alternate must be 'G,CT,T'");
    cp_hashtable_remove_all(alleles_table);
    
    // Merge (0,3) = (T,T) -> "T"
    links[0]->record = input[0];
    links[1]->record = input[3];
    fail_if(strcmp(merge_alternate_field(links, 2, alleles_table), "T"), "After merging (0,3), the alternate must be 'T'");
    cp_hashtable_remove_all(alleles_table);
}
END_TEST

START_TEST (merge_quality_test) {
    vcf_record_t *input[4];
    input[0] = create_example_record_0();
    input[1] = create_example_record_1();
    input[2] = create_example_record_2();
    input[3] = create_example_record_3();
    
    vcf_record_file_link **links = calloc (4, sizeof(vcf_record_file_link*));
    for (int i = 0; i < 4; i++) {
        links[i] = malloc(sizeof(vcf_record_file_link));
        links[i]->file = files[i];
        links[i]->record = input[i];
    }
    
    // Merge (0,1,2,3) = (20*3,30*3,10*1,.*2) -> ~17.778
    fail_if(merge_quality_field(links, 4) - 17.778 > 0.1, "After merging (0,1,2,3), the quality must be ~17.778");
    
    // Merge (0,1,2) = (20*3,30*3,10*1) -> ~22.857
    fail_if(merge_quality_field(links, 3) - 22.857 > 0.1, "After merging (0,1,2), the quality must be ~22.857");
    
    // Merge (0,1) = (20*3,30*3) -> 25
    fail_if(merge_quality_field(links, 2) - 25 > 0.1, "After merging (0,1), the quality must be 25");
    
    // Merge (0,1,3) = (20*3,30*3,.*2) -> 18.75
    links[0]->record = input[0];
    links[1]->record = input[1];
    links[2]->file   = files[3];
    links[2]->record = input[3];
    fail_if(merge_quality_field(links, 3) - 18.75 > 0.1, "After merging (0,1,3), the quality must be 18.75");
    
    // Merge (0,3) = (20*3,.*2) -> 12
    links[0]->record = input[0];
    links[1]->file   = files[3];
    links[1]->record = input[3];
    fail_if(merge_quality_field(links, 2) - 12 > 0.1, "After merging (0,3), the quality must be 12");
}
END_TEST

START_TEST (merge_filter_test) {
    vcf_record_t *input[4];
    input[0] = create_example_record_0();
    input[1] = create_example_record_1();
    input[2] = create_example_record_2();
    input[3] = create_example_record_3();
    
    vcf_record_file_link **links = calloc (4, sizeof(vcf_record_file_link*));
    for (int i = 0; i < 4; i++) {
        links[i] = malloc(sizeof(vcf_record_file_link));
        links[i]->file = files[i];
        links[i]->record = input[i];
    }
    
    // Merge (0,1,2,3) = STD_FILTER,q10
    fail_if(strcmp(merge_filter_field(links, 4), "STD_FILTER,q10"), "After merging (0,1,2,3), the filter must be STD_FILTER,q10");
    
    // Merge (0,1,3) = STD_FILTER
    links[0]->record = input[0];
    links[1]->record = input[1];
    links[2]->record = input[3];
    fail_if(strcmp(merge_filter_field(links, 3), "STD_FILTER"), "After merging (0,1,3), the filter must be STD_FILTER");
    
    // Merge (0,3) = PASS
    links[0]->record = input[0];
    links[1]->record = input[3];
    fail_if(strcmp(merge_filter_field(links, 2), "PASS"), "After merging (0,3), the filter must be PASS");
    
    // Merge (1,2) = STD_FILTER,q10
    links[0]->record = input[1];
    links[1]->record = input[2];
    fail_if(strcmp(merge_filter_field(links, 2), "STD_FILTER,q10"), "After merging (1,2), the filter must be STD_FILTER,q10");
    
    // Merge (2,1) = q10,STD_FILTER
    links[0]->record = input[2];
    links[1]->record = input[1];
    fail_if(strcmp(merge_filter_field(links, 2), "q10,STD_FILTER"), "After merging (2,1), the filter must be q10,STD_FILTER");
}
END_TEST

START_TEST (merge_info_test) {
    
}
END_TEST

START_TEST (merge_format_test) {
    vcf_record_t *input[4];
    input[0] = create_example_record_0();
    input[1] = create_example_record_1();
    input[2] = create_example_record_2();
    input[3] = create_example_record_3();
    
    vcf_record_file_link **links = calloc (4, sizeof(vcf_record_file_link*));
    for (int i = 0; i < 4; i++) {
        links[i] = malloc(sizeof(vcf_record_file_link));
        links[i]->file = files[i];
        links[i]->record = input[i];
    }
    
    array_list_t *format_fields = array_list_new(8, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    format_fields->compare_fn = strcasecmp;
    
    // Merge (0,1,2,3) -> GT:GQ:DP:HQ:RD
    fail_if(strcmp(merge_format_field(links, 4, format_fields), "GT:GQ:DP:HQ:RD"), "After merging (0,1,2,3), the format must be GT:GQ:DP:HQ:RD");
    fail_if(strcmp(format_fields->items[0], "GT"), "After merging (0,1,2,3), field #0 in format must be GT");
    fail_if(strcmp(format_fields->items[1], "GQ"), "After merging (0,1,2,3), field #1 in format must be GQ");
    fail_if(strcmp(format_fields->items[2], "DP"), "After merging (0,1,2,3), field #2 in format must be DP");
    fail_if(strcmp(format_fields->items[3], "HQ"), "After merging (0,1,2,3), field #3 in format must be HQ");
    fail_if(strcmp(format_fields->items[4], "RD"), "After merging (0,1,2,3), field #4 in format must be RD");
    array_list_clear(format_fields, free);
    
    // Merge (1,2,3) -> GT:RD:HQ:GQ
    links[0]->record = input[1];
    links[1]->record = input[2];
    links[2]->record = input[3];
    fail_if(strcmp(merge_format_field(links, 3, format_fields), "GT:RD:HQ:GQ"), "After merging (1,2,3), the format must be GT:RD:HQ:GQ");
    fail_if(strcmp(format_fields->items[0], "GT"), "After merging (1,2,3), field #0 in format must be GT");
    fail_if(strcmp(format_fields->items[1], "RD"), "After merging (1,2,3), field #1 in format must be RD");
    fail_if(strcmp(format_fields->items[2], "HQ"), "After merging (1,2,3), field #2 in format must be HQ");
    fail_if(strcmp(format_fields->items[3], "GQ"), "After merging (1,2,3), field #3 in format must be GQ");
    array_list_clear(format_fields, free);
    
    // Merge (1,3) -> GT:RD
    links[0]->record = input[1];
    links[1]->record = input[3];
    fail_if(strcmp(merge_format_field(links, 2, format_fields), "GT:RD"), "After merging (1,3), the format must be GT:RD");
    fail_if(strcmp(format_fields->items[0], "GT"), "After merging (1,3), field #0 in format must be GT");
    fail_if(strcmp(format_fields->items[1], "RD"), "After merging (1,3), field #1 in format must be RD");
    array_list_clear(format_fields, free);
    
    // Merge (3,1) -> GT:RD
    links[0]->record = input[3];
    links[1]->record = input[1];
    fail_if(strcmp(merge_format_field(links, 2, format_fields), "GT:RD"), "After merging (3,1), the format must be GT:RD");
    fail_if(strcmp(format_fields->items[0], "GT"), "After merging (3,1), field #0 in format must be GT");
    fail_if(strcmp(format_fields->items[1], "RD"), "After merging (3,1), field #1 in format must be RD");
    array_list_clear(format_fields, free);
    
    // Merge (2,1) -> RD:HQ:GT:GQ
    links[0]->record = input[2];
    links[1]->record = input[1];
    fail_if(strcmp(merge_format_field(links, 2, format_fields), "RD:HQ:GT:GQ"), "After merging (2,1), the format must be RD:HQ:GT:GQ");
    fail_if(strcmp(format_fields->items[0], "RD"), "After merging (2,1), field #0 in format must be RD");
    fail_if(strcmp(format_fields->items[1], "HQ"), "After merging (2,1), field #1 in format must be HQ");
    fail_if(strcmp(format_fields->items[2], "GT"), "After merging (2,1), field #2 in format must be GT");
    fail_if(strcmp(format_fields->items[3], "GQ"), "After merging (2,1), field #3 in format must be GQ");
    array_list_clear(format_fields, free);
    
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
    input->filter = "STD_FILTER";
    input->filter_len = strlen(input->filter);
    input->info = "DP=10;NS=4;AF=0.5;H2";
    input->info_len = strlen(input->info);
    input->format = "GT:RD";
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
    input->alternate = "CT";
    input->alternate_len = strlen(input->alternate);
    input->quality = 10;
    input->filter = "q10";
    input->filter_len = strlen(input->filter);
    input->info = "AF=0.5;NS=3;DP=14;DB;H2";
    input->info_len = strlen(input->info);
    input->format = "RD:HQ:GT:GQ";
    input->format_len = strlen(input->format);
    add_vcf_record_sample("20:40:1/1:30", 12, input);
    
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
    input->quality = -1;
    input->filter = ".";
    input->filter_len = strlen(input->filter);
    input->info = "DB;H2";
    input->info_len = strlen(input->info);
    input->format = "GT";
    input->format_len = strlen(input->format);
    add_vcf_record_sample("1/1", 3, input);
    add_vcf_record_sample("0/1", 3, input);
    
    return input;
}
