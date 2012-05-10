#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <check.h>

#include "../effect/effect.h"
#include "../effect/effect_runner.h"


Suite *create_test_suite(void);


effect_options_data_t opts_data;


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/

void setup_effect_ws(void) {
    initialize_ws_output(4, "./");
}

void teardown_effect_ws(void) {
    free_ws_output(4);
}


/* ******************************
 *          Unit tests         *
 * ******************************/

START_TEST (url_composition) {
    // All null arguments
    fail_unless(compose_effect_ws_request(&opts_data) == NULL, "The resulting URL must be null (all args are null)");
    
    // Some not-null arguments
    opts_data.species = "hsa";
    opts_data.version = "v1";
    fail_unless(compose_effect_ws_request(&opts_data) == NULL, "The resulting URL must be null (some args are null)");
    
    // None null argument
    opts_data.host_url = "http://localhost:8080";
    opts_data.species = "hsa";
    opts_data.version = "v1";
    
    char *url = compose_effect_ws_request(&opts_data);
    fail_if(strcmp(url, "http://localhost:8080/cellbase/rest/v1/hsa/genomic/variant/consequence_type"),
            "The resulting URL must be 'http://localhost:8080/cellbase/rest/v1/hsa/genomic/variant/consequence_type'"); 
}
END_TEST


START_TEST (effect_ws_request) {
    char *url = "http://localhost:8080/cellbase/rest/v1/hsa/genomic/variant/consequence_type";
    vcf_batch_t *batch = (vcf_batch_t*) malloc (sizeof(vcf_batch_t));
    list_item_t *item = NULL;
    
    list_init("test batch", 1, 3, batch);
    
    vcf_record_t *record_1 = (vcf_record_t*) malloc (sizeof(vcf_record_t));
    record_1->chromosome = "11";
    record_1->position = 423423;
    record_1->reference = "G";
    record_1->alternate = "C";
    item = list_item_new(1, 1, record_1);
    list_insert_item(item, batch);
    
    vcf_record_t *record_2 = (vcf_record_t*) malloc (sizeof(vcf_record_t));
    record_2->chromosome = "11";
    record_2->position = 4234230;
    record_2->reference = "AC";
    record_2->alternate = "C";
    item = list_item_new(2, 1, record_2);
    list_insert_item(item, batch);
    
    vcf_record_t *record_3 = (vcf_record_t*) malloc (sizeof(vcf_record_t));
    record_3->chromosome = "12";
    record_3->position = 4234230;
    record_3->reference = "A";
    record_3->alternate = "TC";
    item = list_item_new(3, 1, record_3);
    list_insert_item(item, batch);
    
    fail_unless(invoke_effect_ws(url, batch->first_p, 4) == 0, "The web service request was not successfully performed");
    
    write_summary_file();
}
END_TEST


START_TEST (effect_ws_response) {
    char buf[100];
    int i = 0;
    FILE *p;
    
    // Check summary
    p = popen("/usr/bin/wc -l summary.txt","r");
    if (p) {
        i = 0;
        while (!feof(p) && (i < 99) ) {
            fread(&buf[i],1,1,p);
            i++;
        }
        buf[i] = 0;
        fail_unless(atoi(strtok(buf, "\t")) == 3, "There must be 3 entries in summary.txt");
        pclose(p);
    }
}
END_TEST


START_TEST (whole_test) {
    char buf[100];
    int i = 0;
    FILE *p;
    
    // Invoke hpg-variant/effect
    int tdt_ret = system("../hpg-variant effect --vcf-file effect_files/variants_marta_head_3K.vcf --outdir ./ --region 1:10000-400000");
    fail_unless(tdt_ret == 0, "hpg-variant exited with errors");
    
    // Check all variants
    p = popen("/usr/bin/wc -l all_variants.txt","r");
    if (p) {
        while (!feof(p) && (i < 99) ) {
            fread(&buf[i],1,1,p);
            i++;
        }
        buf[i] = 0;
        fail_unless(atoi(strtok(buf, "\t")) == 499, "There must be 499 entries in all_variants.txt");
        pclose(p);
    }
    
    // Check 5KB downstream
    p = popen("/usr/bin/wc -l 5KB_downstream_variant.txt","r");
    if (p) {
        i = 0;
        while (!feof(p) && (i < 99) ) {
            fread(&buf[i],1,1,p);
            i++;
        }
        buf[i] = 0;
        fail_unless(atoi(strtok(buf, "\t")) == 34, "There must be 34 entries in 5KB_downstream_variant.txt");
        pclose(p);
    }
    
    // Check intron
    p = popen("/usr/bin/wc -l intron_variant.txt","r");
    if (p) {
        i = 0;
        while (!feof(p) && (i < 99) ) {
            fread(&buf[i],1,1,p);
            i++;
        }
        buf[i] = 0;
        fail_unless(atoi(strtok(buf, "\t")) == 47, "There must be 47 entries in intron_variant.txt");
        pclose(p);
    }
    
    // Check regulatory region
    p = popen("/usr/bin/wc -l regulatory_region_variant.txt","r");
    if (p) {
        i = 0;
        while (!feof(p) && (i < 99) ) {
            fread(&buf[i],1,1,p);
            i++;
        }
        buf[i] = 0;
        fail_unless(atoi(strtok(buf, "\t")) == 256, "There must be 256 entries in regulatory_region_variant.txt");
        pclose(p);
    }
    
    // Check summary
    p = popen("/usr/bin/wc -l summary.txt","r");
    if (p) {
        i = 0;
        while (!feof(p) && (i < 99) ) {
            fread(&buf[i],1,1,p);
            i++;
        }
        buf[i] = 0;
        fail_unless(atoi(strtok(buf, "\t")) == 13, "There must be 13 entries in summary.txt");
        pclose(p);
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


Suite *create_test_suite(void)
{
    TCase *tc_composition = tcase_create("URL composition");
    tcase_add_test(tc_composition, url_composition);
    
    TCase *tc_web_service = tcase_create("Web service request");
    tcase_add_unchecked_fixture(tc_web_service, setup_effect_ws, teardown_effect_ws);
    tcase_add_test(tc_web_service, effect_ws_request);
    tcase_add_test(tc_web_service, effect_ws_response);
    
    TCase *tc_system = tcase_create("System test");
    tcase_add_test(tc_system, whole_test);
    tcase_set_timeout(tc_system, 0);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Effect (consequence type) web service");
    suite_add_tcase(fs, tc_composition);
    suite_add_tcase(fs, tc_web_service);
    suite_add_tcase(fs, tc_system);
    
    return fs;
}
