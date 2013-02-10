#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <check.h>

#include <containers/array_list.h>

#include "effect/effect.h"
#include "effect/effect_runner.h"


Suite *create_test_suite(void);


shared_options_data_t *global_data;
effect_options_data_t *opts_data;


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/

void setup_effect_ws(void) {
    global_data = (shared_options_data_t*) calloc (1, sizeof(shared_options_data_t));
    global_data->vcf_filename = "effect_files/variants_marta_head_500.vcf";
    global_data->output_directory = "/tmp/variant-test/";
    global_data->num_threads = 4;
    global_data->max_batches = 10;
    global_data->batch_lines = 4000;
    global_data->entries_per_thread = 1000;
    
    opts_data = (effect_options_data_t*) calloc (1, sizeof(effect_options_data_t));
    
    initialize_ws_output(global_data, opts_data);
}

void teardown_effect_ws(void) {
//     free_ws_output(4);
}


/* ******************************
 *          Unit tests         *
 * ******************************/

START_TEST (url_composition) {
    // All null arguments
    fail_if(compose_effect_ws_request("genomic/variant", "consequence_type", global_data), 
            "The resulting URL must be null (all args are null)");
    
    // Some not-null arguments
    global_data->species = "hsa";
    global_data->version = "latest";
    fail_if(compose_effect_ws_request("genomic/variant", "consequence_type", global_data), 
            "The resulting URL must be null (some args are null)");
    
    // None null argument
    global_data->host_url = "http://localhost:8080";
    global_data->species = "hsa";
    global_data->version = "latest";
    
    char *url = compose_effect_ws_request("genomic/variant", "consequence_type", global_data);
    fail_if(strcmp(url, "http://localhost:8080/cellbase/rest/latest/hsa/genomic/variant/consequence_type?header=false"),
            "The resulting URL must be 'http://localhost:8080/cellbase/rest/latest/hsa/genomic/variant/consequence_type'"); 
}
END_TEST


START_TEST (effect_ws_request) {
    char *url = "http://mem16:8080/cellbase/rest/latest/hsa/genomic/variant/consequence_type?header=false";
    
    vcf_record_t *record_1 = (vcf_record_t*) malloc (sizeof(vcf_record_t));
    record_1->chromosome = "11";
    record_1->position = 423423;
    record_1->reference = "G";
    record_1->alternate = "C";
    
    vcf_record_t *record_2 = (vcf_record_t*) malloc (sizeof(vcf_record_t));
    record_2->chromosome = "11";
    record_2->position = 4234230;
    record_2->reference = "AC";
    record_2->alternate = "C";
    
    vcf_record_t *record_3 = (vcf_record_t*) malloc (sizeof(vcf_record_t));
    record_3->chromosome = "12";
    record_3->position = 4234230;
    record_3->reference = "A";
    record_3->alternate = "TC";
    
    vcf_record_t *records[3];
    records[0] = record_1;
    records[1] = record_2;
    records[2] = record_3;
    
    fail_if(invoke_effect_ws(url, records, 3, ""), "The web service request was not successfully performed");
}
END_TEST


START_TEST (whole_test) {
    char buf[100];
    int i = 0;
    FILE *p;
    
    // Invoke hpg-variant/effect
    int tdt_ret = system("../bin/hpg-var-effect --vcf-file effect_files/variants_marta_head_500.vcf \
                                                --outdir ./ --region 1:13000-13500");
    fail_unless(tdt_ret == 0, "hpg-variant exited with errors");
    
    // Check all variants
    p = popen("/usr/bin/wc -l all_variants.txt","r");
    if (p) {
        while (!feof(p) && (i < 99) ) {
            fread(&buf[i],1,1,p);
            i++;
        }
        buf[i] = 0;
        fail_unless(atoi(strtok(buf, "\t")) == 100, "There must be 100 entries in all_variants.txt");
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
        fail_unless(atoi(strtok(buf, "\t")) == 21, "There must be 34 entries in 5KB_downstream_variant.txt");
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
        fail_unless(atoi(strtok(buf, "\t")) == 3, "There must be 3 entries in intron_variant.txt");
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
        fail_unless(atoi(strtok(buf, "\t")) == 45, "There must be 45 entries in regulatory_region_variant.txt");
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
        fail_unless(atoi(strtok(buf, "\t")) == 8, "There must be 8 entries in summary.txt");
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
    tcase_add_unchecked_fixture(tc_composition, setup_effect_ws, teardown_effect_ws);
    tcase_add_test(tc_composition, url_composition);
    
    TCase *tc_web_service = tcase_create("Web service request");
    tcase_add_unchecked_fixture(tc_web_service, setup_effect_ws, teardown_effect_ws);
    tcase_add_test(tc_web_service, effect_ws_request);
    tcase_set_timeout(tc_web_service, 0);
    
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
