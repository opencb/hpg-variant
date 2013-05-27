#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <check.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>

#include "vcf-tools/merge/merge.c"


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

void setup_merge_headers(void) {
    options = calloc(1, sizeof(merge_options_data_t));
    options->missing_mode = MISSING;
    options->num_info_fields = 4;
    options->info_fields = malloc (options->num_info_fields * sizeof(char*));
    options->info_fields[0] = "DP";
    options->info_fields[1] = "AN";
    options->info_fields[2] = "MQ";
    options->info_fields[3] = "H3";
    options->config_search_paths = array_list_new(2, 2, COLLECTION_MODE_ASYNCHRONIZED);
    array_list_insert(".", options->config_search_paths);
    
    vcf_header_entry_t *entry;
    char *value;
    
    /* ************* File 0 **************/
    
    files[0] = vcf_file_new("input0.vcf", INT_MAX);
    
    // GT:GQ:DP:HQ
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    assert(add_vcf_header_entry(entry, files[0]));
    
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[0]);
    
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[0]);
    
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=HQ,Number=1,Type=String,Description=\"Haplotype Qualities\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[0]);
    
    // All filters PASS. Example filters: somefilter, bias1
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FILTER", 6, entry);
    value = "<ID=somefilter,Description=\"FS > 200.0\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[0]);
    
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FILTER", 6, entry);
    value = "<ID=bias1,Description=\"Non-biased read\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[0]);
    
    // Should be ignored
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("INFO", 4, entry);
    value = "<ID=PL,Number=.,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[0]);
    
    // Other headers
    entry = vcf_header_entry_new();
    value = "Example for a test of VCF merging";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[0]);
    
    
    /* ************* File 1 **************/
    
    files[1] = vcf_file_new("input1.vcf", INT_MAX);
    
    // GT:RD
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[1]);
    
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=RD,Number=1,Type=Integer,Description=\"Read Depth\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[1]);
    
    // Filter STD_FILTER
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FILTER", 6, entry);
    value = "<ID=STD_FILTER,Description=\"QD < 2.0\", \"MQ < 40.0\", \"FS > 60.0\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[1]);
    
    // Other headers
    entry = vcf_header_entry_new();
    value = "analysis_type=CombineVariants input_file=[] read_buffer_size=null phone_home=STANDARD read_filter=[]";
    set_vcf_header_entry_name("CombineVariants", strlen("CombineVariants"), entry);
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[1]);
    
    
    /* ************* File 2 **************/
    
    files[2] = vcf_file_new("input2.vcf", INT_MAX);
    
    // RD:HQ:GT:GQ
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=RD,Number=1,Type=Integer,Description=\"Read Depth\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[2]);
    
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=HQ,Number=1,Type=String,Description=\"Haplotype Qualities\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[2]);
    
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[2]);
    
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[2]);
    
    // Filter q10
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FILTER", 6, entry);
    value = "<ID=q10,Description=\"QUAL > 10\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[2]);
    
    // Should be ignored
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("INFO", 4, entry);
    value = "<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[2]);
    
    
    /* ************* File 3 **************/
    
    files[3] = vcf_file_new("input3.vcf", INT_MAX);
    
    // GT
    entry = vcf_header_entry_new();
    set_vcf_header_entry_name("FORMAT", 6, entry);
    value = "<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    add_vcf_header_entry_value(value, strlen(value), entry);
    add_vcf_header_entry(entry, files[3]);
}

void setup_merge_positions(void) {
    options = calloc(1, sizeof(merge_options_data_t));
    options->missing_mode = MISSING_CONDITION;
    options->num_info_fields = 1;
    options->info_fields = malloc (sizeof(char*));
    options->info_fields[0] = "DP";
    options->copy_filter = 0;
    options->copy_info = 0;
    options->config_search_paths = array_list_new(2, 2, COLLECTION_MODE_ASYNCHRONIZED);
    array_list_insert(".", options->config_search_paths);
    
    files[0] = vcf_file_new("input0.vcf", INT_MAX);
    add_vcf_sample_name("S01", 3, files[0]);
    add_vcf_sample_name("S02", 3, files[0]);
    add_vcf_sample_name("S03", 3, files[0]);
    
    files[1] = vcf_file_new("input1.vcf", INT_MAX);
    add_vcf_sample_name("S11", 3, files[1]);
    add_vcf_sample_name("S12", 3, files[1]);
    add_vcf_sample_name("S13", 3, files[1]);
    
    files[2] = vcf_file_new("input2.vcf", INT_MAX);
    add_vcf_sample_name("S21", 3, files[2]);
    
    files[3] = vcf_file_new("input3.vcf", INT_MAX);
    add_vcf_sample_name("S31", 3, files[3]);
    add_vcf_sample_name("S32", 3, files[3]);
}

void teardown_merge_headers(void) { }

void teardown_merge_positions(void) { }



/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST (merge_headers_test) {
    list_t *output_list = malloc (sizeof(list_t));
    list_init("headers", 1, INT_MAX, output_list);
    merge_vcf_headers(files, 4, options, output_list);
    
    vcf_header_entry_t *entry;
    
    fail_unless(output_list->length == 16, "There must be 16 header entries");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("FORMAT", entry->name), "The name of entry #0 must be FORMAT");
    fail_if(strcmp("<ID=GT,Number=1,Type=String,Description=\"Genotype\">", array_list_get(0, entry->values)), 
            "The value of entry #0 must be <ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("FORMAT", entry->name), "The name of entry #1 must be FORMAT");
    fail_if(strcmp("<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">", array_list_get(0, entry->values)), 
            "The value of entry #1 must be <ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("FORMAT", entry->name), "The name of entry #2 must be FORMAT");
    fail_if(strcmp("<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", array_list_get(0, entry->values)), 
            "The value of entry #2 must be <ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("FORMAT", entry->name), "The name of entry #3 must be FORMAT");
    fail_if(strcmp("<ID=HQ,Number=1,Type=String,Description=\"Haplotype Qualities\">", array_list_get(0, entry->values)), 
            "The value of entry #3 must be <ID=HQ,Number=1,Type=String,Description=\"Haplotype Qualities\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("FILTER", entry->name), "The name of entry #4 must be FILTER");
    fail_if(strcmp("<ID=somefilter,Description=\"FS > 200.0\">", array_list_get(0, entry->values)), 
            "The value of entry #4 must be <ID=somefilter,Description=\"FS > 200.0\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("FILTER", entry->name), "The name of entry #5 must be FILTER");
    fail_if(strcmp("<ID=bias1,Description=\"Non-biased read\">", array_list_get(0, entry->values)), 
            "The value of entry #5 must be <ID=bias1,Description=\"Non-biased read\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("INFO", entry->name), "The name of entry #6 must be INFO");
    fail_if(strcmp("<ID=PL,Number=.,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes\">", array_list_get(0, entry->values)), 
            "The value of entry #6 must be <ID=PL,Number=.,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(entry->name, "Entry #7 must have no name");
    fail_if(strcmp("Example for a test of VCF merging", array_list_get(0, entry->values)), 
            "The value of entry #7 must be 'Example for a test of VCF merging'");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("FORMAT", entry->name), "The name of entry #8 must be FORMAT");
    fail_if(strcmp("<ID=RD,Number=1,Type=Integer,Description=\"Read Depth\">", array_list_get(0, entry->values)), 
            "The value of entry #8 must be <ID=RD,Number=1,Type=Integer,Description=\"Read Depth\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("FILTER", entry->name), "The name of entry #9 must be FILTER");
    fail_if(strcmp("<ID=STD_FILTER,Description=\"QD < 2.0\", \"MQ < 40.0\", \"FS > 60.0\">", array_list_get(0, entry->values)), 
            "The value of entry #9 must be <ID=STD_FILTER,Description=\"QD < 2.0\", \"MQ < 40.0\", \"FS > 60.0\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("CombineVariants", entry->name), "The name of entry #10 must be CombineVariants");
    fail_if(strcmp("analysis_type=CombineVariants input_file=[] read_buffer_size=null phone_home=STANDARD read_filter=[]", array_list_get(0, entry->values)), 
            "The value of entry #10 must be 'analysis_type=CombineVariants input_file=[] read_buffer_size=null phone_home=STANDARD read_filter=[]'");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("FILTER", entry->name), "The name of entry #11 must be FILTER");
    fail_if(strcmp("<ID=q10,Description=\"QUAL > 10\">", array_list_get(0, entry->values)), 
            "The value of entry #11 must be <ID=q10,Description=\"QUAL > 10\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("INFO", entry->name), "The name of entry #12 must be INFO");
    fail_if(strcmp("<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">", array_list_get(0, entry->values)), 
            "The value of entry #12 must be <ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("INFO", entry->name), "The name of entry #13 must be INFO");
    fail_if(strcmp("<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", array_list_get(0, entry->values)), 
            "The value of entry #13 must be <ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("INFO", entry->name), "The name of entry #14 must be INFO");
    fail_if(strcmp("<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">", array_list_get(0, entry->values)), 
            "The value of entry #14 must be <ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    
    entry = list_remove_item(output_list)->data_p;
    fail_if(strcmp("INFO", entry->name), "The name of entry #15 must be INFO");
    fail_if(strcmp("<ID=H3,Number=0,Type=Flag,Description=\"HapMap3 Membership\">", array_list_get(0, entry->values)), 
            "The value of entry #15 must be <ID=H3,Number=0,Type=Flag,Description=\"HapMap3 Membership\">");
    
    fail_if(output_list->length > 0, "The list must be empty");
}
END_TEST


START_TEST (merge_position_in_one_file) {
    vcf_record_file_link **links = calloc (1, sizeof(vcf_record_file_link*));
    links[0] = malloc(sizeof(vcf_record_file_link));
    links[0]->file = files[1];
    links[0]->record = create_example_record_0();
    
    int err_code;
    vcf_record_t *input = links[0]->record;
    vcf_record_t *result = merge_position(links, 1, files, 3, options, &err_code);
    
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
    
    // Merge (0,1,2,3) = STD_FILTER;q10
    fail_if(strcmp(merge_filter_field(links, 4), "STD_FILTER;q10"), "After merging (0,1,2,3), the filter must be STD_FILTER,q10");
    
    // Merge (0,1,3) = STD_FILTER
    links[0]->record = input[0];
    links[1]->record = input[1];
    links[2]->record = input[3];
    fail_if(strcmp(merge_filter_field(links, 3), "STD_FILTER"), "After merging (0,1,3), the filter must be STD_FILTER");
    
    // Merge (0,3) = PASS
    links[0]->record = input[0];
    links[1]->record = input[3];
    fail_if(strcmp(merge_filter_field(links, 2), "PASS"), "After merging (0,3), the filter must be PASS");
    
    // Merge (1,2) = STD_FILTER;q10
    links[0]->record = input[1];
    links[1]->record = input[2];
    fail_if(strcmp(merge_filter_field(links, 2), "STD_FILTER;q10"), "After merging (1,2), the filter must be STD_FILTER,q10");
    
    // Merge (2,1) = q10;STD_FILTER
    links[0]->record = input[2];
    links[1]->record = input[1];
    fail_if(strcmp(merge_filter_field(links, 2), "q10;STD_FILTER"), "After merging (2,1), the filter must be q10,STD_FILTER");
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
    fail_if(strcmp(merge_format_field(links, 4, options, format_fields), "GT:GQ:DP:HQ:RD"), "After merging (0,1,2,3), the format must be GT:GQ:DP:HQ:RD");
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
    fail_if(strcmp(merge_format_field(links, 3, options, format_fields), "GT:RD:HQ:GQ"), "After merging (1,2,3), the format must be GT:RD:HQ:GQ");
    fail_if(strcmp(format_fields->items[0], "GT"), "After merging (1,2,3), field #0 in format must be GT");
    fail_if(strcmp(format_fields->items[1], "RD"), "After merging (1,2,3), field #1 in format must be RD");
    fail_if(strcmp(format_fields->items[2], "HQ"), "After merging (1,2,3), field #2 in format must be HQ");
    fail_if(strcmp(format_fields->items[3], "GQ"), "After merging (1,2,3), field #3 in format must be GQ");
    array_list_clear(format_fields, free);
    
    // Merge (1,3) -> GT:RD
    links[0]->record = input[1];
    links[1]->record = input[3];
    fail_if(strcmp(merge_format_field(links, 2, options, format_fields), "GT:RD"), "After merging (1,3), the format must be GT:RD");
    fail_if(strcmp(format_fields->items[0], "GT"), "After merging (1,3), field #0 in format must be GT");
    fail_if(strcmp(format_fields->items[1], "RD"), "After merging (1,3), field #1 in format must be RD");
    array_list_clear(format_fields, free);
    
    // Merge (3,1) -> GT:RD
    links[0]->record = input[3];
    links[1]->record = input[1];
    fail_if(strcmp(merge_format_field(links, 2, options, format_fields), "GT:RD"), "After merging (3,1), the format must be GT:RD");
    fail_if(strcmp(format_fields->items[0], "GT"), "After merging (3,1), field #0 in format must be GT");
    fail_if(strcmp(format_fields->items[1], "RD"), "After merging (3,1), field #1 in format must be RD");
    array_list_clear(format_fields, free);
    
    // Merge (2,1) -> RD:HQ:GT:GQ
    links[0]->record = input[2];
    links[1]->record = input[1];
    fail_if(strcmp(merge_format_field(links, 2, options, format_fields), "RD:HQ:GT:GQ"), "After merging (2,1), the format must be RD:HQ:GT:GQ");
    fail_if(strcmp(format_fields->items[0], "RD"), "After merging (2,1), field #0 in format must be RD");
    fail_if(strcmp(format_fields->items[1], "HQ"), "After merging (2,1), field #1 in format must be HQ");
    fail_if(strcmp(format_fields->items[2], "GT"), "After merging (2,1), field #2 in format must be GT");
    fail_if(strcmp(format_fields->items[3], "GQ"), "After merging (2,1), field #3 in format must be GQ");
    array_list_clear(format_fields, free);
    
}
END_TEST

START_TEST (merge_samples_test) {
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
    
    printf("The coordinates of incorrect samples will be presented as (file,index in file)\n");
    
    // Merge samples of position in all files
    cp_hashtable *alleles_table = cp_hashtable_create(8, cp_hash_istring, (cp_compare_fn) strcasecmp);
    merge_alternate_field(links, 4, alleles_table);
    
    array_list_t *format_fields = array_list_new(8, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    format_fields->compare_fn = strcasecmp;
    char *format = merge_format_field(links, 4, options, format_fields);
    
    int *format_indices = get_format_indices_per_file(links, 4, files, 4, format_fields);
    
    int gt_pos = get_field_position_in_format("GT", strdup(format));
    int filter_pos = get_field_position_in_format("SFT", strdup(format));
    int info_pos = get_field_position_in_format("IN", strdup(format));
    options->missing_mode = MISSING_CONDITION;
    char *empty_sample = get_empty_sample(format_fields->size, gt_pos, options);
    
    array_list_t *samples = merge_samples(links, 4, files, 4, alleles_table, format_fields, format_indices, empty_sample, gt_pos, filter_pos, info_pos, options);
    printf("Samples:\n");
    array_list_print(samples);
    printf("\n------------------------\n");
    
    fail_if(strcmp(array_list_get(0, samples), "1/1:20:40:30:."), "Sample (0,0) must be 1/1:20:40:30:.");
    fail_if(strcmp(array_list_get(1, samples), "0/1:10:60:50:."), "Sample (0,1) must be 0/1:10:60:50:.");
    fail_if(strcmp(array_list_get(2, samples), "0/0:30:50:70:."), "Sample (0,2) must be 0/0:30:50:70:.");
    
    fail_if(strcmp(array_list_get(3, samples), "2/2:.:.:.:40"), "Sample (1,0) must be 2/2:.:.:.:40");
    fail_if(strcmp(array_list_get(4, samples), "0/2:.:.:.:60"), "Sample (1,1) must be 0/2:.:.:.:60");
    fail_if(strcmp(array_list_get(5, samples), "0/0:.:.:.:50"), "Sample (1,2) must be 0/0:.:.:.:50");
    
    fail_if(strcmp(array_list_get(6, samples), "3/3:30:.:40:20"), "Sample (2,0) must be 3/3:30:.:40:20");
    
    fail_if(strcmp(array_list_get(7, samples), "1/1:.:.:.:."), "Sample (3,0) must be 1/1:.:.:.:.");
    fail_if(strcmp(array_list_get(8, samples), "0/1:.:.:.:."), "Sample (3,1) must be 0/1:.:.:.:.");
    
    
    // Merge samples of position only in files (0,2)
    links[0]->file = files[0];
    links[0]->record = input[0];
    links[1]->file = files[2];
    links[1]->record = input[2];
    
    alleles_table = cp_hashtable_create(8, cp_hash_istring, (cp_compare_fn) strcasecmp);
    merge_alternate_field(links, 2, alleles_table);
    
    format_fields = array_list_new(8, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    format_fields->compare_fn = strcasecmp;
    format = merge_format_field(links, 2, options, format_fields);
    
    format_indices = get_format_indices_per_file(links, 2, files, 4, format_fields);
    
    gt_pos = get_field_position_in_format("GT", strdup(format));
    filter_pos = get_field_position_in_format("SFT", strdup(format));
    info_pos = get_field_position_in_format("IN", strdup(format));
    options->missing_mode = REFERENCE;
    empty_sample = get_empty_sample(format_fields->size, gt_pos, options);
    
    samples = merge_samples(links, 2, files, 4, alleles_table, format_fields, format_indices, empty_sample, gt_pos, filter_pos, info_pos, options);
    printf("Samples in (0,2):\n");
    array_list_print(samples);
    printf("\n------------------------\n");
    
    fail_if(strcmp(array_list_get(0, samples), "1/1:20:40:30:."), "Sample (0,0) must be 1/1:20:40:30:.");
    fail_if(strcmp(array_list_get(1, samples), "0/1:10:60:50:."), "Sample (0,1) must be 0/1:10:60:50:.");
    fail_if(strcmp(array_list_get(2, samples), "0/0:30:50:70:."), "Sample (0,2) must be 0/0:30:50:70:.");
    
    fail_if(strcmp(array_list_get(3, samples), "0/0:.:.:.:."), "Sample (1,0) must be 0/0:.:.:.:.");
    fail_if(strcmp(array_list_get(4, samples), "0/0:.:.:.:."), "Sample (1,1) must be 0/0:.:.:.:.");
    fail_if(strcmp(array_list_get(5, samples), "0/0:.:.:.:."), "Sample (1,2) must be 0/0:.:.:.:.");
    
    fail_if(strcmp(array_list_get(6, samples), "2/2:30:.:40:20"), "Sample (2,0) must be 2/2:30:.:40:20");
    
    fail_if(strcmp(array_list_get(7, samples), "0/0:.:.:.:."), "Sample (3,0) must be 0/0:.:.:.:.");
    fail_if(strcmp(array_list_get(8, samples), "0/0:.:.:.:."), "Sample (3,1) must be 0/0:.:.:.:.");
    
    
    // Merge samples of position only in files (2,0)
    links[0]->file = files[2];
    links[0]->record = input[2];
    links[1]->file = files[0];
    links[1]->record = input[0];
    
    alleles_table = cp_hashtable_create(8, cp_hash_istring, (cp_compare_fn) strcasecmp);
    merge_alternate_field(links, 2, alleles_table);
    
    format_fields = array_list_new(8, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    format_fields->compare_fn = strcasecmp;
    format = merge_format_field(links, 2, options, format_fields);
    
    format_indices = get_format_indices_per_file(links, 2, files, 4, format_fields);
    
    gt_pos = get_field_position_in_format("GT", strdup(format));
    filter_pos = get_field_position_in_format("SFT", strdup(format));
    info_pos = get_field_position_in_format("IN", strdup(format));
    options->missing_mode = MISSING_CONDITION;
    empty_sample = get_empty_sample(format_fields->size, gt_pos, options);
    
    samples = merge_samples(links, 2, files, 4, alleles_table, format_fields, format_indices, empty_sample, gt_pos, filter_pos, info_pos, options);
    printf("Samples in (2,0):\n");
    array_list_print(samples);
    printf("\n------------------------\n");
    
    fail_if(strcmp(array_list_get(0, samples), ".:30:2/2:20:40"), "Sample (0,0) must be .:30:2/2:20:40");
    fail_if(strcmp(array_list_get(1, samples), ".:50:0/2:10:60"), "Sample (0,1) must be .:50:0/2:10:60");
    fail_if(strcmp(array_list_get(2, samples), ".:70:0/0:30:50"), "Sample (0,2) must be .:70:0/0:30:50");
    
    fail_if(strcmp(array_list_get(3, samples), ".:.:./.:.:."), "Sample (1,0) must be .:.:./.:.:.");
    fail_if(strcmp(array_list_get(4, samples), ".:.:./.:.:."), "Sample (1,1) must be .:.:./.:.:.");
    fail_if(strcmp(array_list_get(5, samples), ".:.:./.:.:."), "Sample (1,2) must be .:.:./.:.:.");
    
    fail_if(strcmp(array_list_get(6, samples), "20:40:1/1:30:."), "Sample (2,0) must be 20:40:1/1:30:.");
    
    fail_if(strcmp(array_list_get(7, samples), ".:.:./.:.:."), "Sample (3,0) must be .:.:./.:.:.");
    fail_if(strcmp(array_list_get(8, samples), ".:.:./.:.:."), "Sample (3,1) must be .:.:./.:.:.");
}
END_TEST

START_TEST (merge_info_test) {
    vcf_record_t *input[4];
    input[0] = create_example_record_0();
    input[1] = create_example_record_1();
    input[2] = create_example_record_2();
    input[3] = create_example_record_3();
    
    vcf_record_t *result = vcf_record_new();
    
    vcf_record_file_link **links = calloc (4, sizeof(vcf_record_file_link*));
    for (int i = 0; i < 4; i++) {
        links[i] = malloc(sizeof(vcf_record_file_link));
        links[i]->file = files[i];
        links[i]->record = input[i];
    }
    
    set_vcf_record_id(input[0]->id, input[0]->id_len, result);
    set_vcf_record_filter(input[0]->filter, input[0]->filter_len, result);
    set_vcf_record_quality(merge_quality_field(links, 4), result);
    
    cp_hashtable *alleles_table = cp_hashtable_create(8, cp_hash_istring, (cp_compare_fn) strcasecmp);
    char *alternate = merge_alternate_field(links, 4, alleles_table);
    set_vcf_record_alternate(alternate, strlen(alternate), result);
    
    array_list_t *format_fields = array_list_new(8, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    format_fields->compare_fn = strcasecmp;
    char *format = merge_format_field(links, 4, options, format_fields);
    set_vcf_record_format(format, strlen(format), result);
    
    int *format_indices = get_format_indices_per_file(links, 4, files, 4, format_fields);
    
    int gt_pos = get_field_position_in_format("GT", strndup(result->format, result->format_len));
    int filter_pos = get_field_position_in_format("SFT", strndup(result->format, result->format_len));
    int info_pos = get_field_position_in_format("IN", strndup(result->format, result->format_len));
    options->missing_mode = MISSING_CONDITION;
    char *empty_sample = get_empty_sample(format_fields->size, gt_pos, options);
    
    result->samples = merge_samples(links, 4, files, 4, alleles_table, format_fields, format_indices, empty_sample, gt_pos, filter_pos, info_pos, options);
    
    int num_info_fields = 13;
    options->info_fields = malloc (num_info_fields * sizeof(char*));
    options->info_fields[0] = "AC";
    options->info_fields[1] = "AF";
    options->info_fields[2] = "AN";
    options->info_fields[3] = "DB";
    options->info_fields[4] = "DP";
    options->info_fields[5] = "H2";
    options->info_fields[6] = "H3";
    options->info_fields[7] = "MQ";
    options->info_fields[8] = "MQ0";
    options->info_fields[9] = "QD";
    options->info_fields[10] = "SOMATIC";
    options->info_fields[11] = "VALIDATED";
    options->info_fields[12] = "NS";
    
    char *info = merge_info_field(links, 4, options->info_fields, num_info_fields, result, alleles_table, empty_sample);
    
    printf("info = %s\n", info);
    
    fail_if(!strstr(info, "DB"), "INFO/DB must be present");
    fail_if(!strstr(info, "H2"), "INFO/H2 must be present");
    fail_if(strstr(info, "H3"), "INFO/H3 must not be present");
    fail_if(strstr(info, "SOMATIC"), "INFO/SOMATIC must not be present");
    fail_if(strstr(info, "VALIDATED"), "INFO/VALIDATED must not be present");
    
    char **split_info = split(info, ";", &num_info_fields);
    
    fail_if(strcmp(split_info[0], "AC=6,3,2"), "INFO/AC value must be 6,3,2");
    fail_if(strcmp(split_info[1], "AF=0.545,0.273,0.182"), "INFO/AC value must be 0.545,0.273,0.182");
    fail_if(strcmp(split_info[2], "AN=4"), "INFO/AN value must be 4");
    fail_if(strcmp(split_info[3], "DB"), "INFO/DB must have no value");
    fail_if(strcmp(split_info[4], "DP=150"), "INFO/DP value must be 150");
    fail_if(strcmp(split_info[5], "H2"), "INFO/H2 must have no value");
    fail_if(strcmp(split_info[6], "MQ=15.986"), "INFO/MQ value must be 15.986"); // sqrt( (10^2+20^2+30^2+30^2) / 9 )
    fail_if(strcmp(split_info[7], "MQ0=5"), "INFO/MQ0 value must be 5");
    fail_if(strcmp(split_info[8], "QD=0.119"), "INFO/QD value must be 0.119"); // QUAL = (20*3+30*3+10)/9, DP = 150, QD = 0.119
    fail_if(strcmp(split_info[9], "NS=9"), "INFO/NS value must be 9");
    
}
END_TEST

START_TEST (add_filter_test) {
    vcf_record_t *input[4];
    input[0] = create_example_record_0();
    input[1] = create_example_record_1();
    input[2] = create_example_record_2();
    input[3] = create_example_record_3();
    
    vcf_record_t *result = vcf_record_new();
    
    vcf_record_file_link **links = calloc (4, sizeof(vcf_record_file_link*));
    for (int i = 0; i < 4; i++) {
        links[i] = malloc(sizeof(vcf_record_file_link));
        links[i]->file = files[i];
        links[i]->record = input[i];
    }
    
    set_vcf_record_id(input[0]->id, input[0]->id_len, result);
    set_vcf_record_filter(input[0]->filter, input[0]->filter_len, result);
    set_vcf_record_quality(merge_quality_field(links, 4), result);
    
    cp_hashtable *alleles_table = cp_hashtable_create(8, cp_hash_istring, (cp_compare_fn) strcasecmp);
    char *alternate = merge_alternate_field(links, 4, alleles_table);
    set_vcf_record_alternate(alternate, strlen(alternate), result);
    
    array_list_t *format_fields = array_list_new(8, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    format_fields->compare_fn = strcasecmp;
    
    // Insert FILTER
    options->copy_filter = 1;
    options->copy_info = 0;
    char *format = merge_format_field(links, 4, options, format_fields);
    fail_if(strcmp(format, "GT:GQ:DP:HQ:RD:SFT"), "Format with filter must be GT:GQ:DP:HQ:RD:SFT");
    set_vcf_record_format(format, strlen(format), result);
    
    int *format_indices = get_format_indices_per_file(links, 4, files, 4, format_fields);
    
    int gt_pos = get_field_position_in_format("GT", strndup(result->format, result->format_len));
    int filter_pos = get_field_position_in_format("SFT", strndup(result->format, result->format_len));
    int info_pos = get_field_position_in_format("IN", strndup(result->format, result->format_len));
    options->missing_mode = MISSING_CONDITION;
    char *empty_sample = get_empty_sample(format_fields->size, gt_pos, options);
    
    array_list_t* samples = merge_samples(links, 4, files, 4, alleles_table, format_fields, format_indices, empty_sample, gt_pos, filter_pos, info_pos, options);
    
    fail_if(strcmp(array_list_get(0, samples), "1/1:20:40:30:.:PASS"), "Sample (0,0) must be 1/1:20:40:30:.:PASS");
    fail_if(strcmp(array_list_get(1, samples), "0/1:10:60:50:.:PASS"), "Sample (0,1) must be 0/1:10:60:50:.:PASS");
    fail_if(strcmp(array_list_get(2, samples), "0/0:30:50:70:.:PASS"), "Sample (0,2) must be 0/0:30:50:70:.:PASS");
    
    fail_if(strcmp(array_list_get(3, samples), "2/2:.:.:.:40:STD_FILTER"), "Sample (1,0) must be 2/2:.:.:.:40:STD_FILTER");
    fail_if(strcmp(array_list_get(4, samples), "0/2:.:.:.:60:STD_FILTER"), "Sample (1,1) must be 0/2:.:.:.:60:STD_FILTER");
    fail_if(strcmp(array_list_get(5, samples), "0/0:.:.:.:50:STD_FILTER"), "Sample (1,2) must be 0/0:.:.:.:50:STD_FILTER");
    
    fail_if(strcmp(array_list_get(6, samples), "3/3:30:.:40:20:q10"), "Sample (2,0) must be 3/3:30:.:40:20:q10");
    
    fail_if(strcmp(array_list_get(7, samples), "1/1:.:.:.:.:."), "Sample (3,0) must be 1/1:.:.:.:.:.");
    fail_if(strcmp(array_list_get(8, samples), "0/1:.:.:.:.:."), "Sample (3,1) must be 0/1:.:.:.:.:.");
    
}
END_TEST

START_TEST (add_info_test) {
    vcf_record_t *input[4];
    input[0] = create_example_record_0();
    input[1] = create_example_record_1();
    input[2] = create_example_record_2();
    input[3] = create_example_record_3();
    
    vcf_record_t *result = vcf_record_new();
    
    vcf_record_file_link **links = calloc (4, sizeof(vcf_record_file_link*));
    for (int i = 0; i < 4; i++) {
        links[i] = malloc(sizeof(vcf_record_file_link));
        links[i]->file = files[i];
        links[i]->record = input[i];
    }
    
    set_vcf_record_id(input[0]->id, input[0]->id_len, result);
    set_vcf_record_filter(input[0]->filter, input[0]->filter_len, result);
    set_vcf_record_quality(merge_quality_field(links, 4), result);
    
    cp_hashtable *alleles_table = cp_hashtable_create(8, cp_hash_istring, (cp_compare_fn) strcasecmp);
    char *alternate = merge_alternate_field(links, 4, alleles_table);
    set_vcf_record_alternate(alternate, strlen(alternate), result);
    
    array_list_t *format_fields = array_list_new(8, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    format_fields->compare_fn = strcasecmp;
        
    // Insert INFO
    options->copy_filter = 0;
    options->copy_info = 1;
    char *format = merge_format_field(links, 4, options, format_fields);
    fail_if(strcmp(format, "GT:GQ:DP:HQ:RD:IN"), "Format with info must be GT:GQ:DP:HQ:RD:IN");
    set_vcf_record_format(format, strlen(format), result);
    
    int *format_indices = get_format_indices_per_file(links, 4, files, 4, format_fields);
    
    int gt_pos = get_field_position_in_format("GT", strndup(result->format, result->format_len));
    int filter_pos = get_field_position_in_format("SFT", strndup(result->format, result->format_len));
    int info_pos = get_field_position_in_format("IN", strndup(result->format, result->format_len));
    options->missing_mode = MISSING_CONDITION;
    char *empty_sample = get_empty_sample(format_fields->size, gt_pos, options);
    
    array_list_t* samples = merge_samples(links, 4, files, 4, alleles_table, format_fields, format_indices, empty_sample, gt_pos, filter_pos, info_pos, options);
    
    fail_if(strcmp(array_list_get(0, samples), "1/1:20:40:30:.:NS=3;DP=14;H2"), "Sample (0,0) must be 1/1:20:40:30:.:NS=3;DP=14;H2");
    fail_if(strcmp(array_list_get(1, samples), "0/1:10:60:50:.:NS=3;DP=14;H2"), "Sample (0,1) must be 0/1:10:60:50:.:NS=3;DP=14;H2");
    fail_if(strcmp(array_list_get(2, samples), "0/0:30:50:70:.:NS=3;DP=14;H2"), "Sample (0,2) must be 0/0:30:50:70:.:NS=3;DP=14;H2");
    
    fail_if(strcmp(array_list_get(3, samples), "2/2:.:.:.:40:DP=10;NS=4;AF=0.5;H2"), "Sample (1,0) must be 2/2:.:.:.:40:DP=10;NS=4;AF=0.5;H2");
    fail_if(strcmp(array_list_get(4, samples), "0/2:.:.:.:60:DP=10;NS=4;AF=0.5;H2"), "Sample (1,1) must be 0/2:.:.:.:60:DP=10;NS=4;AF=0.5;H2");
    fail_if(strcmp(array_list_get(5, samples), "0/0:.:.:.:50:DP=10;NS=4;AF=0.5;H2"), "Sample (1,2) must be 0/0:.:.:.:50:DP=10;NS=4;AF=0.5;H2");
    
    fail_if(strcmp(array_list_get(6, samples), "3/3:30:.:40:20:AF=0.5;NS=3;DP=14;DB"), "Sample (2,0) must be 3/3:30:.:40:20:AF=0.5;NS=3;DP=14;DB");
    
    fail_if(strcmp(array_list_get(7, samples), "1/1:.:.:.:.:DB;H2"), "Sample (3,0) must be 1/1:.:.:.:.:DB;H2");
    fail_if(strcmp(array_list_get(8, samples), "0/1:.:.:.:.:DB;H2"), "Sample (3,1) must be 0/1:.:.:.:.:DB;H2");
}
END_TEST

START_TEST (add_info_filter_test) {
    vcf_record_t *input[4];
    input[0] = create_example_record_0();
    input[1] = create_example_record_1();
    input[2] = create_example_record_2();
    input[3] = create_example_record_3();
    
    vcf_record_t *result = vcf_record_new();
    
    vcf_record_file_link **links = calloc (4, sizeof(vcf_record_file_link*));
    for (int i = 0; i < 4; i++) {
        links[i] = malloc(sizeof(vcf_record_file_link));
        links[i]->file = files[i];
        links[i]->record = input[i];
    }
    
    set_vcf_record_id(input[0]->id, input[0]->id_len, result);
    set_vcf_record_filter(input[0]->filter, input[0]->filter_len, result);
    set_vcf_record_quality(merge_quality_field(links, 4), result);
    
    cp_hashtable *alleles_table = cp_hashtable_create(8, cp_hash_istring, (cp_compare_fn) strcasecmp);
    char *alternate = merge_alternate_field(links, 4, alleles_table);
    set_vcf_record_alternate(alternate, strlen(alternate), result);
    
    array_list_t *format_fields = array_list_new(8, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    format_fields->compare_fn = strcasecmp;
        
    // Insert FILTER and INFO
    options->copy_filter = 1;
    options->copy_info = 1;
    char *format = merge_format_field(links, 4, options, format_fields);
    fail_if(strcmp(format, "GT:GQ:DP:HQ:RD:SFT:IN"), "Format with filter must be GT:GQ:DP:HQ:RD:SFT:IN");
    set_vcf_record_format(format, strlen(format), result);
    
    int *format_indices = get_format_indices_per_file(links, 4, files, 4, format_fields);
    
    char *format_bak = strndup(result->format, result->format_len);
    int gt_pos = get_field_position_in_format("GT", strndup(result->format, result->format_len));
    int filter_pos = get_field_position_in_format("SFT", strndup(result->format, result->format_len));
    int info_pos = get_field_position_in_format("IN", strndup(result->format, result->format_len));
    options->missing_mode = MISSING_CONDITION;
    char *empty_sample = get_empty_sample(format_fields->size, gt_pos, options);
    
    array_list_t* samples = merge_samples(links, 4, files, 4, alleles_table, format_fields, format_indices, empty_sample, gt_pos, filter_pos, info_pos, options);
    
    fail_if(strcmp(array_list_get(0, samples), "1/1:20:40:30:.:PASS:NS=3;DP=14;H2"), "Sample (0,0) must be 1/1:20:40:30:.:PASS:NS=3;DP=14;H2");
    fail_if(strcmp(array_list_get(1, samples), "0/1:10:60:50:.:PASS:NS=3;DP=14;H2"), "Sample (0,1) must be 0/1:10:60:50:.:PASS:NS=3;DP=14;H2");
    fail_if(strcmp(array_list_get(2, samples), "0/0:30:50:70:.:PASS:NS=3;DP=14;H2"), "Sample (0,2) must be 0/0:30:50:70:.:PASS:NS=3;DP=14;H2");
    
    fail_if(strcmp(array_list_get(3, samples), "2/2:.:.:.:40:STD_FILTER:DP=10;NS=4;AF=0.5;H2"), "Sample (1,0) must be 2/2:.:.:.:40:STD_FILTER:DP=10;NS=4;AF=0.5;H2");
    fail_if(strcmp(array_list_get(4, samples), "0/2:.:.:.:60:STD_FILTER:DP=10;NS=4;AF=0.5;H2"), "Sample (1,1) must be 0/2:.:.:.:60:STD_FILTER:DP=10;NS=4;AF=0.5;H2");
    fail_if(strcmp(array_list_get(5, samples), "0/0:.:.:.:50:STD_FILTER:DP=10;NS=4;AF=0.5;H2"), "Sample (1,2) must be 0/0:.:.:.:50:STD_FILTER:DP=10;NS=4;AF=0.5;H2");
    
    fail_if(strcmp(array_list_get(6, samples), "3/3:30:.:40:20:q10:AF=0.5;NS=3;DP=14;DB"), "Sample (2,0) must be 3/3:30:.:40:20:q10:AF=0.5;NS=3;DP=14;DB");
    
    fail_if(strcmp(array_list_get(7, samples), "1/1:.:.:.:.:.:DB;H2"), "Sample (3,0) must be 1/1:.:.:.:.:.:DB;H2");
    fail_if(strcmp(array_list_get(8, samples), "0/1:.:.:.:.:.:DB;H2"), "Sample (3,1) must be 0/1:.:.:.:.:.:DB;H2");
}
END_TEST


START_TEST (get_format_indices_per_file_test) {
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
    char *format = merge_format_field(links, 4, options, format_fields);
    
    // Get from position in all files
    int file_idx;
    int *indices = get_format_indices_per_file(links, 4, files, 4, format_fields);
    printf("The coordinates of incorrect format indices will be presented as (file,index in file)\n");
    
//     printf("------------\n");
//     for (int i = 0; i < 20; i++) {
//         printf("%d ", indices[i]);
//     }
//     printf("------------\n");
    
    file_idx = 0;
    fail_if(indices[file_idx*5 + 0] != 0, "(0,0) must be 0");
    fail_if(indices[file_idx*5 + 1] != 1, "(0,1) must be 1");
    fail_if(indices[file_idx*5 + 2] != 2, "(0,2) must be 2");
    fail_if(indices[file_idx*5 + 3] != 3, "(0,3) must be 3");
    fail_if(indices[file_idx*5 + 4] != -1, "(0,4) must be -1");
    
    file_idx = 1;
    fail_if(indices[file_idx*5 + 0] != 0, "(1,0) must be 0");
    fail_if(indices[file_idx*5 + 1] != -1, "(1,1) must be -1");
    fail_if(indices[file_idx*5 + 2] != -1, "(1,2) must be -1");
    fail_if(indices[file_idx*5 + 3] != -1, "(1,3) must be -1");
    fail_if(indices[file_idx*5 + 4] != 1, "(1,4) must be 1");
    
    file_idx = 2;
    fail_if(indices[file_idx*5 + 0] != 2, "(2,0) must be 2");
    fail_if(indices[file_idx*5 + 1] != 3, "(2,1) must be 3");
    fail_if(indices[file_idx*5 + 2] != -1, "(2,2) must be -1");
    fail_if(indices[file_idx*5 + 3] != 1, "(2,3) must be 1");
    fail_if(indices[file_idx*5 + 4] != 0, "(2,4) must be 0");
    
    file_idx = 3;
    fail_if(indices[file_idx*5 + 0] != 0, "(3,0) must be 0");
    fail_if(indices[file_idx*5 + 1] != -1, "(3,1) must be -1");
    fail_if(indices[file_idx*5 + 2] != -1, "(3,2) must be -1");
    fail_if(indices[file_idx*5 + 3] != -1, "(3,3) must be -1");
    fail_if(indices[file_idx*5 + 4] != -1, "(3,4) must be -1");
    
    free(indices);
    
    // Get from position only in files 0,2
    links[0]->file = files[0];
    links[0]->record = input[0];
    links[1]->file = files[2];
    links[1]->record = input[2];
    indices = get_format_indices_per_file(links, 2, files, 4, format_fields);
    
//     printf("------------\n");
//     for (int i = 0; i < 20; i++) {
//         printf("%d ", indices[i]);
//     }
//     printf("------------\n");
    
    file_idx = 0;
    fail_if(indices[file_idx*5 + 0] != 0, "(0,0) must be 0");
    fail_if(indices[file_idx*5 + 1] != 1, "(0,1) must be 1");
    fail_if(indices[file_idx*5 + 2] != 2, "(0,2) must be 2");
    fail_if(indices[file_idx*5 + 3] != 3, "(0,3) must be 3");
    fail_if(indices[file_idx*5 + 4] != -1, "(0,4) must be -1");
    
    file_idx = 1;
    fail_if(indices[file_idx*5 + 0] != -1, "(1,0) must be -1");
    fail_if(indices[file_idx*5 + 1] != -1, "(1,1) must be -1");
    fail_if(indices[file_idx*5 + 2] != -1, "(1,2) must be -1");
    fail_if(indices[file_idx*5 + 3] != -1, "(1,3) must be -1");
    fail_if(indices[file_idx*5 + 4] != -1, "(1,4) must be -1");
    
    file_idx = 2;
    fail_if(indices[file_idx*5 + 0] != 2, "(2,0) must be 2");
    fail_if(indices[file_idx*5 + 1] != 3, "(2,1) must be 3");
    fail_if(indices[file_idx*5 + 2] != -1, "(2,2) must be -1");
    fail_if(indices[file_idx*5 + 3] != 1, "(2,3) must be 1");
    fail_if(indices[file_idx*5 + 4] != 0, "(2,4) must be 0");
    
    file_idx = 3;
    fail_if(indices[file_idx*5 + 0] != -1, "(3,0) must be -1");
    fail_if(indices[file_idx*5 + 1] != -1, "(3,1) must be -1");
    fail_if(indices[file_idx*5 + 2] != -1, "(3,2) must be -1");
    fail_if(indices[file_idx*5 + 3] != -1, "(3,3) must be -1");
    fail_if(indices[file_idx*5 + 4] != -1, "(3,4) must be -1");
    
    free(indices);
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
    TCase *tc_headers = tcase_create("Merge headers from several files");
    tcase_add_checked_fixture(tc_headers, setup_merge_headers, teardown_merge_headers);
    tcase_add_test(tc_headers, merge_headers_test);
    
    TCase *tc_unique = tcase_create("Merge position in one file");
    tcase_add_checked_fixture(tc_unique, setup_merge_positions, teardown_merge_positions);
    tcase_add_test(tc_unique, merge_position_in_one_file);
    
    TCase *tc_auxiliary = tcase_create("Auxiliary functions");
    tcase_add_checked_fixture(tc_auxiliary, setup_merge_positions, teardown_merge_positions);
    tcase_add_test(tc_auxiliary, get_format_indices_per_file_test);
    
    TCase *tc_repeated = tcase_create("Merge position in several files");
    tcase_add_checked_fixture(tc_repeated, setup_merge_positions, teardown_merge_positions);
    tcase_add_test(tc_repeated, merge_id_test);
    tcase_add_test(tc_repeated, merge_alternate_test);
    tcase_add_test(tc_repeated, merge_quality_test);
    tcase_add_test(tc_repeated, merge_filter_test);
    tcase_add_test(tc_repeated, merge_format_test);
    tcase_add_test(tc_repeated, merge_samples_test);
    tcase_add_test(tc_repeated, merge_info_test);
    
    TCase *tc_extra = tcase_create("Adding extra fields to samples");
    tcase_add_checked_fixture(tc_extra, setup_merge_positions, teardown_merge_positions);
    tcase_add_test(tc_extra, add_filter_test);
    tcase_add_test(tc_extra, add_info_test);
    tcase_add_test(tc_extra, add_info_filter_test);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Check for hpg-vcf/merge");
    suite_add_tcase(fs, tc_headers);
    suite_add_tcase(fs, tc_unique);
    suite_add_tcase(fs, tc_auxiliary);
    suite_add_tcase(fs, tc_repeated);
    suite_add_tcase(fs, tc_extra);
    
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
    input->info = "AF=0.5;NS=3;DP=14;DB";
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
