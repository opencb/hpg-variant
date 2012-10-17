#include "split.h"

/* ******************************
 *      Splitting per variant   *
 * ******************************/
 
split_result_t *new_split_result(vcf_record_t *record, char *split_name) {
    split_result_t *result = (split_result_t*) malloc (sizeof(split_result_t));
    result->record = record;
    result->split_name = split_name;
    return result;
}

void free_split_result(split_result_t* split_result) {
    free(split_result->split_name);
    free(split_result);
}

int split_by_chromosome(vcf_record_t **variants, int num_variants, list_t* output_list) {
    char *output_prefix;
    vcf_record_t *record;
    split_result_t *split_result;
    
    // For each variant, its output filename will be 'chromosome<#chr>_<original_filename>.vcf'
    for (int i = 0; i < num_variants; i++) {
        record = variants[i];
        output_prefix = (char*) calloc (record->chromosome_len + 12, sizeof(char));
        strncat(output_prefix, "chromosome_", 11);
        strncat(output_prefix, record->chromosome, record->chromosome_len);
        split_result = new_split_result(vcf_record_copy(record), output_prefix);
        
        // Insert results in output list
        list_item_t *item = list_item_new(i, 0, split_result);
        list_insert_item(item, output_list);
    }
    
    return 0;
}
