#include "split.h"

/* ******************************
 *      Splitting per variant   *
 * ******************************/
 
variant_split_result_t *new_variant_split_result(vcf_record_t *record, char *split_name) {
    variant_split_result_t *result = (variant_split_result_t*) malloc (sizeof(variant_split_result_t));
    result->record = record;
    result->split_name = split_name;
    return result;
}

void free_variant_split_result(variant_split_result_t* split_result) {
    free(split_result->split_name);
    free(split_result);
}


int split_by_chromosome(list_item_t* variants, int num_variants, list_t* output_list) {
    char *output_filename;
    vcf_record_t *record;
    variant_split_result_t *split_result;
    list_item_t *cur_variant = variants;
    
    // For each variant, its output filename will be 'chromosome<#chr>_<original_filename>.vcf'
    for (int i = 0; i < num_variants && cur_variant != NULL; i++, cur_variant = cur_variant->next_p) {
        record = (vcf_record_t*) cur_variant->data_p;
        output_filename = (char*) calloc (strlen(record->chromosome) + 12, sizeof(char));
        strncat(output_filename, "chromosome_", 11);
        strncat(output_filename, record->chromosome, strlen(record->chromosome);
        split_result = new_variant_split_result(record->chromosome, output_filename);
        
        // Insert results in output list
        list_item_t *item = list_item_new(i, 0, split_result);
        list_insert_item(item, output_list);
    }
    
    return 0;
}
