#include "merge.h"

/* ******************************
 *      Splitting per variant   *
 * ******************************/
 
// merge_result_t *new_merge_result(vcf_record_t *record, char *merge_name) {
//     merge_result_t *result = (merge_result_t*) malloc (sizeof(merge_result_t));
//     result->record = record;
//     result->merge_name = merge_name;
//     return result;
// }
// 
// void free_merge_result(merge_result_t* merge_result) {
//     free(merge_result->merge_name);
//     free(merge_result);
// }

int merge(vcf_record_t **variants, int num_variants, list_t* output_list) {
    // TODO
    printf("Merging %d variants...\n", num_variants);
    
    return 0;
}
