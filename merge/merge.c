#include "merge.h"


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

// int merge(vcf_record_t **variants, int num_variants, list_t* output_list) {
int merge(array_list_t **records_by_position, int num_positions, vcf_file_t **files, int num_files, merge_options_data_t *options, list_t *output_list) {
    // TODO get samples list, sorted by file (as provided as input)
    array_list_t *samples = get_global_samples(files, num_files);
    array_list_print(samples);
    
    // TODO
    vcf_record_t *merged;
    for (int i = 0; i < num_positions; i++) {
        if (records_by_position[i]->size == 1) {
            // TODO merge position present just in one file
            merged = merge_unique_position(records_by_position[i]->items[0], files, num_files, options);
        } else if (records_by_position[i]->size > 1) {
            // TODO merge position present in more than one file
        }
    }
    
    printf("Merging %d variants...\n", num_positions);
    
    array_list_free(samples, free);
    
    return 0;
}

vcf_record_t *merge_unique_position(vcf_record_file_link *position, vcf_file_t **files, int num_files, merge_options_data_t *options) {
    vcf_record_t *input = position->record;
    vcf_record_t *result = create_record();
    result->chromosome = strdup(input->chromosome);
    result->position = input->position;
    result->id = strdup(input->id);
    result->reference = strdup(input->reference);
    result->alternate = strdup(input->alternate);
    
    result->format = strdup(input->format);
    int num_format_fields = 1;
    for (int i = 0; i < strlen(result->format); i++) {
        if (result->format[i] == ':') {
            num_format_fields++;
        }
    }
    char *format_bak = strdup(result->format);
    int gt_pos = get_field_position_in_format("gt", format_bak);
    
    for (int i = 0; i < num_files; i++) {
        if(!strcmp(files[i]->filename, position->file->filename)) {
            // Samples of the file where the position has been read are directly copied
            array_list_insert_all(input->samples->items, files[i]->num_samples, result->samples);
        } else {
            // Samples in the rest of files must be filled according to the specified format
            int sample_len = num_format_fields * 2 + 1; // Each field + ':' = 2 chars, except for GT which is 1 char more
            char *sample[sample_len];
            memset(sample, 0, sample_len * sizeof(char));
            for (int j = 0; j < num_format_fields; j++) {
                if (j > 0) {
                    strncat(sample, ":", 1);
                }
                if (j != gt_pos) {
                    strncat(sample, ".", 1);
                } else {
                    if (options->missing_mode == MISSING) {
                        strncat(sample, "./.", 3);
                    } else if (options->missing_mode == REFERENCE) {
                        strncat(sample, "0/0", 3);
                    }
                }
                array_list_insert(sample, result->samples);
            }
        }
    }
    
    // INFO field must be calculated based on statistics about the new list of samples
    result->info = recalculate_info_field(result->samples->items, result->samples->size);
    
    return result;
}


array_list_t *get_global_samples(vcf_file_t **files, int num_files) {
    array_list_t *samples = array_list_new(files[0]->num_samples * 2, 2, COLLECTION_MODE_ASYNCHRONIZED);
    for (int i = 0; i < num_files; i++) {
        array_list_insert_all(files[i]->samples_names->items, files[i]->samples_names->size, samples);
    }
    return samples;
}

char *recalculate_info_field(char **samples, size_t num_samples) {
    // TODO
    return strdup(samples[0]);
}
