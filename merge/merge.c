#include "merge.h"
#include <assert.h>

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
    result->quality = input->quality;
    result->filter = strdup(input->filter);
    
    result->format = strdup(input->format);
    int num_format_fields = 1;
    for (int i = 0; i < strlen(result->format); i++) {
        if (result->format[i] == ':') {
            num_format_fields++;
        }
    }
    char *format_bak = strdup(result->format);
    int gt_pos = get_field_position_in_format("GT", format_bak);
    
    // Create the text for empty samples
    char *empty_sample = get_empty_sample(num_format_fields, gt_pos, options->missing_mode);
    
    // Fill list of samples
    for (int i = 0; i < num_files; i++) {
        if(!strcmp(files[i]->filename, position->file->filename)) {
            // Samples of the file where the position has been read are directly copied
            array_list_insert_all(input->samples->items, files[i]->num_samples, result->samples);
            LOG_DEBUG_F("%d samples of file %s inserted\n", files[i]->num_samples, files[i]->filename);
        } else {
            // Samples in the rest of files must be filled according to the specified format
            for (int j = 0; j < files[i]->num_samples; j++) {
                array_list_insert(strdup(empty_sample), result->samples);
            }
            LOG_DEBUG_F("file %s\t%d empty samples inserted\n", files[i]->filename, files[i]->num_samples);
        }
    }
    
    // INFO field must be calculated based on statistics about the new list of samples
    // TODO add info-fields as argument (AF, NS and so on)
    result->info = generate_info_field(result->samples->items, result->samples->size);
    
    return result;
}

vcf_record_t *merge_shared_position(vcf_record_file_link **position_in_files, int position_occurrences, 
                                    vcf_file_t **files, int num_files, merge_options_data_t *options, int *info) {
    if (position_occurrences < 1) {
        return NULL;
    }
    
    vcf_file_t *file;
    vcf_record_t *input;
    vcf_record_t *result = create_record();
    
    // Check consistency among chromosome, position and reference
    for (int i = 0; i < position_occurrences-1; i++) {
        if (strcmp(position_in_files[i]->record->chromosome, position_in_files[i+1]->record->chromosome)) {
            LOG_ERROR("Positions can't be merged: Discordant chromosome\n");
            *info = DISCORDANT_CHROMOSOME;
            return NULL;
        } else if (position_in_files[i]->record->position != position_in_files[i+1]->record->position) {
            LOG_ERROR("Positions can't be merged: Discordant position\n");
            *info = DISCORDANT_POSITION;
            return NULL;
        } else if (strcmp(position_in_files[i]->record->reference, position_in_files[i+1]->record->reference)) {
            LOG_ERROR("Positions can't be merged: Discordant reference allele\n");
            *info = DISCORDANT_REFERENCE;
            return NULL;
        }
    }
    
    // Once checked, copy fields that are constant in all files
    result->chromosome = strdup(position_in_files[0]->record->chromosome);
    result->position = position_in_files[0]->record->position;
    result->reference = strdup(position_in_files[0]->record->reference);
    
    // Get first non-dot ID
    for (int i = 0; i < position_occurrences; i++) {
        input = position_in_files[i]->record;
        if (strcmp(".", input->id)) {
            result->id = strdup(input->id);
            break;
        }
    }
    if (result->id == NULL) {
        result->id = strdup(".");
    }
    
    // Calculate weighted mean of the quality
    float accum_quality = 0.0f;
    int total_samples = 0;
    for (int i = 0; i < position_occurrences; i++) {
        file = position_in_files[i]->file;
        input = position_in_files[i]->record;
        
        accum_quality += input->quality * file->num_samples;
        total_samples += file->num_samples;
    }
    result->quality = accum_quality / total_samples;
    
    // TODO concatenate alternates and set their order (used later to assign samples' alleles number)
    size_t max_len = 0, concat_len = 0;
    char *alternate = NULL, *aux = NULL;
    for (int i = 0; i < position_occurrences; i++) {
        file = position_in_files[i]->file;
        input = position_in_files[i]->record;
        
        if (!alternate) {
            alternate = strdup(input->alternate);
            max_len = strlen(input->alternate);
        } else {
            if (!strstr(alternate, input->alternate)) {
                concat_len = strlen(input->alternate);
                aux = realloc(alternate, max_len + concat_len + 1);
                if (aux) {
                    strncat(aux, ",", 1);
                    strncat(aux, input->alternate, concat_len);
                    alternate = aux;
                    max_len += concat_len + 1;
                } else {
                    LOG_FATAL_F("Can't allocate memory for alternate alleles in position %s:%ld\n", 
                                input->chromosome, input->position);
                }
            }
        }
    }
    
    // TODO Get the union of all FORMAT fields
    
    
    return result;
}


array_list_t *get_global_samples(vcf_file_t **files, int num_files) {
    array_list_t *samples = array_list_new(files[0]->num_samples * 2, 2, COLLECTION_MODE_ASYNCHRONIZED);
    for (int i = 0; i < num_files; i++) {
        array_list_insert_all(files[i]->samples_names->items, files[i]->samples_names->size, samples);
    }
    return samples;
}

char *generate_info_field(char **samples, size_t num_samples) {
    // TODO
    return strdup(samples[0]);
}

char *get_empty_sample(int num_format_fields, int gt_pos, enum missing_mode mode) {
    int sample_len = num_format_fields * 2 + 2; // Each field + ':' = 2 chars, except for GT which is 1 char more
    char *sample = (char*) calloc (sample_len, sizeof(char));
    for (int j = 0; j < num_format_fields; j++) {
        if (j > 0) {
            strncat(sample, ":", 1);
        }
        if (j != gt_pos) {
            strncat(sample, ".", 1);
        } else {
            if (mode == MISSING) {
                strncat(sample, "./.", 3);
            } else if (mode == REFERENCE) {
                strncat(sample, "0/0", 3);
            }
        }
    }
    return sample;
}