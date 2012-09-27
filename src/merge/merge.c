#include "merge.h"

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
    vcf_record_t *result = vcf_record_new();
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
        int file_num_samples = get_num_vcf_samples(files[i]);
        if(!strcmp(files[i]->filename, position->file->filename)) {
            // Samples of the file where the position has been read are directly copied
            array_list_insert_all(input->samples->items, file_num_samples, result->samples);
            LOG_DEBUG_F("%d samples of file %s inserted\n", file_num_samples, files[i]->filename);
        } else {
            // Samples in the rest of files are empty and must be filled according to the specified format
            for (int j = 0; j < file_num_samples; j++) {
                array_list_insert(strdup(empty_sample), result->samples);
            }
            LOG_DEBUG_F("file %s\t%d empty samples inserted\n", files[i]->filename, file_num_samples);
        }
    }
    
    // INFO field must be calculated based on statistics about the new list of samples
    // TODO add info-fields as argument (AF, NS and so on)
    result->info = merge_info_field((char**) result->samples->items, result->samples->size);
    
    return result;
}

vcf_record_t *merge_shared_position(vcf_record_file_link **position_in_files, int position_occurrences, 
                                    vcf_file_t **files, int num_files, merge_options_data_t *options, int *info) {
    vcf_file_t *file;
    vcf_record_t *input;
    vcf_record_t *result = vcf_record_new();
    
    // Check consistency among chromosome, position and reference
    for (int i = 0; i < position_occurrences-1; i++) {
        assert(position_in_files[i]);
        assert(position_in_files[i+1]);
        if (strncmp(position_in_files[i]->record->chromosome, position_in_files[i+1]->record->chromosome, position_in_files[i]->record->chromosome_len)) {
            LOG_ERROR("Positions can't be merged: Discordant chromosome\n");
            *info = DISCORDANT_CHROMOSOME;
            return NULL;
        } else if (position_in_files[i]->record->position != position_in_files[i+1]->record->position) {
            LOG_ERROR("Positions can't be merged: Discordant position\n");
            *info = DISCORDANT_POSITION;
            return NULL;
        } else if (strncmp(position_in_files[i]->record->reference, position_in_files[i+1]->record->reference, position_in_files[i]->record->reference_len)) {
            LOG_ERROR("Positions can't be merged: Discordant reference allele\n");
            *info = DISCORDANT_REFERENCE;
            return NULL;
        }
    }
    
    // Once checked, copy constant fields
    result->chromosome = strndup(position_in_files[0]->record->chromosome, position_in_files[0]->record->chromosome_len);
    result->position = position_in_files[0]->record->position;
    result->reference = strndup(position_in_files[0]->record->reference, position_in_files[0]->record->reference_len);
    
    // Get first non-dot ID
    // TODO what can we do when having several ID?
    result->id = merge_id_field(position_in_files, position_occurrences);
    result->id_len = strlen(result->id);
    
    // Calculate weighted mean of the quality
    result->quality = merge_quality_field(position_in_files, position_occurrences);
    
    // Concatenate alternates and set their order (used later to assign samples' alleles number)
    cp_hashtable *alleles_table = cp_hashtable_create(8, cp_hash_istring, (cp_compare_fn) strcasecmp);
    result->alternate = merge_alternate_field(position_in_files, position_occurrences, alleles_table);
    result->alternate_len = strlen(result->alternate);
    
    // Get the union of all FORMAT fields with their corresponding position
    // TODO include INFO fields in FORMAT
    array_list_t *format_fields = array_list_new(16, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    format_fields->compare_fn = strcasecmp;
    result->format = merge_format_field(position_in_files, position_occurrences, format_fields);
    
    // TODO Get indexes for the format in each file
    // TODO Illustrate with an example
    size_t *format_indices = (int**) calloc (num_files, sizeof(int*));
    size_t *file_format_indices;
    int num_fields;
    for (int i = 0; i < num_files; i++) {
        file_format_indices = (size_t*) malloc (format_fields->size * sizeof(size_t));
        char **file_fields = split(strndup(files[i]->format, files[i]->format_len), ":", &num_fields);
        for (int j = 0; j < format_fields->size; j++) {
            size_t idx = array_list_index_of(array_list_get(j, format_fields), file_fields);
            file_format_indices[j] = (idx < ULONG_MAX) ? idx : -1;
        }
        format_indices[i] = file_format_indices;
    }
    
    
    // Create the text for empty samples
    char *format_bak = strdup(result->format);
    int gt_pos = get_field_position_in_format("GT", format_bak);
    char *empty_sample = get_empty_sample(format_fields->size, gt_pos, options->missing_mode);
    size_t empty_sample_len = strlen(empty_sample);
    printf("empty sample = %s\n", empty_sample);
    
    
    // TODO Generate samples using the reordered FORMAT fields and the new alleles' numerical values
    for (int i = 0; i < num_files; i++) {
        vcf_record_file_link *link = NULL;
        for (int j = 0; j < position_occurrences; j++) {
            if (position_in_files[j]->file == files[i]) {
                link = position_in_files[j];
                break;
            }
        }
        
        if (link) {
            vcf_record_t *record = link->record;
            int num_sample_fields;
            char sample[256];
            for (int j = 0; j < get_num_vcf_samples(files[i]); j++) {
                memset(sample, 0, 256 * sizeof(char));
                char **split_sample = split(array_list_get(j, record->samples), ":", &num_sample_fields);
                for (int k = 0; k < format_fields->size; k++) {
                    if (k == gt_pos) {
                        // TODO manage genotype (now it copies contents)
                        size_t idx = format_indices[k];
                        strcat(sample, split_sample[idx]);
                    } else {
                        size_t idx = format_indices[k];
                        strcat(sample, split_sample[idx]);
                    }
                    
                    if (k < format_fields - 1) {
                        strncat(sample, ":", 1);
                    }
                }
                
                printf("sample = %s\n", sample);
            }
            
        } else {
            // If the file has no samples in that position, fill with empty samples
            for (int j = 0; j < get_num_vcf_samples(files[i]); j++) {
                add_vcf_record_sample(empty_sample, empty_sample_len, result);
            }
        }
    }
    
    
    array_list_free(format_fields, free);
    cp_hashtable_destroy(alleles_table);
//     cp_hashtable_destroy(format_fields);
    
    return result;
}


/* ******************************
 *    Field merging functions   *
 * ******************************/

char* merge_id_field(vcf_record_file_link** position_in_files, int position_occurrences) {
    vcf_record_t *input;
    char *result = NULL;
    
    for (int i = 0; i < position_occurrences; i++) {
        input = position_in_files[i]->record;
        if (strncmp(".", input->id, input->id_len)) {
            result = strndup(input->id, input->id_len);
            break;
        }
    }
    
    if (result == NULL) {
        result = strdup(".");
    }
    
    return result;
}

float merge_quality_field(vcf_record_file_link** position_in_files, int position_occurrences) {
    vcf_file_t *file;
    vcf_record_t *input;
    float accum_quality = 0.0f;
    int total_samples = 0;
    int file_num_samples;
    float result;
    
    for (int i = 0; i < position_occurrences; i++) {
        file = position_in_files[i]->file;
        input = position_in_files[i]->record;
        file_num_samples = get_num_vcf_samples(file);
        
        accum_quality += input->quality * file_num_samples;
        total_samples += file_num_samples;
    }
    if (total_samples > 0) {
        result = accum_quality / total_samples;
    } else {
        result = -1;
    }
    
    return result;
}

char* merge_alternate_field(vcf_record_file_link** position_in_files, int position_occurrences, cp_hashtable* alleles_table) {
    vcf_file_t *file;
    vcf_record_t *input;
    size_t max_len = 0, concat_len = 0;
    char *alternates = NULL, *cur_alternate = NULL, *aux = NULL;
    
    int *allele_index = (int*) calloc (1, sizeof(int)); *allele_index = 0;
    cp_hashtable_put(alleles_table,
                     strndup(position_in_files[0]->record->reference, position_in_files[0]->record->reference_len), 
                     allele_index);
    
    for (int i = 0; i < position_occurrences; i++) {
        file = position_in_files[i]->file;
        input = position_in_files[i]->record;
        
        cur_alternate = strndup(input->alternate, input->alternate_len);
        
        if (!alternates) {
            alternates = cur_alternate;
            max_len = input->alternate_len;
            
            allele_index = (int*) calloc (1, sizeof(int)); *allele_index = 1;
            cp_hashtable_put(alleles_table, strndup(cur_alternate, input->alternate_len), allele_index);
        } else {
            if (!strstr(alternates, cur_alternate)) {
                concat_len = input->alternate_len;
                aux = realloc(alternates, max_len + concat_len + 1);
                if (aux) {
                    // Concatenate alternate value to the existing list
                    strncat(aux, ",", 1);
                    strncat(aux, cur_alternate, concat_len);
                    alternates = aux;
                    max_len += concat_len + 1;
                    
                    // In case the allele is not in the hashtable, insert it with a new index
//                     if (!cp_hashtable_contains(alleles_table, cur_alternate)) {
                        allele_index = (int*) calloc (1, sizeof(int)); (*allele_index)++;
                        cp_hashtable_put(alleles_table, cur_alternate, allele_index);
//                     }
                } else {
                    LOG_FATAL_F("Can't allocate memory for alternate alleles in position %s:%ld\n", 
                                input->chromosome, input->position);
                }
            }
        }
    }

    return alternates;
}


char *merge_info_field(char **samples, size_t num_samples) {
    // TODO
    return strdup(samples[0]);
}

char* merge_format_field(vcf_record_file_link** position_in_files, int position_occurrences, array_list_t* format_fields) {
    char *result;
    int format_text_len = 0;
    int *field_index = (int*) calloc (1, sizeof(int)); *field_index = 0;
    // TODO doesn't work, must be checked against a record of each file
//     for (int i = 0; i < num_files; i++) {
//         int num_fields;
//         char **fields = split(strndup(files[i]->format, files[i]->format_len), ":", &num_fields);
//         for (int j = 0; j < num_fields; j++) {
//             if (!array_list_contains(fields[j], format_fields)) {
//                 array_list_insert(fields[j], format_fields);
//                 format_text_len += strlen(fields[j]) + 1; // concat field + ":"
//             }
//         }
//     }
    result = calloc (format_text_len + 1, sizeof(char));
    strcat(result, (char*) array_list_get(0, format_fields));
    for (int i = 1; i < format_fields->size; i++) {
        strncat(result, ":", 1);
        strcat(result, (char*) array_list_get(i, format_fields));
    }
    
    return result;
}


/* ******************************
 *      Auxiliary functions     *
 * ******************************/

array_list_t *get_global_samples(vcf_file_t **files, int num_files) {
    array_list_t *samples = array_list_new(get_num_vcf_samples(files[0]) * 2, 2, COLLECTION_MODE_ASYNCHRONIZED);
    for (int i = 0; i < num_files; i++) {
        array_list_insert_all(files[i]->samples_names->items, files[i]->samples_names->size, samples);
    }
    return samples;
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
