/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * hpg-variant is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-variant. If not, see <http://www.gnu.org/licenses/>.
 */

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
