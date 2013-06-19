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
#include <malloc.h>

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
        output_prefix = (char*) calloc (12 + record->chromosome_len, sizeof(char));
        strncat(output_prefix, "chromosome_", 11);
        strncat(output_prefix, record->chromosome, record->chromosome_len);
        split_result = new_split_result(vcf_record_copy(record), output_prefix);
        
        // Insert results in output list
        list_item_t *item = list_item_new(i, 0, split_result);
        list_insert_item(item, output_list);
    }
    
    return 0;
}

int split_by_coverage(vcf_record_t **variants, int num_variants, long *intervals, int num_intervals, list_t* output_list) {
    char *output_prefix;
    vcf_record_t *record;
    split_result_t *split_result;
    
    // Get length of the string containing the value of the last (and greatest) interval
    char limit_buf[64];
    sprintf(limit_buf, "%ld", intervals[num_intervals-1]);
    size_t last_limit_len = strlen(limit_buf);
    
    // For each variant, its output filename will be 'coverage_<#limit1>_<#limit2>_<original_filename>.vcf'
    for (int i = 0; i < num_variants; i++) {
        record = variants[i];
        output_prefix = (char*) calloc (11 + last_limit_len * 2, sizeof(char));
	
        char *info = strndup(record->info, record->info_len);
        char *value_str = get_field_value_in_info("DP", info);
        int value = atoi(value_str);
        free(info);
        
        long limit_lo = 0, limit_hi = 0;
        
        if (value <= intervals[0]) {
            // Within first interval
            limit_hi = intervals[0];
        } else {
            int interval_found = 0;
            // Within intermediate interval
            for (int j = 1; j < num_intervals; j++) {
                if (value <= intervals[j]) {
                    limit_lo = intervals[j-1];
                    limit_hi = intervals[j];
                    interval_found = 1;
                    break;
                }
            }
            
            // Within last interval
            if (!interval_found) {
                limit_lo = intervals[num_intervals-1];
                limit_hi = LONG_MAX;
            }
        }
        
/*
        free(value_str);
*/
        
        if (limit_hi < LONG_MAX) {
            sprintf(output_prefix, "coverage_%ld_%ld", limit_lo, limit_hi);
        } else {
            sprintf(output_prefix, "coverage_%ld_N", limit_lo);
        }

        split_result = new_split_result(vcf_record_copy(record), output_prefix);

        // Insert results in output list
        list_item_t *item = list_item_new(i, 0, split_result);
        list_insert_item(item, output_list);
    }
    
    return 0;
}