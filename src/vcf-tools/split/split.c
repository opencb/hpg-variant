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

int split_by_coverage(vcf_record_t **variants, int num_variants, int *intervals, int num_intervals, list_t* output_list) {
    char *output_prefix;
    vcf_record_t *record;
    split_result_t *split_result;
    
    // For each variant, its output filename will be 'coverage_<#limit1>_<#limit2>_<original_filename>.vcf'
    for (int i = 0; i < num_variants; i++) {
        record = variants[i];
	//char buffer[64]
	char* buffer = (char*)malloc (64 * sizeof(char));
	int size_prefix;
	//size_prefix = sprintf(buffer, "coverage_%i_%i", intervals[0], intervals[1]) + 1;
        size_prefix = 25;
	output_prefix = (char*) malloc (size_prefix*sizeof(char));
        strncat(output_prefix, "coverage_", 9);
	
	char* value_str = get_field_value_in_info("DP", strndup(record->info, record->info_len));
	int value = atoi(value_str);
	
	if (value != NULL){  
	  int limit1 = 0;
	  int aux_index = num_intervals-1;
	  int limit2 = intervals[aux_index];
	  int aux_limit1 = intervals[0];
	  int aux_limit2 = intervals[aux_index];
	  
	  printf("\nBEFORE limit1 = %i limint2 = %i\n", limit1,limit2);
	  
	  //get limits
	  for (int j = 0; j < num_intervals; j++){
	    
	    if ( value > aux_limit1){
		limit1 = aux_limit1;
	    }
	      
	    if (value < aux_limit2){
	      limit2 = aux_limit2;
	    }
	    
	    if (j != num_intervals - 1){
		aux_limit1 = intervals[j+1];
	    }
	    if (aux_index != 0){
	      aux_limit2 = intervals[aux_index-1];
	      aux_index--;
	    }
	    
	    // if both limits are equal to intervals[num_intervals - 1] the output filename must be "coverage_<#limit1>_N_<original_filename>.vcf"
	  }//end for
	  
	  printf("AFTER limit1 = %i limint2 = %i\n", limit1,limit2);
	  
	  strncat(output_prefix, limit1, sizeof(limit1));
	  strncat(output_prefix, "_", 1);
	  if (limit1 != limit2){
	    strncat(output_prefix, limit2, sizeof(limit2));
	  }
	  else{
	    strncat(output_prefix, "N", 1);
	  }
	  
	  printf("output_prefix = %s\n", output_prefix);
	  
	  split_result = new_split_result(vcf_record_copy(record), output_prefix);
	  
	  // Insert results in output list
	  list_item_t *item = list_item_new(i, 0, split_result);
	  list_insert_item(item, output_list);
	  }//end if (value)
	  
	else{
	  return SPLIT_OPTION_DOES_NOT_EXIST;
	}//end else
    }//end for (variant)
    
    return 0;
}//end function