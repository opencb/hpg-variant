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

#include "assoc.h"

void assoc_test(enum ASSOC_task test_type, vcf_record_t **variants, int num_variants, individual_t **samples, int num_samples,
                const void *opt_input, list_t *output_list) {
    int tid = omp_get_thread_num();

    vcf_record_t *record;
    individual_t *individual;
    char *sample_data, *format;
    
    int gt_position;
    int allele1, allele2;

    // Affection counts
    int A1 = 0, A2 = 0, U1 = 0, U2 = 0;
    
    // Perform analysis for each variant
    for (int i = 0; i < num_variants; i++) {
        record = variants[i];
//         LOG_DEBUG_F("[%d] Checking variant %.*s:%ld\n", tid, record->chromosome_len, record->chromosome, record->position);
        
        A1 = 0; A2 = 0;
        U1 = 0; U2 = 0;

        format = strndup(record->format, record->format_len);
        gt_position = get_field_position_in_format("GT", format);
        free(format);
    
        // Count over individuals
        for (int j = 0; j < num_samples; j++) {
        	individual = samples[j];
        	sample_data = strdup(array_list_get(j, record->samples));
        	if (get_alleles(sample_data, gt_position, &allele1, &allele2) == ALLELES_OK) {
                    assoc_count_individual(individual, record, allele1, allele2, &A1, &A2, &U1, &U2);
                }
        	free(sample_data);
        }
        
        // Finished counting: now compute the statistics
        if (test_type == CHI_SQUARE) {
            double assoc_basic_chisq = assoc_basic_test(A1, U1, A2, U2);
            assoc_basic_result_t *result = assoc_basic_result_new(record->chromosome, record->chromosome_len, 
                                                                  record->position, record->id, record->id_len, 
                                                                  record->reference, record->reference_len,
                                                                  record->alternate, record->alternate_len,
                                                                  A1, A2, U1, U2, assoc_basic_chisq);
            list_item_t *output_item = list_item_new(tid, 0, result);
            list_insert_item(output_item, output_list);
        } else if (test_type == FISHER) {
            double p_value = assoc_fisher_test(A1, A2, U1, U2, (double*) opt_input);
            assoc_fisher_result_t *result = assoc_fisher_result_new(record->chromosome, record->chromosome_len, 
                                                                    record->position, record->id, record->id_len, 
                                                                    record->reference, record->reference_len,
                                                                    record->alternate, record->alternate_len,
                                                                    A1, A2, U1, U2, p_value);
            list_item_t *output_item = list_item_new(tid, 0, result);
            list_insert_item(output_item, output_list);
        }
        
//         LOG_DEBUG_F("[%d] after adding %.*s:%ld\n", tid, record->chromosome_len, record->chromosome, record->position);
        
    } // next variant

}


void assoc_count_individual(individual_t *individual, vcf_record_t *record, int allele1, int allele2, 
                           int *affected1, int *affected2, int *unaffected1, int *unaffected2) {
    int A1 = 0, A2 = 0, A0 = 0;
    int U1 = 0, U2 = 0, U0 = 0;
    
    assert(individual);
    
    if (!strncmp("X", record->chromosome, record->chromosome_len)) {
        if (individual->condition == AFFECTED) { // if affected 
            if (!allele1 && !allele2) {
                A1++;
            } else if (allele1 && allele2) {
                A2++;
            }
        } else if (individual->condition == UNAFFECTED) { // unaffected if not missing
            if (!allele1 && !allele2) {
                U1++;
            } else if (allele1 && allele2) {
                U2++;
            }
        }
    } else {
        if (individual->condition == AFFECTED) { // if affected
            if (!allele1 && !allele2) {
                A1 += 2;
            } else if (allele1 && allele2) {
                A2 += 2;
            } else if (allele1 != allele2) {
                A1++; A2++;
            }
        } else if (individual->condition == UNAFFECTED) { // unaffected if not missing
            if (!allele1 && !allele2) {
                U1 += 2;
            } else if (allele1 && allele2) {
                U2 += 2;
            } else if (allele1 != allele2) {
                U1++; U2++;
            }
        }
          
    }
    
    // Set output values
    *affected1 += A1;
    *affected2 += A2;
    *unaffected1 += U1;
    *unaffected2 += U2;
}
