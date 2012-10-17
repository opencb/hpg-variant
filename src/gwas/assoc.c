/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
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

void assoc_test(enum GWAS_task test_type, vcf_record_t **variants, int num_variants, family_t **families, int num_families, 
                cp_hashtable *sample_ids, const void *opt_input, list_t *output_list) {
    int ret_code = 0;
    int tid = omp_get_thread_num();
    int num_samples = cp_hashtable_count(sample_ids);
    
    char **sample_data;
    
    int gt_position;
    int allele1, allele2;

    // Affection counts
    int A1 = 0, A2 = 0, U1 = 0, U2 = 0;
    int num_read = 0, num_analyzed = 0;
    
    ///////////////////////////////////
    // Perform analysis for each variant
    
    vcf_record_t *record;
    for (int i = 0; i < num_variants; i++) {
        record = variants[i];
//         LOG_DEBUG_F("[%d] Checking variant %.*s:%ld\n", tid, record->chromosome_len, record->chromosome, record->position);
        
        A1 = 0; A2 = 0;
        U1 = 0; U2 = 0;
        num_read = 0;
        num_analyzed = 0;

        sample_data = (char**) record->samples->items;
        gt_position = get_field_position_in_format("GT", strndup(record->format, record->format_len));
    
        // Count over families
        family_t *family;
        
        for (int f = 0; f < num_families; f++) {
            family = families[f];
            individual_t *father = family->father;
            individual_t *mother = family->mother;
            cp_list *children = family->children;

//             LOG_DEBUG_F("Read = %d\tAnalyzed = %d\n", num_read, num_analyzed);
//             LOG_DEBUG_F("Family = %d (%s)\n", f, family->id);
            
            // Perform test with father
            if (father != NULL && father->condition != MISSING) {
                num_read++;
                int *father_pos = cp_hashtable_get(sample_ids, father->id);
//                 if (father_pos != NULL) {
//                     LOG_DEBUG_F("[%d] Father %s is in position %d\n", tid, father->id, *father_pos);
//                 } else {
//                     LOG_DEBUG_F("[%d] Father %s is not positioned\n", tid, father->id);
//                     continue;
//                 }
                if (!father_pos) {
                    continue;
                }
                
                char *father_sample = strdup(sample_data[*father_pos]);
                if (!get_alleles(father_sample, gt_position, &allele1, &allele2)) {
                    num_analyzed++;
                    assoc_count_individual(father, record, allele1, allele2, &A1, &A2, &U1, &U2);
                }
                free(father_sample);
            }
            
            // Perform test with mother
            if (mother != NULL && mother->condition != MISSING) {
                num_read++;
                int *mother_pos = cp_hashtable_get(sample_ids, mother->id);
//                 if (mother_pos != NULL) {
//                     LOG_DEBUG_F("[%d] Mother %s is in position %d\n", tid, mother->id, *mother_pos);
//                 } else {
//                     LOG_DEBUG_F("[%d] Mother %s is not positioned\n", tid, mother->id);
//                     continue;
//                 }
                if (!mother_pos) {
                    continue;
                }
                
                char *mother_sample = strdup(sample_data[*mother_pos]);
                if (!get_alleles(mother_sample, gt_position, &allele1, &allele2)) {
                    num_analyzed++;
                    assoc_count_individual(mother, record, allele1, allele2, &A1, &A2, &U1, &U2);
                }
                free(mother_sample);
            }
            
            // Perform test with children
            cp_list_iterator *children_iterator = cp_list_create_iterator(family->children, COLLECTION_LOCK_READ);
            individual_t *child = NULL;
            while ((child = cp_list_iterator_next(children_iterator)) != NULL) {
                int *child_pos = cp_hashtable_get(sample_ids, child->id);
//                 if (child_pos != NULL) {
//                     LOG_DEBUG_F("[%d] Child %s is in position %d\n", tid, child->id, *child_pos);
//                 } else {
//                     LOG_DEBUG_F("[%d] Child %s is not positioned\n", tid, child->id);
//                     continue;
//                 }
                if (!child_pos) {
                    continue;
                }
                
                num_read++;
                char *child_sample = strdup(sample_data[*child_pos]);
                if (!get_alleles(child_sample, gt_position, &allele1, &allele2)) {
                    num_analyzed++;
                    assoc_count_individual(child, record, allele1, allele2, &A1, &A2, &U1, &U2);
                }
                free(child_sample);
            
            } // next offspring in family
            cp_list_iterator_destroy(children_iterator);
        }  // next nuclear family

        /////////////////////////////
        // Finished counting: now compute
        // the statistics
        
        
//         LOG_DEBUG_F("[%d] before adding %.*s:%ld\n", tid, record->chromosome_len, record->chromosome, record->position);
        
        if (test_type == ASSOCIATION_BASIC) {
            double assoc_basic_chisq = assoc_basic_test(A1, U1, A2, U2);
            assoc_basic_result_t *result = assoc_basic_result_new(record->chromosome, record->chromosome_len, 
                                                                  record->position, 
                                                                  record->reference, record->reference_len,
                                                                  record->alternate, record->alternate_len,
                                                                  A1, A2, U1, U2, assoc_basic_chisq);
            list_item_t *output_item = list_item_new(tid, 0, result);
            list_insert_item(output_item, output_list);
        } else if (test_type == FISHER) {
            double p_value = assoc_fisher_test(A1, A2, U1, U2, (double*) opt_input);
            assoc_fisher_result_t *result = assoc_fisher_result_new(record->chromosome, record->chromosome_len, 
                                                                    record->position, 
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
