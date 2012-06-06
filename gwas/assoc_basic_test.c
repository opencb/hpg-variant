#include "assoc_basic_test.h"

int assoc_basic_test(ped_file_t *ped_file, list_item_t *variants, int num_variants, cp_hashtable *sample_ids, list_t *output_list) {
    int ret_code = 0;
    int tid = omp_get_thread_num();
    cp_hashtable *families = ped_file->families;
    char **families_keys = (char**) cp_hashtable_get_keys(families);
    int num_families = get_num_families(ped_file);
    int num_samples = cp_hashtable_count(sample_ids);
    
    assoc_basic_result_t *result;
    
    char **sample_data;
    
    int gt_position;
    int allele1, allele2;

    ///////////////////////////////////
    // Perform analysis for each variant

    list_item_t *cur_variant = variants;
    for (int i = 0; i < num_variants && cur_variant != NULL; i++) {
        vcf_record_t *record = (vcf_record_t*) cur_variant->data_p;
        LOG_DEBUG_F("[%d] Checking variant %s:%ld\n", tid, record->chromosome, record->position);
        
        // TODO implement arraylist in order to avoid this conversion
        sample_data = (char**) list_to_array(record->samples);
        gt_position = get_field_position_in_format("GT", record->format);
    
        // Affection counts
        int A1 = 0;
        int A2 = 0;
        int U1 = 0;
        int U2 = 0;

        // Count over families
        family_t *family;
        
        for (int f = 0; f < num_families; f++) {
            family = cp_hashtable_get(families, families_keys[f]);
            individual_t *father = family->father;
            individual_t *mother = family->mother;
            cp_list *children = family->children;

            // TODO perform test with father
            if (father != NULL) {
                int *father_pos = cp_hashtable_get(sample_ids, father->id);
                if (father_pos != NULL) {
                    LOG_DEBUG_F("[%d] Father %s is in position %d\n", tid, father->id, *father_pos);
                } else {
                    LOG_DEBUG_F("[%d] Father %s is not positioned\n", tid, father->id);
                    continue;
                }
                
                char *father_sample = sample_data[*father_pos];
                LOG_DEBUG_F("[%d] Father sample = %s\n", tid, father_sample);
                
                if (get_alleles(father_sample, gt_position, &allele1, &allele2)) {
                    continue;
                } else {
                    assoc_basic_individual(father, record, allele1, allele2, &A1, &A2, &U1, &U2);
                }
            }
            
            // TODO perform test with mother
            if (mother != NULL) {
                int *mother_pos = cp_hashtable_get(sample_ids, mother->id);
                if (mother_pos != NULL) {
                    LOG_DEBUG_F("[%d] Mother %s is in position %d\n", tid, mother->id, *mother_pos);
                } else {
                    LOG_DEBUG_F("[%d] Mother %s is not positioned\n", tid, mother->id);
                    continue;
                }
                
                char *mother_sample = sample_data[*mother_pos];
                LOG_DEBUG_F("[%d] Mother sample = %s\n", tid, mother_sample);
                
                if (get_alleles(mother_sample, gt_position, &allele1, &allele2)) {
                    continue;
                } else {
                    assoc_basic_individual(mother, record, allele1, allele2, &A1, &A2, &U1, &U2);
                }
            }
            
            // TODO perform test with children
            cp_list_iterator *children_iterator = cp_list_create_iterator(family->children, COLLECTION_LOCK_READ);
            individual_t *child = NULL;
            while ((child = cp_list_iterator_next(children_iterator)) != NULL) {
                int *child_pos = cp_hashtable_get(sample_ids, child->id);
                if (child_pos != NULL) {
                    LOG_DEBUG_F("[%d] Child %s is in position %d\n", tid, child->id, *child_pos);
                } else {
                    LOG_DEBUG_F("[%d] Child %s is not positioned\n", tid, child->id);
                    continue;
                }
                
                char *child_sample = sample_data[*child_pos];
                LOG_DEBUG_F("[%d] Child sample = %s\n", tid, child_sample);
                
                if (get_alleles(child_sample, gt_position, &allele1, &allele2)) {
                    continue;
                } else {
                    assoc_basic_individual(child, record, allele1, allele2, &A1, &A2, &U1, &U2);
                }
            
            } // next offspring in family
            cp_list_iterator_destroy(children_iterator);
        
        }  // next nuclear family

        /////////////////////////////
        // Finished counting: now compute
        // the statistics
        
        double assoc_basic_chisq = chi_square(A1, U1, A2, U2);
        
        LOG_DEBUG_F("[%d] before adding %s:%ld\n", tid, record->chromosome, record->position);
        result = assoc_basic_result_new(record->chromosome, record->position, record->reference, record->alternate, 
                                        A1, A2, U1, U2, 
                                        assoc_basic_chisq);
        list_item_t *output_item = list_item_new(tid, 0, result);
        list_insert_item(output_item, output_list);
        LOG_DEBUG_F("[%d] after adding %s:%ld\n", tid, record->chromosome, record->position);
        
        cur_variant = cur_variant->next_p;
        
        // Free samples
        // TODO implement arraylist in order to avoid this code
        free(sample_data);
    } // next variant

    // Free families' keys
//     free(families_keys);
    
    return ret_code;
}


void assoc_basic_individual(individual_t *individual, vcf_record_t *record, int allele1, int allele2, 
                           int *affected1, int *affected2, int *unaffected1, int *unaffected2) {
    int A1 = 0, A2 = 0, A0 = 0;
    int U1 = 0, U2 = 0, U0 = 0;
    
//     printf("individual %s, condition %d\n", individual->id, individual->condition);
           
    if (!strcmp("X", record->chromosome)) {
        if (individual->condition == AFFECTED) { // if affected 
            if (!allele1 && !allele2) {
                A1++;
            } else if (allele1 && allele2) {
                A2++;
            } /*else {
                A0++;
            }*/
//             if (!allele1) {
//                 if (!allele2) {  // 0|0
//                     A1++;
//                 }
//             } else {
//                 if (allele2) { 
//                     A2++;
//                 } else {
//                     A0++;
//                 }
//             }
        } else if (individual->condition == UNAFFECTED) { // unaffected if not missing
            if (!allele1 && !allele2) {
                U1++;
            } else if (allele1 && allele2) {
                U2++;
            } /*else {
                U0++;
            }*/
//             if (!allele1) {
//                 if (!allele2) {  
//                     U1++;
//                 }
//             } else {
//                 if (allele2) {
//                     U2++;
//                 } else {
//                     U0++;
//                 }
//             }
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
//             if (!allele1) {
//                 if (!allele2) {  
//                     A1+=2;
//                 } else { 
//                     A1++; A2++;
//                 }
//             } else {
//                 if (allele2) {  
//                     A2+=2;
//                 } else {
//                     A0+=2;
//                 }
//             }
        } else if (individual->condition == UNAFFECTED) { // unaffected if not missing
//             printf("unaffected = %s\n", individual->id);
            if (!allele1 && !allele2) {
                U1 += 2;
            } else if (allele1 && allele2) {
                U2 += 2;
            } else if (allele1 != allele2) {
                U1++; U2++;
            }
//             if (!allele1)  {
//                 if (!allele2) {
//                     U1+=2;
//                 } else { 
//                     U1++; U2++;
//                 } 
//             } else {
//                 if (allele2) {  
//                     U2+=2;
//                 } else {
//                     U0+=2;
//                 }
//             }
        }
          
    }
    
//     if (!strcmp("3737", individual->id)) {
//         printf("[%s] A1 = %d, U1 = %d, A2 = %d, U2 = %d\n", individual->id, A1, U1, A2, U2);
//     }
    
    // Set output values
    *affected1 += A1;
    *affected2 += A2;
    *unaffected1 += U1;
    *unaffected2 += U2;
}


double chi_square(int a, int b, int c, int d) {
    double total_alleles = a + c + b + d;
    
    double total_affected = a + c;
    double total_unaffected = b + d;
    
    double total_allele1 = a + b;
    double total_allele2 = c + d;
    
    double expected_affected_allele1    = (total_affected   * total_allele1) / total_alleles;
    double expected_affected_allele2    = (total_affected   * total_allele2) / total_alleles;
    double expected_unaffected_allele1  = (total_unaffected * total_allele1) / total_alleles;
    double expected_unaffected_allele2  = (total_unaffected * total_allele2) / total_alleles;

    return ((a - expected_affected_allele1)   * (a - expected_affected_allele1))   / expected_affected_allele1 + 
           ((c - expected_affected_allele2)   * (c - expected_affected_allele2))   / expected_affected_allele2 +
           ((b - expected_unaffected_allele1) * (b - expected_unaffected_allele1)) / expected_unaffected_allele1 + 
           ((d - expected_unaffected_allele2) * (d - expected_unaffected_allele2)) / expected_unaffected_allele2 ;
               
}


assoc_basic_result_t* assoc_basic_result_new(char *chromosome, unsigned long int position, char *reference, char *alternate, 
                                             int affected1, int affected2, int unaffected1, int unaffected2, double chi_square) {
    assoc_basic_result_t *result = (assoc_basic_result_t*) malloc (sizeof(assoc_basic_result_t));
    
    result->chromosome = strdup(chromosome);
    result->position = position;
    result->reference = strdup(reference);
    result->alternate = strdup(alternate);
    result->affected1 = affected1;
    result->affected2 = affected2;
    result->unaffected1 = unaffected1;
    result->unaffected2 = unaffected2;
    result->odds_ratio = (affected2 == 0 || unaffected1 == 0) ? NAN : 
                         ((double) affected1 / affected2) * ((double) unaffected2 / unaffected1);
    result->chi_square = chi_square;
    
    return result;
}

void assoc_basic_result_free(assoc_basic_result_t* result) {
    free(result->chromosome);
    free(result->reference);
    free(result->alternate);
    free(result);
}
