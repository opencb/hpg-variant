#include "assoc.h"

void prepare_assoc_counters(ped_file_t *ped_file, list_item_t *variants, int num_variants, cp_hashtable *sample_ids, int *counters) {
    int tid = omp_get_thread_num();
    cp_hashtable *families = ped_file->families;
    char **families_keys = (char**) cp_hashtable_get_keys(families);
    int num_families = get_num_families(ped_file);
    int num_samples = cp_hashtable_count(sample_ids);
    
    char **sample_data;
    
    int gt_position;
    int allele1, allele2;

    ///////////////////////////////////
    // Perform analysis for each variant

    list_item_t *cur_variant = variants;
    for (int i = 0; i < num_variants && cur_variant != NULL; i++) {
        vcf_record_t *record = (vcf_record_t*) cur_variant->data_p;
        LOG_DEBUG_F("[%d] Checking variant %s:%ld\n", tid, record->chromosome, record->position);
        
        int num_read = 0, number_analyzed = 0;
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

            printf("Read = %d\tAnalyzed = %d\n", num_read, number_analyzed);
            printf("Family = %d (%s)\n", f, family->id);
            
            // Perform test with father
            if (father != NULL && father->condition != MISSING) {
                num_read++;
                int *father_pos = cp_hashtable_get(sample_ids, father->id);
                if (father_pos != NULL) {
                    LOG_DEBUG_F("[%d] Father %s is in position %d\n", tid, father->id, *father_pos);
                } else {
                    LOG_DEBUG_F("[%d] Father %s is not positioned\n", tid, father->id);
                    continue;
                }
                
                char *father_sample = sample_data[*father_pos];
                if (!get_alleles(father_sample, gt_position, &allele1, &allele2)) {
                    number_analyzed++;
                    assoc_count_individual(father, record, allele1, allele2, &A1, &A2, &U1, &U2);
                } else {
                    printf("[%s] Father sample = %s\n", family->id, father_sample);
                }
            }
            
            // Perform test with mother
            if (mother != NULL && mother->condition != MISSING) {
                num_read++;
                int *mother_pos = cp_hashtable_get(sample_ids, mother->id);
                if (mother_pos != NULL) {
                    LOG_DEBUG_F("[%d] Mother %s is in position %d\n", tid, mother->id, *mother_pos);
                } else {
                    LOG_DEBUG_F("[%d] Mother %s is not positioned\n", tid, mother->id);
                    continue;
                }
                
                char *mother_sample = sample_data[*mother_pos];
                if (!get_alleles(mother_sample, gt_position, &allele1, &allele2)) {
                    number_analyzed++;
                    assoc_count_individual(mother, record, allele1, allele2, &A1, &A2, &U1, &U2);
                } else {
                    printf("[%s] Mother sample = %s\n", family->id, mother_sample);
                }
            }
            
            // Perform test with children
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
                
                num_read++;
                char *child_sample = sample_data[*child_pos];
                if (!get_alleles(child_sample, gt_position, &allele1, &allele2)) {
                    number_analyzed++;
                    assoc_count_individual(child, record, allele1, allele2, &A1, &A2, &U1, &U2);
                } else {
                    printf("[%s] Child sample = %s\n", family->id, child_sample);
                }
            
            } // next offspring in family
            cp_list_iterator_destroy(children_iterator);
        }  // next nuclear family
        
        counters[i * 4]     = A1;
        counters[i * 4 + 1] = A2;
        counters[i * 4 + 2] = U1;
        counters[i * 4 + 3] = U2;
    }

}

void assoc_count_individual(individual_t *individual, vcf_record_t *record, int allele1, int allele2, 
                           int *affected1, int *affected2, int *unaffected1, int *unaffected2) {
    int A1 = 0, A2 = 0, A0 = 0;
    int U1 = 0, U2 = 0, U0 = 0;
    
    if (!strcmp("X", record->chromosome)) {
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
