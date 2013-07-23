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

#include "tdt.h"

int tdt_test(vcf_record_t **variants, int num_variants, family_t **families, int num_families, khash_t(ids) *sample_ids, list_t *output_list) {
    double start = omp_get_wtime();
    
    int ret_code = 0;
    int tid = omp_get_thread_num();
    
    tdt_result_t *result;
    char **sample_data;
    
    int gt_position;
    int father_allele1, father_allele2;
    int mother_allele1, mother_allele2;
    int child_allele1, child_allele2;

    ///////////////////////////////////
    // Perform analysis for each variant

    vcf_record_t *record;
    for (int i = 0; i < num_variants; i++) {
        record = variants[i];
        LOG_DEBUG_F("[%d] Checking variant %.*s:%ld\n", tid, record->chromosome_len, record->chromosome, record->position);
        
        sample_data = (char**) record->samples->items;
        char *format_dup = strndup(record->format, record->format_len);
        gt_position = get_field_position_in_format("GT", format_dup);
        free(format_dup);
        
        // Transmission counts
        int t1 = 0;
        int t2 = 0;
        
        
        // Count over families
        family_t *family;
        for (int f = 0; f < num_families; f++) {
            family = families[f];
            individual_t *father = family->father;
            individual_t *mother = family->mother;

//           LOG_DEBUG_F("[%d] Checking suitability of family %s\n", tid, family->id);
            
            if (father == NULL || mother == NULL) {
                continue;
            }
            
            int father_pos = -1, mother_pos = -1, child_pos = -1;
            
            khiter_t iter = kh_get(ids, sample_ids, father->id);
            if (iter != kh_end(sample_ids)) {
                father_pos = kh_value(sample_ids, iter);
            } else {
                continue;
            }
            
            iter = kh_get(ids, sample_ids, mother->id);
            if (iter != kh_end(sample_ids)) {
                mother_pos = kh_value(sample_ids, iter);
            } else {
                continue;
            }
            
            char *father_sample = strdup(sample_data[father_pos]);
            char *mother_sample = strdup(sample_data[mother_pos]);
            
//           LOG_DEBUG_F("[%d] Samples: Father = %s\tMother = %s\n", tid, father_sample, mother_sample);
            
            // If any parent's alleles can't be read or is missing, go to next family
            if (get_alleles(father_sample, gt_position, &father_allele1, &father_allele2) ||
                get_alleles(mother_sample, gt_position, &mother_allele1, &mother_allele2)) {
                free(father_sample);
                free(mother_sample);
                continue;
            }
            
//           LOG_DEBUG_F("[%d] Alleles: Father = %d/%d\tMother = %d/%d\n", tid, father_allele1, father_allele2, mother_allele1, mother_allele2);
            
            // We need two genotyped parents, with at least one het
            if (father_allele1 == father_allele2 && mother_allele1 == mother_allele2) {
                free(father_sample);
                free(mother_sample);
                continue;
            }
            
            if ((father_allele1 && !father_allele2) || (mother_allele1 && !mother_allele2)) {
                free(father_sample);
                free(mother_sample);
                continue;
            }

//           LOG_DEBUG_F("[%d] Proceeding to analyse family %s...\n", tid, family->id);

            
            int trA = 0;  // transmitted allele from first het parent
            int unA = 0;  // untransmitted allele from first het parent
            
            int trB = 0;  // transmitted allele from second het parent
            int unB = 0;  // untransmitted allele from second het parent
            
            // Consider all offspring in nuclear family
            linked_list_iterator_t *children_iterator = linked_list_iterator_new(family->children);
            individual_t *child = NULL;
            
            while (child = linked_list_iterator_curr(children_iterator)) {
                // Only consider affected children
                if (child->condition != AFFECTED) { 
                    linked_list_iterator_next(children_iterator);
                    continue;
                }
                
                iter = kh_get(ids, sample_ids, child->id);
                if (iter != kh_end(sample_ids)) {
                    child_pos = kh_value(sample_ids, iter);
                } else {
                    linked_list_iterator_next(children_iterator);
                    continue;
                }
                char *child_sample = strdup(sample_data[child_pos]);
    //           LOG_DEBUG_F("[%d] Samples: Child = %s\n", tid, child_sample);
                
                // Skip if offspring has missing genotype
                if (get_alleles(child_sample, gt_position, &child_allele1, &child_allele2)) {
                    free(child_sample);
                    linked_list_iterator_next(children_iterator);
                    continue;
                }
                
                // Exclude mendelian errors
                char *aux_chromosome = strndup(record->chromosome, record->chromosome_len);
                if (check_mendel(aux_chromosome, father_allele1, father_allele2, mother_allele1, mother_allele2, 
                    child_allele1, child_allele2, child->sex)) {
                    free(child_sample);
                    free(aux_chromosome);
                    linked_list_iterator_next(children_iterator);
                    continue;
                }
                free(aux_chromosome);
                
                
                // We've now established: no missing genotypes
                // and at least one heterozygous parent

                // Kid is 00

                if (!child_allele1 && !child_allele2) {
                    if ( ( (!father_allele1) && father_allele2 ) && 
                        ( (!mother_allele1) && mother_allele2 ) )
                    { trA=1; unA=2; trB=1; unB=2; }
                    else 
                    { trA=1; unA=2; } 
                }
                else if ( (!child_allele1) && child_allele2 )  // Kid is 01
                {
                    // het dad
                    if (father_allele1 != father_allele2 )
                    {
                        // het mum
                        if ( mother_allele1 != mother_allele2 )
                    { trA=1; trB=2; unA=2; unB=1; }
                        else if ( !mother_allele1 ) 
                    { trA=2; unA=1; }
                        else { trA=1; unA=2; }
                    }
                    else if ( !father_allele1 ) 
                    {
                        trA=2; unA=1; 
                    }           
                    else
                    {
                        trA=1; unA=2;
                    }
                }
                else // kid is 1/1
                {
                    
                    if ( ( (!father_allele1) && father_allele2 ) && 
                        ( (!mother_allele1) && mother_allele2 ) )
                    { trA=2; unA=1; trB=2; unB=1; }
                    else 
                    { 
                        trA=2; unA=1;
                    }
                }
                
                // We have now populated trA (first transmission) 
                // and possibly trB also 
                
                ////////////////////////////////////////
                // Permutation? 50:50 flip (precomputed)
                
//                 if (permute) {
//                     if (flipA[f])
//                     {
//                     int t = trA;
//                     trA = unA;
//                     unA = t;
//                     
//                     t = trB;
//                     trB = unB;
//                     unB = t;
//                     }
//                 }
                
                // Increment transmission counts
                if (trA==1) { t1++; }
                else if (trA==2) { t2++; }
                
                if (trB==1) { t1++; }
                else if (trB==2) { t2++; }
                
//     //           LOG_DEBUG_F("TDT\t%.*s %s : %d %d - %d %d - %d %d - F %d/%d - M %d/%d - C %d/%d\n", 
//                             record->id_len, record->id, family->id, trA, unA, trB, unB, t1, t2, 
//                             father_allele1, father_allele2, mother_allele1, mother_allele2, child_allele1, child_allele2);
                free(child_sample);
                linked_list_iterator_next(children_iterator);
            } // next offspring in family
            
            linked_list_iterator_free(children_iterator);
            free(father_sample);
            free(mother_sample);
        }  // next nuclear family

        /////////////////////////////
        // Finished counting: now compute
        // the statistics
        
        double tdt_chisq = -1;
        
        // Basic TDT test
        if (t1+t2 > 0) {
            tdt_chisq = ((double) ((t1-t2) * (t1-t2))) / (t1+t2);
        }
        
//         LOG_DEBUG_F("[%d] before adding %s:%ld\n", tid, record->chromosome, record->position);
        result = tdt_result_new(record->chromosome, record->chromosome_len, 
                                record->position, record->id, record->id_len,
                                record->reference, record->reference_len, 
                                record->alternate, record->alternate_len,
                                t1, t2, tdt_chisq);
        list_item_t *output_item = list_item_new(tid, 0, result);
        list_insert_item(output_item, output_list);
//         LOG_DEBUG_F("[%d] after adding %.*s:%ld\n", tid, record->chromosome_len, record->chromosome, record->position);
        
    } // next variant

    double end = omp_get_wtime();
    
    return ret_code;
}


tdt_result_t* tdt_result_new(char *chromosome, int chromosome_len, unsigned long int position, char *id, int id_len, 
                             char *reference, int reference_len, char *alternate, int alternate_len, double t1, double t2, double chi_square) {
    tdt_result_t *result = (tdt_result_t*) malloc (sizeof(tdt_result_t));
    
    result->chromosome = strndup(chromosome, chromosome_len);
    result->position = position;
    result->id = strndup(id, id_len);
    result->reference = strndup(reference, reference_len);
    result->alternate = strndup(alternate, alternate_len);
    result->t1 = t1;
    result->t2 = t2;
    result->odds_ratio = (t2 == 0.0) ? NAN : ((double) t1/t2);
    result->chi_square = chi_square;
    result->p_value = 1 - gsl_cdf_chisq_P(chi_square, 1);
    
    return result;
}

void tdt_result_free(tdt_result_t* result) {
    free(result->chromosome);
    free(result->reference);
    free(result->alternate);
    free(result);
}
