#include "model.h"


/* **************************
 *       Main pipeline      *
 * **************************/

risky_combination *get_model_from_combination_in_fold(int order, int comb[order], uint8_t **val, unsigned int num_affected_in_training, unsigned int num_unaffected_in_training,
                                                      int num_genotype_combinations, uint8_t **genotype_combinations, int num_counts, masks_info info, double *masks_time, double *counts_time) {
    risky_combination *risky_comb = NULL;
    
    // Get counts for the provided genotypes
    double start_masks = omp_get_wtime();
    uint8_t *masks = get_masks(order, val, num_affected_in_training, num_unaffected_in_training, info); // Grouped by SNP
    *masks_time += omp_get_wtime() - start_masks;
    
//     printf("masks (%d) = {\n", info.num_masks);
//     printf("%d ", masks[0]);
//     for (int i = 1; i < info.num_masks; i++) {
//         if (i % info.num_samples_per_mask == 0) {
//             printf("\n");
//         }
//         printf("%d ", masks[i]);
//     }
//     printf("}\n");
    
    double start_counts = omp_get_wtime();
    int *counts = get_counts(order, masks, genotype_combinations, num_genotype_combinations, 
                             info.num_affected_with_padding, info.num_unaffected_with_padding, info.num_samples_per_mask, num_counts);
    *counts_time += omp_get_wtime() - start_counts;
    
//     _mm_free(masks);
    
//     printf("counts = {\n");
//     for (int j = 0; j < 3; j++) {
//         printf("  ");
//         for (int k = 0; k < 6; k++) {
//             printf("%d ", counts[j * 6 + k]);
//         }
//         printf("\n");
//     }
//     printf("}\n");
    
    // Get high risk pairs for those counts
    void *aux_info;
    int num_risky;
    int *risky_idx = get_high_risk_combinations(counts, num_counts, num_affected_in_training, num_unaffected_in_training, 
                                                &num_risky, &aux_info, mdr_high_risk_combinations);
    
    // Filter non-risky SNP combinations
    if (num_risky > 0) {
        // Put together the info about the SNP combination and its genotype combinations
        risky_comb = risky_combination_new(order, comb, genotype_combinations, num_risky, risky_idx, aux_info);
        
//         printf("risky combination = {\n  SNP: ");
//         print_combination(risky_comb->combination, 0, order);
//         printf("  GT: ");
//         for (int j = 0; j < num_risky * 2; j++) {
//             if (j % 2) {
//                 printf("%d), ", risky_comb->genotypes[j]);
//             } else {
//                 printf("(%d ", risky_comb->genotypes[j]);
//             }
//         }
//         printf("\n}\n");
    }
    
    free(counts);
    free(risky_idx);
    
    return risky_comb;
}


double test_model(int order, risky_combination *risky_comb, uint8_t **val, 
                  unsigned int num_affected, unsigned int num_unaffected, double *confusion_time) {
    // Step 5 -> Check against a testing partition
    // Get the matrix containing {FP,FN,TP,TN}
    double start_conf = omp_get_wtime();
    unsigned int *confusion_matrix = get_confusion_matrix(order, risky_comb, num_affected, num_unaffected, val);
    *confusion_time += omp_get_wtime() - start_conf;
    
//     printf("confusion matrix = { ");
//     for (int k = 0; k < 4; k++) {
//         printf("%u ", confusion_matrix[k]);
//     }
//     printf("}\n");
    
    // Evaluate the model, basing on the confusion matrix
    double eval = evaluate_model(confusion_matrix, BA);
    
//     printf("risky combination = {\n  SNP: ");
//     print_combination(risky_comb->combination, 0, order);
//     printf("  GT: ");
//     for (int j = 0; j < risky_comb->num_risky * 2; j++) {
//         if (j % 2) {
//             printf("%d), ", risky_comb->genotypes[j]);
//         } else {
//             printf("(%d ", risky_comb->genotypes[j]);
//         }
//     }
//     printf("\n}\n", eval);
    
    risky_comb->accuracy = eval;
    
    free(confusion_matrix);
    
    return eval;
}


int add_to_model_ranking(risky_combination *risky_comb, int max_ranking_size, linked_list_t *ranking_risky) {
    // Step 6 -> Ellaborate a ranking of the best N combinations
    risky_combination *last_element = (linked_list_size(ranking_risky) > 0) ? linked_list_get_last(ranking_risky) : NULL;
    size_t current_ranking_size = ranking_risky->size;
    
    linked_list_iterator_t* iter = linked_list_iterator_new(ranking_risky);
    risky_combination *element = NULL;
//     printf("Ranking (size %zu) = { ", current_ranking_size);
//     while(element = linked_list_iterator_next(iter)) {
//         printf("(%d %d - %.3f) ", element->combination[0], element->combination[1], element->accuracy);
//     }
//     printf("}\n");
//     linked_list_iterator_first(iter);
    
    if (current_ranking_size > 0) {
        if (last_element) {
            LOG_DEBUG_F("To insert %.3f\tRanking's last is %.3f\n", risky_comb->accuracy, last_element->accuracy);
        } else {
            LOG_DEBUG_F("To insert %.3\n", risky_comb->accuracy );
        }
        
        // If accuracy is not greater than the last element, don't bother inserting
        if (risky_comb->accuracy > last_element->accuracy) {
            int position = 0;
            while (element = linked_list_iterator_curr(iter)) {
                LOG_DEBUG_F("To insert %.3f\tIn ranking (pos #%d) %.3f\n", risky_comb->accuracy, position, element->accuracy);
                if (risky_comb->accuracy > element->accuracy) {
                    linked_list_iterator_insert(risky_comb, iter);
                    
                    if (current_ranking_size >= max_ranking_size) {
                        linked_list_iterator_last(iter);
                        risky_combination *removed = linked_list_iterator_remove(iter);
                        risky_combination_free(removed);
                    }
                    
                    linked_list_iterator_free(iter);
                    return position;
                }
                element = linked_list_iterator_next(iter);
                position++;
            }
        }
        
        if (current_ranking_size < max_ranking_size) {
            LOG_DEBUG_F("To insert %.3f at the end", risky_comb->accuracy);
            linked_list_insert_last(risky_comb, ranking_risky);
            linked_list_iterator_free(iter);
            return ranking_risky->size - 1;
        }
    } else {
        linked_list_insert_last(risky_comb, ranking_risky);
        linked_list_iterator_free(iter);
        return ranking_risky->size - 1;
    }
    
    linked_list_iterator_free(iter);
    
    return -1;
}


/* **************************
 *          Counts          *
 * **************************/

int* get_counts(int order, uint8_t *masks, uint8_t **genotype_combinations, int num_genotype_combinations, 
                int num_affected, int num_unaffected, int num_samples_per_mask, int num_counts) {
    // TODO now that num_affected, num_unaffected are padded, num_samples_per_mask should not be neccesary
    int *counts = malloc(num_counts * sizeof(int)); // Affected and unaffected
    
    uint8_t *comb;
    int flag = 1, count = 0;
    
    __m128i snp_and, snp_cmp;
    
    for (int c = 0; c < num_genotype_combinations; c++) {
        comb = genotype_combinations[c];
//         print_gt_combination(comb, c, order);
        
        count = 0;
       
        for (int i = 0; i < num_affected; i += 16) {
            // Aligned loading
            snp_and = _mm_load_si128(masks + comb[0] * num_samples_per_mask + i);
            
            // Perform AND operation with all SNPs in the combination
            for (int j = 0; j < order; j++) {
                snp_cmp = _mm_load_si128(masks + j * NUM_GENOTYPES * num_samples_per_mask + comb[j] * num_samples_per_mask + i);
                snp_and = _mm_and_si128(snp_and, snp_cmp);
            }
            
            count += _mm_popcnt_u64(_mm_extract_epi64(snp_and, 0)) + 
                     _mm_popcnt_u64(_mm_extract_epi64(snp_and, 1));
        }
        
        LOG_DEBUG_F("aff comb idx (%d) = %d, %d\n", c * 2, count);
        counts[c * 2] = count;
        
        count = 0;
        
        for (int i = 0; i < num_unaffected; i += 16) {
            // Aligned loading
            snp_and = _mm_load_si128(masks + comb[0] * num_samples_per_mask + num_affected + i);
            
            // Perform AND operation with all SNPs in the combination
            for (int j = 0; j < order; j++) {
                snp_cmp = _mm_load_si128(masks + j * NUM_GENOTYPES * num_samples_per_mask + comb[j] * num_samples_per_mask + num_affected + i);
                snp_and = _mm_and_si128(snp_and, snp_cmp);
            }
            
            count += _mm_popcnt_u64(_mm_extract_epi64(snp_and, 0)) + 
                     _mm_popcnt_u64(_mm_extract_epi64(snp_and, 1));
        }
        
        LOG_DEBUG_F("unaff comb idx (%d) = %d\n", c * 2 + 1, count);
        counts[c * 2 + 1] = count;
    }
    
    return counts;
}

uint8_t* get_masks(int order, uint8_t **genotypes, int num_affected, int num_unaffected, masks_info info) {
    /* 
     * Structure: Genotypes of a SNP in each 'row'
     * 
     * SNP(0) - Mask genotype 0 (all samples)
     * SNP(0) - Mask genotype 1 (all samples)
     * SNP(0) - Mask genotype 2 (all samples)
     * 
     * SNP(1) - Mask genotype 0 (all samples)
     * SNP(1) - Mask genotype 1 (all samples)
     * SNP(1) - Mask genotype 2 (all samples)
     * 
     * ...
     * 
     * SNP(order-1) - Mask genotype 0 (all samples)
     * SNP(order-1) - Mask genotype 1 (all samples)
     * SNP(order-1) - Mask genotype 2 (all samples)
     */
    int num_samples = num_affected + num_unaffected;
    uint8_t *masks = info.masks;
    
    assert(masks);
    
    for (int j = 0; j < order; j++) {
        
//         printf("snp %d = {\n", j);
        
        // Genotypes in the range (0,2)
        for (int i = 0; i < NUM_GENOTYPES; i++) {
            int k = 0;
//             printf("(1) ");
            for (; k < num_affected; k++) {
//                 printf("%d ", genotypes[j][k]);
                masks[j * NUM_GENOTYPES * (info.num_samples_per_mask) + i * (info.num_samples_per_mask) + k] = (genotypes[j][k] == i);
            }
//             printf("(2-");
            for (; k < info.num_affected_with_padding; k++) {
                masks[j * NUM_GENOTYPES * (info.num_samples_per_mask) + i * (info.num_samples_per_mask) + k] = 0;
            }
//             printf("3) ");
            for (k = 0; k < num_unaffected; k++) {
//                 printf("%d ", genotypes[j][num_affected + k]);
                masks[j * NUM_GENOTYPES * (info.num_samples_per_mask) + i * (info.num_samples_per_mask) + info.num_affected_with_padding + k] = (genotypes[j][num_affected + k] == i);
            }
//             printf("(4) ");
            for (; k < info.num_unaffected_with_padding; k++) {
                masks[j * NUM_GENOTYPES * (info.num_samples_per_mask) + i * (info.num_samples_per_mask) + info.num_affected_with_padding + k] = 0;
            }
//             printf("\n");
        }
    }
    
    return masks;
}

void masks_info_new(int order, int num_affected, int num_unaffected, masks_info *info) {
    info->num_affected_with_padding = 16 * (int) ceil(((double) num_affected) / 16);
    info->num_unaffected_with_padding = 16 * (int) ceil(((double) num_unaffected) / 16);
    info->num_samples_per_mask = info->num_affected_with_padding + info->num_unaffected_with_padding;
    info->num_masks = NUM_GENOTYPES * order * info->num_samples_per_mask;
    info->masks = _mm_malloc(info->num_masks * sizeof(uint8_t), 16);
    
    assert(info->masks);
    assert(info->num_affected_with_padding);
    assert(info->num_unaffected_with_padding);
}


/* **************************
 *         High risk        *
 * **************************/

int* get_high_risk_combinations(unsigned int* counts, unsigned int num_counts, unsigned int num_affected, unsigned int num_unaffected, 
                                unsigned int *num_risky, void** aux_ret, 
                                bool (*test_func)(unsigned int, unsigned int, unsigned int, unsigned int, void **)) {
    int *risky = malloc ((num_counts / 2) * sizeof(int));
    *num_risky = 0;
    
    for (int i = 0; i < num_counts; i += 2) {
        void *test_return_values = NULL;
        bool is_high_risk = test_func(counts[i], counts[i+1], num_affected, num_unaffected, &test_return_values);
        
        if (is_high_risk) {
            risky[*num_risky] = i / 2;
            if (test_return_values) { *aux_ret = test_return_values; }
            (*num_risky)++;
        }
    }
    
    return risky;
}

risky_combination* risky_combination_new(int order, int comb[order], uint8_t** possible_genotypes_combinations, int num_risky, int* risky_idx, void *aux_info) {
    risky_combination *risky = malloc(sizeof(risky_combination));
    risky->order = order;
    risky->combination = malloc(order * sizeof(int));
    // TODO for SSE, could it be useful to create N arrays, where N = order?
    risky->genotypes = malloc(num_risky * order * sizeof(uint8_t));
    risky->num_risky_genotypes = num_risky;
    risky->auxiliary_info = aux_info; // TODO improvement: set this using a method-dependant (MDR, MB-MDR) function
    
    memcpy(risky->combination, comb, order * sizeof(int));
    
    for (int i = 0; i < num_risky; i++) {
        memcpy(risky->genotypes + (order * i), possible_genotypes_combinations[risky_idx[i]], order * sizeof(uint8_t));
    }
    
    return risky;
}

void risky_combination_free(risky_combination* combination) {
    free(combination->combination);
    free(combination->genotypes);
    free(combination);
}


/* **************************
 *  Evaluation and ranking  *
 * **************************/

unsigned int *get_confusion_matrix(int order, risky_combination *combination, int num_affected_in_fold, int num_unaffected_in_fold, uint8_t *genotypes) {
    // TP, FN, FP, TN
    unsigned int *rates = calloc(4, sizeof(unsigned int));
    int num_samples = num_affected_in_fold + num_unaffected_in_fold;
    
    for (int i = 0; i < num_samples; i++) {
        bool marked_affected = 0;
        
        // Search through all the possible combinations until one of them applies
        for (int j = 0; j < combination->num_risky_genotypes && !marked_affected; j++) {
            marked_affected = 1;
            for (int k = 0; k < order && marked_affected; k++) {
                // If some of the genotypes in a combination does not match, don't keep checking it
//                 LOG_DEBUG_F("[%d,%d,%d] %d == %d\n", i, j, k, combination->genotypes[j * order + k], genotypes[k * num_samples + i]);
                marked_affected = combination->genotypes[j * order + k] == genotypes[k * num_samples + i];
//                 marked_affected = combination->genotypes[j * order + k] == genotypes[i][k];
            }
        }
        
//         LOG_DEBUG_F("marked affected? %d\n", marked_affected);
        
        if (marked_affected) {
            if (i < num_affected_in_fold) {
                (rates[0])++; // TP++
            } else {
                (rates[2])++; // FP++
            }
        } else {
            if (i < num_affected_in_fold) {
                (rates[1])++; // FN++
            } else {
                (rates[3])++; // TN++
            }
        }
    }
    
    assert(rates[0] + rates[1] + rates[2] + rates[3] == num_samples);
    
    return rates;
}


double evaluate_model(unsigned int *confusion_matrix, enum eval_function function) {
    double TP = confusion_matrix[0], FN = confusion_matrix[1], FP = confusion_matrix[2], TN = confusion_matrix[3];
    
    if (!function) {
        function = BA;
    }
    
    switch(function) {
        case CA:
            return (TP + TN) / (TP + FN + TN + FP);
        case BA:
            return ((TP / (TP + FN)) + (TN / (TN + FP))) / 2;
        case GAMMA:
            return (TP * TN - FP * FN) / (TP * TN + FP * FN);
        case TAU_B:
            return (TP * TN - FP * FN) / sqrt((TP + FN) * (TN + FP) * (TP + FP) * (TN + FN));
    }
}

