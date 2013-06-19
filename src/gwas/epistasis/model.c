#include "model.h"


/* **************************
 *          Counts          *
 * **************************/

uint8_t* set_genotypes_masks(int order, uint8_t **genotypes, int num_combinations, masks_info info) {
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
    __m128i reference_genotype; // The genotype to compare for generating a mask (of the form {0 0 0 0 ... }, {1 1 1 1 ... })
    __m128i input_genotypes;    // Genotypes from the input dataset
    __m128i mask;               // Comparison between the reference genotype and input genotypes
    
    uint8_t *masks = info.masks;
    
    for (int c = 0; c < num_combinations; c++) {
        masks = info.masks + c * info.num_masks;
        uint8_t **combination_genotypes = genotypes + c * order;

        for (int j = 0; j < order; j++) {
            for (int i = 0; i < NUM_GENOTYPES; i++) {
                reference_genotype = _mm_set1_epi8(i);

                // Set value of masks
                for (int k = 0; k < info.num_samples_with_padding; k += 16) {
                    input_genotypes = _mm_load_si128(combination_genotypes[j] + k);
                    mask = _mm_cmpeq_epi8(input_genotypes, reference_genotype);
                    _mm_store_si128(masks + j * NUM_GENOTYPES * (info.num_samples_with_padding) + i * (info.num_samples_with_padding) + k, mask);
                }

                // Set padding with zeroes
                memset(masks + j * NUM_GENOTYPES * (info.num_samples_with_padding) + i * (info.num_samples_with_padding) + info.num_affected,
                    0, info.num_affected_with_padding - info.num_affected);
                memset(masks + j * NUM_GENOTYPES * (info.num_samples_with_padding) + i * (info.num_samples_with_padding) +
                    info.num_affected_with_padding + info.num_unaffected,
                    0, info.num_unaffected_with_padding - info.num_unaffected);
            }
        }
    }

    masks = info.masks;

    return masks;
}

void combination_counts(int order, uint8_t *masks, uint8_t **genotype_permutations, int num_genotype_permutations,
                        int *counts_aff, int *counts_unaff, masks_info info) {
    uint8_t *permutation;
    int count = 0;

    __m128i snp_and, snp_cmp;

    for (int rc = 0; rc < info.num_combinations_in_a_row; rc++) {
        uint8_t *rc_masks = info.masks + rc * order * NUM_GENOTYPES * info.num_samples_with_padding;
        for (int c = 0; c < num_genotype_permutations; c++) {
            permutation = genotype_permutations[c];
            // print_gt_combination(permutation, c, order);
            count = 0;

            for (int i = 0; i < info.num_affected; i += 16) {
                // Aligned loading
                snp_and = _mm_load_si128(rc_masks + permutation[0] * info.num_samples_with_padding + i);

                // Perform AND operation with all SNPs in the combination
                for (int j = 1; j < order; j++) {
                    snp_cmp = _mm_load_si128(rc_masks + j * NUM_GENOTYPES * info.num_samples_with_padding +
                                             permutation[j] * info.num_samples_with_padding + i);
                    snp_and = _mm_and_si128(snp_and, snp_cmp);
                }

                count += _mm_popcnt_u64(_mm_extract_epi64(snp_and, 0)) +
                         _mm_popcnt_u64(_mm_extract_epi64(snp_and, 1));
            }

            LOG_DEBUG_F("aff comb idx (%d) = %d\n", c, count / 8);
            counts_aff[rc * info.num_cell_counts_per_combination + c] = count / 8;

            count = 0;

            for (int i = 0; i < info.num_unaffected; i += 16) {
                // Aligned loading
                snp_and = _mm_load_si128(rc_masks + permutation[0] * info.num_samples_with_padding + info.num_affected_with_padding + i);

                // Perform AND operation with all SNPs in the combination
                for (int j = 1; j < order; j++) {
                    snp_cmp = _mm_load_si128(rc_masks + j * NUM_GENOTYPES * info.num_samples_with_padding +
                                             permutation[j] * info.num_samples_with_padding + info.num_affected_with_padding + i);
                    snp_and = _mm_and_si128(snp_and, snp_cmp);
                }

                count += _mm_popcnt_u64(_mm_extract_epi64(snp_and, 0)) +
                         _mm_popcnt_u64(_mm_extract_epi64(snp_and, 1));
            }

            LOG_DEBUG_F("unaff comb idx (%d) = %d\n", c, count / 8);
            counts_unaff[rc * info.num_cell_counts_per_combination + c] = count / 8;
        }
    }
}

void combination_counts_all_folds(int order, uint8_t *fold_masks, int num_folds,
                                  uint8_t **genotype_permutations, masks_info info, 
                                  int *counts_aff, int *counts_unaff) {
    uint8_t *permutation;
    int count[num_folds];

    __m128i snp_and, snp_cmp, snp_result;

    for (int rc = 0; rc < info.num_combinations_in_a_row; rc++) {
        uint8_t *rc_masks = info.masks + rc * order * NUM_GENOTYPES * info.num_samples_with_padding;
        for (int c = 0; c < info.num_cell_counts_per_combination; c++) {
            permutation = genotype_permutations[c];
            // print_gt_combination(permutation, c, order);
            
            memset(count, 0, num_folds * sizeof(int));

            for (int i = 0; i < info.num_affected; i += 16) {
                // Aligned loading
                snp_and = _mm_load_si128(rc_masks + permutation[0] * info.num_samples_with_padding + i);

                // Perform AND operation with all SNPs in the combination
                for (int j = 1; j < order; j++) {
                    snp_cmp = _mm_load_si128(rc_masks + j * NUM_GENOTYPES * info.num_samples_with_padding +
                                             permutation[j] * info.num_samples_with_padding + i);
                    snp_and = _mm_and_si128(snp_and, snp_cmp);
                }

                // Final AND with fold_masks
                for (int f = 0; f < num_folds; f++) {
                    snp_cmp = _mm_load_si128(fold_masks + f * info.num_samples_with_padding + i);
                    snp_result = _mm_and_si128(snp_and, snp_cmp);
                    
                    count[f] += _mm_popcnt_u64(_mm_extract_epi64(snp_result, 0)) +
                                _mm_popcnt_u64(_mm_extract_epi64(snp_result, 1));
                }
            }

            // Assign to count in fold
            for (int f = 0; f < num_folds; f++) {
                LOG_DEBUG_F("%d) aff comb idx (%d) = %d\n", f, c, count[f]);
                counts_aff[f * info.num_combinations_in_a_row * info.num_cell_counts_per_combination +
                           rc * info.num_cell_counts_per_combination + c] = count[f];
            }

            memset(count, 0, num_folds * sizeof(int));

            for (int i = 0; i < info.num_unaffected; i += 16) {
                // Aligned loading
                snp_and = _mm_load_si128(rc_masks + permutation[0] * info.num_samples_with_padding + info.num_affected_with_padding + i);

                // Perform AND operation with all SNPs in the combination
                for (int j = 1; j < order; j++) {
                    snp_cmp = _mm_load_si128(rc_masks + j * NUM_GENOTYPES * info.num_samples_with_padding +
                                             permutation[j] * info.num_samples_with_padding + info.num_affected_with_padding + i);
                    snp_and = _mm_and_si128(snp_and, snp_cmp);
                }

                // Final AND with fold_masks
                for (int f = 0; f < num_folds; f++) {
                    snp_cmp = _mm_load_si128(fold_masks + f * info.num_samples_with_padding + info.num_affected_with_padding + i);
                    snp_result = _mm_and_si128(snp_cmp, snp_and);
                    
                    count[f] += _mm_popcnt_u64(_mm_extract_epi64(snp_result, 0)) +
                                _mm_popcnt_u64(_mm_extract_epi64(snp_result, 1));
                }
            }

            // Assign to count in fold
            for (int f = 0; f < num_folds; f++) {
                LOG_DEBUG_F("%d) unaff comb idx (%d) = %d\n", f, c, count[f]);
                counts_unaff[f * info.num_combinations_in_a_row * info.num_cell_counts_per_combination +
                           rc * info.num_cell_counts_per_combination + c] = count[f];
            }
        }
    }
}

void masks_info_init(int order, int num_combinations_in_a_row, int num_affected, int num_unaffected, masks_info *info) {
    info->num_affected = num_affected;
    info->num_unaffected = num_unaffected;
    info->num_affected_with_padding = 16 * (int) ceil(((double) num_affected) / 16);
    info->num_unaffected_with_padding = 16 * (int) ceil(((double) num_unaffected) / 16);
    info->num_combinations_in_a_row = num_combinations_in_a_row;
    info->num_cell_counts_per_combination = pow(NUM_GENOTYPES, order);
    info->num_samples_with_padding = info->num_affected_with_padding + info->num_unaffected_with_padding;
    info->num_masks = NUM_GENOTYPES * order * info->num_samples_with_padding;
    info->masks = _mm_malloc(info->num_combinations_in_a_row * info->num_masks * sizeof(uint8_t), 16);

    assert(info->masks);
    assert(info->num_affected_with_padding);
    assert(info->num_unaffected_with_padding);
}


/* **************************
 *         High risk        *
 * **************************/

int* choose_high_risk_combinations2(unsigned int* counts_aff, unsigned int* counts_unaff, 
                                   unsigned int num_combinations, unsigned int num_counts_per_combination,
                                   unsigned int num_affected, unsigned int num_unaffected, 
                                   unsigned int *num_risky, void** aux_ret, 
                                   int* (*test_func)(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, void **)) {
    int num_counts = num_combinations * num_counts_per_combination;
    
    void *test_return_values = NULL;
    // Check high risk for all combinations
    int *is_high_risk = test_func(counts_aff, counts_unaff, num_counts, num_affected, num_unaffected, &test_return_values);
        
    int *risky = malloc (num_counts * sizeof(int)); // Put all risky indexes together
    
    int total_risky = 0;
    for (int i = 0; i < num_counts; i++) {
        if (is_high_risk[i]) {
            int c = i / num_counts_per_combination;
            int idx = i % num_counts_per_combination;
            
            risky[total_risky] = idx;
            num_risky[c]++;
            total_risky++;
        }
    }
    
    free(is_high_risk);
    
    return risky;
}

int* choose_high_risk_combinations(unsigned int* counts_aff, unsigned int* counts_unaff, unsigned int num_counts, 
                                   unsigned int num_affected, unsigned int num_unaffected, 
                                   unsigned int *num_risky, void** aux_ret, 
                                   bool (*test_func)(unsigned int, unsigned int, unsigned int, unsigned int, void **)) {
    int *risky = malloc (num_counts * sizeof(int));
    *num_risky = 0;
    
    for (int i = 0; i < num_counts; i++) {
        void *test_return_values = NULL;
        bool is_high_risk = test_func(counts_aff[i], counts_unaff[i], num_affected, num_unaffected, &test_return_values);
        
        if (is_high_risk) {
            risky[*num_risky] = i;
            if (test_return_values) { *aux_ret = test_return_values; }
            (*num_risky)++;
        }
    }
    
    return risky;
}

risky_combination* risky_combination_new(int order, int comb[order], uint8_t** possible_genotypes_combinations, 
                                         int num_risky, int* risky_idx, void *aux_info) {
    risky_combination *risky = malloc(sizeof(risky_combination));
    risky->order = order;
    risky->combination = malloc(order * sizeof(int));
    risky->accuracy = 0.0f;
    risky->genotypes = malloc(pow(NUM_GENOTYPES, order) * order * sizeof(uint8_t)); // Maximum possible
    risky->num_risky_genotypes = num_risky;
    risky->auxiliary_info = aux_info; // TODO improvement: set this using a method-dependant (MDR, MB-MDR) function
    
    memcpy(risky->combination, comb, order * sizeof(int));
    
    for (int i = 0; i < num_risky; i++) {
        memcpy(risky->genotypes + (order * i), possible_genotypes_combinations[risky_idx[i]], order * sizeof(uint8_t));
    }
    
    return risky;
}

risky_combination* risky_combination_copy(int order, int comb[order], uint8_t** possible_genotypes_combinations, 
                                          int num_risky, int* risky_idx, void *aux_info, risky_combination* risky) {
    assert(risky);
    risky->num_risky_genotypes = num_risky;
    risky->auxiliary_info = aux_info; // TODO improvement: set this using a method-dependant (MDR, MB-MDR) function
    risky->accuracy = 0.0f;
    
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

double test_model(int order, risky_combination *risky_comb, uint8_t **genotypes, masks_info info, unsigned int *conf_matrix) {
    // Get the matrix containing {FP,FN,TP,TN}
    confusion_matrix(order, risky_comb, info, genotypes, conf_matrix);

    // Evaluate the model, basing on the confusion matrix
    double eval = evaluate_model(conf_matrix, BA);
    risky_comb->accuracy = eval;

    return eval;
}

void confusion_matrix(int order, risky_combination *combination, masks_info info, uint8_t **genotypes, unsigned int *matrix) {
    int num_samples = info.num_samples_with_padding;
    uint8_t confusion_masks[combination->num_risky_genotypes * num_samples];
    memset(confusion_masks, 0, combination->num_risky_genotypes * num_samples * sizeof(uint8_t));
    
    __m128i comb_genotypes;     // The genotype to compare for generating a mask (of the form {0 0 0 0 ... }, {1 1 1 1 ... })
    __m128i input_genotypes;    // Genotypes from the input dataset
    __m128i mask;               // Comparison between the reference genotype and input genotypes
    
    
    // Check whether the input genotypes can be combined in any of the risky combinations
    for (int i = 0; i < combination->num_risky_genotypes; i++) {
        // First SNP in the combination
        comb_genotypes = _mm_set1_epi8(combination->genotypes[i * order]);
        
        for (int k = 0; k < info.num_samples_with_padding; k += 16) {
            input_genotypes = _mm_load_si128(genotypes[0] + k);
            mask = _mm_cmpeq_epi8(input_genotypes, comb_genotypes);
            _mm_store_si128(confusion_masks + i * num_samples + k, mask);
        }
        
        // Next SNPs in the combination
        for (int j = 1; j < order; j++) {
            comb_genotypes = _mm_set1_epi8(combination->genotypes[i * order + j]);
            
            for (int k = 0; k < info.num_samples_with_padding; k += 16) {
                input_genotypes = _mm_load_si128(genotypes[j] + k);
                mask = _mm_load_si128(confusion_masks + i * num_samples + k);
                mask = _mm_and_si128(mask, _mm_cmpeq_epi8(input_genotypes, comb_genotypes));
                _mm_store_si128(confusion_masks + i * num_samples + k, mask);
            }
        }
    }
    
/*
    printf("confusion masks sse = {\n");
    for (int j = 0; j < combination->num_risky_genotypes; j++) {
        printf(" comb %d = { ", j);
        for (int k = 0; k < num_samples; k++) {
            printf("%03d ", confusion_masks[j * num_samples + k]);
        }
        printf("}\n");
    }
    printf("}\n");
*/
   
    uint8_t final_masks[num_samples];
    __m128i final_or, other_mask;
    
    for (int k = 0; k < num_samples; k += 16) {
        final_or = _mm_load_si128(confusion_masks + k); // TODO first mask
        
        // Merge all positives (1) and negatives (0)
        for (int j = 1; j < combination->num_risky_genotypes; j++) {
            other_mask = _mm_load_si128(confusion_masks + j * num_samples + k);
            final_or = _mm_or_si128(final_or, other_mask);
        }
        
        _mm_store_si128(final_masks + k, final_or);
    }
    
/*
    printf("final masks sse = {\n");
    for (int k = 0; k < num_samples; k++) {
        printf("%d ", final_masks[k]);
    }
    printf("}\n");
*/
   
    // Get the counts (popcount is the number of 1s -> popcount / 8 is the number of positives)
    int popcount0 = 0, popcount1 = 0;
    __m128i snp_and;
    
    memset(final_masks + info.num_affected, 0, info.num_affected_with_padding - info.num_affected);
    memset(final_masks + info.num_affected_with_padding + info.num_unaffected, 0, info.num_unaffected_with_padding - info.num_unaffected);
    
    for (int k = 0; k < info.num_affected; k += 16) {
        snp_and = _mm_load_si128(final_masks + k);
        popcount0 += _mm_popcnt_u64(_mm_extract_epi64(snp_and, 0)) + 
                     _mm_popcnt_u64(_mm_extract_epi64(snp_and, 1));
    }
    
    for (int k = 0; k < info.num_unaffected; k += 16) {
        snp_and = _mm_load_si128(final_masks + info.num_affected_with_padding + k);
        popcount1 += _mm_popcnt_u64(_mm_extract_epi64(snp_and, 0)) + 
                     _mm_popcnt_u64(_mm_extract_epi64(snp_and, 1));
    }
    
    matrix[0] = popcount0 / 8;
    matrix[1] = info.num_affected - popcount0 / 8;
    matrix[2] = popcount1 / 8;
    matrix[3] = info.num_unaffected - popcount1 / 8;
    
/*
    assert(matrix[0] + matrix[1] + matrix[2] + matrix[3] == info.num_affected + info.num_unaffected);
*/
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

int add_to_model_ranking(risky_combination *risky_comb, int max_ranking_size, struct heap *ranking_risky,
                         int (*priority_func) (struct heap_node* a, struct heap_node* b)) {
    // Step 6 -> Construct ranking of the best N combinations
    size_t current_ranking_size = ranking_risky->size;

    if (current_ranking_size > 0) {
        struct heap_node *last_node = heap_peek(priority_func, ranking_risky);
        risky_combination *last_element = last_node->value;

        // If accuracy is not greater than the last element, don't bother inserting
        if (risky_comb->accuracy > last_element->accuracy) {
            struct heap_node *hn = malloc (sizeof(struct heap_node));
            heap_node_init(hn, risky_comb);
            heap_insert(priority_func, ranking_risky, hn);

            if (current_ranking_size >= max_ranking_size) {
                struct heap_node *removed = heap_take(priority_func, ranking_risky);
                risky_combination_free((risky_combination*) removed->value);
                free(removed);
            }

            return ranking_risky->size - 1;
        }

        if (current_ranking_size < max_ranking_size) {
            LOG_DEBUG_F("To insert %.3f at the end", risky_comb->accuracy);
            struct heap_node *hn = malloc (sizeof(struct heap_node));
            heap_node_init(hn, risky_comb);
            heap_insert(priority_func, ranking_risky, hn);
            return ranking_risky->size - 1;
        }
    } else {
        struct heap_node *hn = malloc (sizeof(struct heap_node));
        heap_node_init(hn, risky_comb);
        heap_insert(priority_func, ranking_risky, hn);
        return ranking_risky->size - 1;

    }

    return -1;
}


int compare_risky_heap_max(struct heap_node* a, struct heap_node* b) {
    risky_combination *r1 = (risky_combination*) a->value;
    risky_combination *r2 = (risky_combination*) b->value;
    return r1->accuracy > r2->accuracy;
}

int compare_risky_heap_min(struct heap_node* a, struct heap_node* b) {
    risky_combination *r1 = (risky_combination*) a->value;
    risky_combination *r2 = (risky_combination*) b->value;
    return r1->accuracy < r2->accuracy;
}
