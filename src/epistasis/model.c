#include "model.h"

#define NUM_GENOTYPES   3


/* **************************
 *       Main pipeline      *
 * **************************/

risky_combination *get_model_from_combination_in_fold(int order, int comb[order], int *fold_samples, unsigned int fold_size, 
                                                      unsigned int num_samples, unsigned int num_affected_in_training, unsigned int num_unaffected_in_training,
                                                      int stride, uint8_t *block_starts[2],
                                                      int num_genotype_combinations, uint8_t **genotype_combinations,
                                                      array_list_t* aux_ret) {
    risky_combination *risky_comb = NULL;
    
    // Get genotypes of that combination
    uint8_t *val = get_genotypes_for_combination_exclude_fold(order, comb, num_samples, fold_size, fold_samples, stride, block_starts);
    
//     printf("* val = { ");
//     for (int i = 0; i < order; i++) {
//         for (int j = 0; j < (num_samples - num_samples_in_fold); j++) {
//             printf("%u ", val[i * (num_samples - num_samples_in_fold) + j]);
//         }
//         printf("\t");
//     }
//     printf("}\n\n");
    
    // Get counts for those genotypes
    int num_counts, num_risky;
    int *counts = get_counts(order, val, genotype_combinations, num_genotype_combinations, num_affected_in_training, num_unaffected_in_training, &num_counts);
    
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
    int *risky_idx = get_high_risk_combinations(counts, num_counts, num_affected_in_training, num_unaffected_in_training, 
                                                &num_risky, aux_ret, mdr_high_risk_combinations);
    
    // Filter non-risky SNP combinations
    if (num_risky > 0) {
        // Put together the info about the SNP combination and its genotype combinations
        risky_comb = risky_combination_new(order, comb, genotype_combinations, num_risky, risky_idx);
        
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
    
    free(val);
    free(counts);
    free(risky_idx);
    
    return risky_comb;
}


double test_model(int order, risky_combination *risky_comb, int *fold_samples, 
                  unsigned int num_samples, unsigned int num_affected_in_testing, unsigned int num_unaffected_in_testing,
                  int stride, uint8_t *block_starts[2]) {
    // Step 5 -> Check against the testing dataset
    uint8_t *val = get_genotypes_for_combination_and_fold(order, risky_comb->combination, 
                                                          num_samples, num_affected_in_testing + num_unaffected_in_testing, 
                                                          fold_samples, stride, block_starts);
    // Get the matrix containing {FP,FN,TP,TN}
    unsigned int *confusion_matrix = get_confusion_matrix(order, risky_comb, num_affected_in_testing, num_affected_in_testing, val);
    
//     printf("confusion matrix = { ");
//     for (int k = 0; k < 4; k++) {
//         printf("%u ", confusion_matrix[k]);
//     }
//     printf("}\n");
    
    // Evaluate the model, basing on the confusion matrix
    double eval = evaluate_model(confusion_matrix, BA);
    
    printf("risky combination = {\n  SNP: ");
    print_combination(risky_comb->combination, 0, order);
    printf("  GT: ");
    for (int j = 0; j < risky_comb->num_risky * 2; j++) {
        if (j % 2) {
            printf("%d), ", risky_comb->genotypes[j]);
        } else {
            printf("(%d ", risky_comb->genotypes[j]);
        }
    }
    printf("\n}\n", eval);
    
    risky_comb->accuracy = eval;
    
    return eval;
}


int add_to_model_ranking(risky_combination *risky_comb, int max_ranking_size, linked_list_t *ranking_risky) {
    // Step 6 -> Ellaborate a ranking of the best N combinations
    risky_combination *last_element = linked_list_get_last(ranking_risky);
    size_t current_ranking_size = ranking_risky->size;
    printf("Ranking (size %zu) = { ", current_ranking_size);
    
    linked_list_iterator_t* iter = linked_list_iterator_new(ranking_risky);
    risky_combination *element = NULL;
    while ((element = linked_list_iterator_curr(iter)) != last_element) {
        printf("(%d %d - %.3f) ", element->combination[0], element->combination[1], element->accuracy);
        element = linked_list_iterator_next(iter);
    }
    printf("}\n");
    linked_list_iterator_first(iter);
    
    if (current_ranking_size > 0) {
        bool inserted = false;
        
        LOG_DEBUG_F("To insert %.3f\tRanking's last is %.3f\n", risky_comb->accuracy, element->accuracy);
        
        // If accuracy is not greater than the last element, don't bother inserting
        if (risky_comb->accuracy > last_element->accuracy) {
            int position = 0;
            while ((element = linked_list_iterator_curr(iter)) != last_element) {
                LOG_DEBUG_F("To insert %.3f\tIn ranking (pos #%d) %.3f\n", risky_comb->accuracy, position, element->accuracy);
                if (risky_comb->accuracy > element->accuracy) {
                    linked_list_iterator_insert(risky_comb, iter);
                    linked_list_remove_last(ranking_risky);
                    inserted = true;
                    return position;
                }
                element = linked_list_iterator_next(iter);
                position++;
            }
        }
        
        if (!inserted && current_ranking_size < max_ranking_size) {
            linked_list_insert_last(risky_comb, ranking_risky);
            return ranking_risky->size - 1;
        }
    } else {
        linked_list_insert_last(risky_comb, ranking_risky);
        return ranking_risky->size - 1;
    }
    
    return -1;
}


/* **************************
 *          Counts          *
 * **************************/

int* get_counts(int order, uint8_t* genotypes, uint8_t **genotype_combinations, int num_genotype_combinations, int num_affected, int num_unaffected, int *num_counts) {
    int num_masks;
    int num_samples = num_affected + num_unaffected;
    *num_counts = 2 * pow(NUM_GENOTYPES, order);
    uint8_t *masks = get_masks(order, genotypes, num_samples, &num_masks); // Grouped by SNP
    int *counts = malloc((*num_counts) * sizeof(int)); // Affected and unaffected
    
    uint8_t *comb;
    int flag = 1, count = 0;
    
    for (int i = 0; i < num_genotype_combinations; i++) {
        comb = genotype_combinations[i];
//         print_combination(comb, i, order);
        
        flag = 1, count = 0;
        for (int i = 0; i < num_affected; i++) {
            for (int j = 0; j < order && flag; j++) {
                flag &= masks[j * NUM_GENOTYPES * num_samples + comb[j] * num_samples + i];
            }
            if (flag) {
                count++;
            }
            flag = 1;
        }
        LOG_DEBUG_F("aff comb idx (%d) = %d\n", i * 2, count);
        counts[i * 2] = count;
        
        flag = 1, count = 0;
        for (int i = num_affected; i < num_samples; i++) {
            for (int j = 0; j < order && flag; j++) {
                flag &= masks[j * NUM_GENOTYPES * num_samples + comb[j] * num_samples + i];
            }
            if (flag) {
                count++;
            }
            flag = 1;
        }
        LOG_DEBUG_F("unaff comb idx (%d) = %d\n", i * 2 + 1, count);
        counts[i * 2 + 1] = count;
    }
    
    return counts;
}

uint8_t* get_masks(int order, uint8_t *genotypes, int num_samples, int *num_masks) {
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
    *num_masks = NUM_GENOTYPES * order * num_samples;
    uint8_t *masks = malloc((*num_masks) * sizeof(uint8_t));
    
    for (int j = 0; j < order; j++) {
        // Genotypes in the range (0,2)
        for (int i = 0; i < NUM_GENOTYPES; i++) {
            for (int k = 0; k < num_samples; k++) {
                // group by SNP (better spatial locality than grouping by genotype (0/1/2))
                // num_samples allows to get the row inside a group
                masks[j * NUM_GENOTYPES * num_samples + i * num_samples + k] = (genotypes[j * num_samples + k] == i);
            }
        }
    }
    
    return masks;
}


/* **************************
 *         High risk        *
 * **************************/

int* get_high_risk_combinations(unsigned int* counts, unsigned int num_counts, unsigned int num_affected, unsigned int num_unaffected, 
                                unsigned int *num_risky, array_list_t* aux_ret, 
                                bool (*test_func)(unsigned int, unsigned int, unsigned int, unsigned int, void **)) {
    int *risky = malloc ((num_counts / 2) * sizeof(int));
    *num_risky = 0;
    
    for (int i = 0; i < num_counts; i += 2) {
        void *test_return_values = NULL;
        bool is_high_risk = test_func(counts[i], counts[i+1], num_affected, num_unaffected, &test_return_values);
        
        if (is_high_risk) {
            risky[*num_risky] = i / 2;
            if (test_return_values) {
                array_list_insert(test_return_values, aux_ret);
            }
            (*num_risky)++;
        }
    }
    
    return risky;
}

risky_combination* risky_combination_new(int order, int comb[order], uint8_t** possible_genotypes_combinations, int num_risky, int* risky_idx) {
    risky_combination *risky = malloc(sizeof(risky_combination));
    risky->combination = malloc(order * sizeof(int));
    // TODO for SSE, could it be useful to create N arrays, where N = order?
    risky->genotypes = malloc(num_risky * order * sizeof(uint8_t));
    risky->num_risky = num_risky;
    
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
        for (int j = 0; j < combination->num_risky && !marked_affected; j++) {
            marked_affected = 1;
            for (int k = 0; k < order && marked_affected; k++) {
                // If some of the genotypes in a combination does not match, don't keep checking it
                LOG_DEBUG_F("[%d,%d,%d] %d == %d\n", i, j, k, combination->genotypes[j * order + k], genotypes[k * num_samples + i]);
                marked_affected = combination->genotypes[j * order + k] == genotypes[k * num_samples + i];
            }
        }
        
        LOG_DEBUG_F("marked affected? %d\n", marked_affected);
        
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

