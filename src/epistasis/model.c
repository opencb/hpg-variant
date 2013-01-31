#include "model.h"

#define NUM_GENOTYPES   3


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

int *get_confusion_matrix(int order, risky_combination *combination, int num_affected_in_fold, int num_unaffected_in_fold, uint8_t *genotypes) {
    // TP, FN, FP, TN
    int *rates = calloc(4, sizeof(int));
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


// TODO change arguments order to: TP, FN, FP, TN
double evaluate_model(unsigned int true_positives, unsigned int true_negatives, unsigned int false_positives, unsigned int false_negatives, enum eval_function function) {
    double TP = true_positives, TN = true_negatives, FP = false_positives, FN = false_negatives;
    
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

