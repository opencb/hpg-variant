#include "model.h"

#define NUM_GENOTYPES   3

int* get_counts(int order, uint8_t* genotypes, int num_affected, int num_unaffected, int *num_counts) {
    int num_masks;
    int num_samples = num_affected + num_unaffected;
    *num_counts = 2 * pow(NUM_GENOTYPES, order);
    uint8_t *masks = get_masks(order, genotypes, num_samples, &num_masks); // Grouped by SNP
    int *counts = malloc((*num_counts) * sizeof(int)); // Affected and unaffected
    
    int comb[order]; memset(comb, 0, order * sizeof(int));
    int comb_idx = 0, flag = 1, count = 0;
    do {
//         print_combination(comb, comb_idx, order);
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
//         printf("aff comb idx (%d) = %d\n", comb_idx * 2, count);
        counts[comb_idx * 2] = count;
        
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
//         printf("unaff comb idx (%d) = %d\n", comb_idx * 2 + 1, count);
        counts[comb_idx * 2 + 1] = count;
        
        comb_idx++;
    } while (get_next_genotype_combination(order, comb));
    // TODO Optimization: precalculate combinations for a given order, replace {order,genotypes} by {combinations}
    
    return counts;
}

uint8_t* get_masks(int order, uint8_t *genotypes, int num_samples, int *num_masks) {
    /* 
     * Structure: Genotypes of a SNP in each 'row'
     * 
     * SNP(0) - Mask genotype 0
     * SNP(0) - Mask genotype 1
     * SNP(0) - Mask genotype 2
     * 
     * SNP(1) - Mask genotype 0
     * SNP(1) - Mask genotype 1
     * SNP(1) - Mask genotype 2
     * 
     * ...
     * 
     * SNP(order-1) - Mask genotype 0
     * SNP(order-1) - Mask genotype 1
     * SNP(order-1) - Mask genotype 2
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


int* get_high_risk_combinations(int* counts, int num_counts, int num_affected, int num_unaffected, 
                                int *num_risky, array_list_t* aux_ret, 
                                bool (*test_func)(int, int, int, int, void **)) {
    int *combinations = malloc ((num_counts / 2) * sizeof(int));
    *num_risky = 0;
    
    for (int i = 0; i <= num_counts; i += 2) {
        void *test_return_values;
        bool is_high_risk = test_func(counts[i], counts[i+1], num_affected, num_unaffected, &test_return_values);
        if (is_high_risk) {
            combinations[*num_risky] = i / 2;
            if (test_return_values) {
                array_list_insert(test_return_values, aux_ret);
            }
            (*num_risky)++;
        }
    }
    
    return combinations;
}
