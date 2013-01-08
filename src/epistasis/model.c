#include "model.h"


int* get_counts(int order, uint8_t* genotypes, int num_affected, int num_unaffected) {
    int num_masks;
    int num_samples = num_affected + num_unaffected;
    uint8_t *masks = get_masks(order, genotypes, num_samples, &num_masks);
    int *counts = malloc(pow(3, order) * sizeof(int));
    
    // TODO
    
    
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
    *num_masks = 3 * order * num_samples;
    uint8_t *masks = malloc(*num_masks * sizeof(uint8_t));
    
    for (int j = 0; j < order; j++) {
        // Genotypes in the range (0,2)
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < num_samples; k++) {
                // group by SNP (better spatial locality than grouping by genotype (0/1/2))
                // num_samples allows to get the row inside a group
                masks[j * 3 * num_samples + i * num_samples + k] = (genotypes[j * num_samples + k] == i);
            }
        }
    }
    
    return masks;
}
