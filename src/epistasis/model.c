#include "model.h"


int* get_counts(int order, uint8_t* gts, int num_affected, int num_unaffected) {
    
}

uint8_t* get_masks(int order, uint8_t *genotypes, int num_samples, int *num_masks) {
    /*
     * Structure: Genotypes of a SNP in each 'row'
     * 
     * SNP(0) - Mask genotype 0
     * SNP(1) - Mask genotype 0
     * ...
     * SNP(order-1) - Mask genotype 0
     * 
     * SNP(0) - Mask genotype 1
     * SNP(1) - Mask genotype 1
     * ...
     * SNP(order-1) - Mask genotype 1
     * 
     * SNP(0) - Mask genotype 2
     * SNP(1) - Mask genotype 2
     * ...
     * SNP(order-1) - Mask genotype 2
     */
    *num_masks = 3 * order * num_samples;
    uint8_t *masks = malloc(*num_masks * sizeof(uint8_t));
    int samples_per_order = order * num_samples;
    
    // Genotypes in the range (0,2)
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < order; j++) {
            for (int k = 0; k < num_samples; k++) {
                // samples_per_order allows to group by mask (0,1,2)
                // num_samples allows to get the row inside a group
                masks[i * samples_per_order + j * num_samples + k] = genotypes[j * num_samples + k] == i;
            }
        }
    }
    
    return masks;
}
