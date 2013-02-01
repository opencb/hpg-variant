#ifndef CROSS_VALIDATION_H
#define CROSS_VALIDATION_H

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include <commons/log.h>
#include <data/array_utils.h>

#include "hpg_variant_utils.h"


int** get_k_folds(unsigned int samples_affected, unsigned int samples_unaffected, unsigned int k, unsigned int **sizes);


uint8_t *get_genotypes_for_combination_and_fold(int order, int comb[order], int num_samples, 
                                                int num_samples_in_fold, int fold_samples[num_samples_in_fold], 
                                                int stride, uint8_t **block_starts);

/**
 * @brief Retrieves the genotypes of all samples except for the ones in the specified fold_samples
 *
 * @param order The number of SNPs combined at the same
 * @param comb The combination of SNPs
 * @param num_samples The total number of samples
 * @param num_samples_in_fold The number of samples in the fold to be excluded
 * @param fold_samples The samples in the fold to be excluded
 * @param stride The number of combinations in a block
 * @param block_starts Pointers to the first genotypes in the block
 * @return The list of genotypes of k-1 folds
 **/
uint8_t *get_genotypes_for_combination_exclude_fold(int order, int comb[order], int num_samples, 
                                                    int num_samples_in_fold, int fold_samples[num_samples_in_fold], 
                                                    int stride, uint8_t **block_starts);


#endif
