#ifndef EPISTASIS_MDR
#define EPISTASIS_MDR

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include <commons/log.h>
#include <data/array_utils.h>

#include "hpg_variant_utils.h"


bool mdr_high_risk_combinations(unsigned int count_affected, unsigned int count_unaffected, 
                                unsigned int samples_affected, unsigned int samples_unaffected, void **aux_return_values);

int** get_k_folds(unsigned int samples_affected, unsigned int samples_unaffected, unsigned int k, unsigned int **sizes);

int *get_genotypes_for_combination_and_fold(int order, int comb[order], int num_samples, int num_samples_in_fold, int fold_samples[num_samples_in_fold], 
                                            int block_coordinates[order], int stride, uint8_t **block_starts);

#endif
