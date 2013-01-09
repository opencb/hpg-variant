#ifndef EPISTASIS_MODEL
#define EPISTASIS_MODEL

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <containers/array_list.h>

#include "dataset.h"

/**
 * @brief Gets the number of ocurrences of each genotype both in affected and unaffected groups.
 * @details Gets the number of ocurrences of each genotype both in affected and unaffected groups. For using 
 * these values in a contingency table, number of not-occurrences can be calculated like the following:
 * not_occur_affected = num_affected - occur_affected
 *
 * @param order Number of SNPs combined
 * @param genotypes Genotypes of all SNPs (grouped by SNP)
 * @param num_affected Number of affected samples
 * @param num_unaffected Number of unaffected samples
 * @param num_counts Number of counts returned
 * @return List of counts, paired in (affected,unaffected)
 **/
int* get_counts(int order, uint8_t *genotypes, int num_affected, int num_unaffected, int *num_counts);

uint8_t* get_masks(int order, uint8_t *genotypes, int num_samples, int *num_masks);


int* get_high_risk_combinations(int *counts, int num_counts, int num_affected, int num_unaffected, 
                                int *num_risky, array_list_t* aux_ret,
                                bool (*test_func)(int, int, int, int, void **));

#endif
