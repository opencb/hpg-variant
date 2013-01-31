#ifndef EPISTASIS_MODEL
#define EPISTASIS_MODEL

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <commons/log.h>
#include <containers/array_list.h>

#include "dataset.h"


typedef struct {
    int *combination;
    uint8_t *genotypes;
    int num_risky;
} risky_combination;

/**
 * @brief Functions available for evaluating the results of a MDR model
 * @details Functions available for evaluating the results of a MDR model. These are:
 * - CA: accuracy
 * - BA: balanced accuracy, used by the original MDR
 * - wBA: weighted balanced accuracy, Namkung et al. (2008) TODO
 * - Gamma: Goodman and Kruskal (1954)
 * - Tau-b: Kendall rank correlation coefficient
 **/
enum eval_function { CA, BA, wBA, GAMMA, TAU_B };



/* **************************
 *          Counts          *
 * **************************/

/**
 * @brief Gets the number of ocurrences of each genotype both in affected and unaffected groups.
 * @details Gets the number of ocurrences of each genotype both in affected and unaffected groups. For using 
 * these values in a contingency table, number of not-occurrences can be calculated like the following:
 * not_occur_affected = num_affected - occur_affected
 *
 * @param order Number of SNPs combined
 * @param genotypes Genotypes of all SNPs (grouped by SNP)
 * @param genotype_combinations Possible genotype combinations for the given order
 * @param num_genotype_combinations Number of genotype combinations for the given order
 * @param num_affected Number of affected samples
 * @param num_unaffected Number of unaffected samples
 * @param num_counts Number of counts returned
 * @return List of counts, paired in (affected,unaffected)
 **/
int* get_counts(int order, uint8_t* genotypes, uint8_t **genotype_combinations, int num_genotype_combinations, int num_affected, int num_unaffected, int *num_counts);

uint8_t* get_masks(int order, uint8_t *genotypes, int num_samples, int *num_masks);


/* **************************
 *         High risk        *
 * **************************/

int* get_high_risk_combinations(unsigned int *counts, unsigned int num_counts, unsigned int num_affected, unsigned int num_unaffected, 
                                unsigned int *num_risky, array_list_t* aux_ret,
                                bool (*test_func)(unsigned int, unsigned int, unsigned int, unsigned int, void **));


risky_combination *risky_combination_new(int order, int comb[order], uint8_t **possible_genotypes_combinations, int num_risky, int *risky_idx);

void risky_combination_free(risky_combination *combination);



/* **************************
 *  Evaluation and ranking  *
 * **************************/

double evaluate_model(unsigned int true_positives, unsigned int true_negatives, unsigned int false_positives, unsigned int false_negatives, enum eval_function function);

#endif
