#ifndef EPISTASIS_MODEL
#define EPISTASIS_MODEL

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <popcntintrin.h>

#include <commons/log.h>
#include <containers/array_list.h>
#include <containers/heap.h>
#include <containers/linked_list.h>

#include "dataset.h"
#include "mdr.h"

#define NUM_GENOTYPES           3

#define COMBINATIONS_ROW_SSE    16
#define COMBINATIONS_ROW_GPU    128

typedef int (*compare_risky_heap_func) (struct heap_node*, struct heap_node*);

typedef struct {
    double accuracy;
    int order;
    int num_risky_genotypes;
    int cross_validation_count;
    uint8_t *genotypes;
    int *combination;
    void *auxiliary_info;
} risky_combination;


typedef struct {
    int num_affected;
    int num_unaffected;
    int num_affected_with_padding;
    int num_unaffected_with_padding;
    int num_samples_with_padding;
    int num_masks;
    int num_combinations_in_a_row;
    int num_cell_counts_per_combination;
    uint8_t *masks;
} masks_info;

enum evaluation_subset { TESTING, TRAINING };
enum evaluation_mode { CV_C, CV_A };

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

void set_genotypes_masks(int order, uint8_t **genotypes, int num_combinations, uint8_t *masks, masks_info info);

/**
 * @brief Gets the number of ocurrences of each genotype both in affected and unaffected groups.
 * @details Gets the number of ocurrences of each genotype both in affected and unaffected groups. For using 
 * these values in a contingency table, number of not-occurrences can be calculated like the following:
 * not_occur_affected = num_affected - occur_affected
 *
 * @param order Number of SNPs combined
 * @param genotype_combinations Possible genotype combinations for the given order
 * @param num_genotype_combinations Number of genotype combinations for the given order
 * @param num_counts Number of counts that will be calculated
 * @return List of counts, paired in (affected,unaffected)
 **/
void combination_counts(int order, uint8_t *masks, uint8_t **genotype_combinations, int num_genotype_combinations, 
                        int *counts_aff, int *counts_unaff, masks_info info);

void combination_counts_all_folds(int order, uint8_t *fold_masks, int num_folds,
                                  uint8_t **genotype_permutations, uint8_t *masks, masks_info info, 
                                  int *counts_aff, int *counts_unaff);

void masks_info_init(int order, int num_combinations_in_a_row, int num_affected, int num_unaffected, masks_info *info);


/* **************************
 *         High risk        *
 * **************************/

int* choose_high_risk_combinations2(unsigned int* counts_aff, unsigned int* counts_unaff, 
                                   unsigned int num_combinations, unsigned int num_counts_per_combination,
                                   unsigned int num_affected, unsigned int num_unaffected, 
                                   unsigned int *num_risky, void** aux_ret, 
                                   int* (*test_func)(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, void **));

int* choose_high_risk_combinations(unsigned int* counts_aff, unsigned int* counts_unaff, unsigned int num_counts, 
                                   unsigned int num_affected, unsigned int num_unaffected, 
                                   unsigned int *num_risky, void** aux_ret, 
                                   bool (*test_func)(unsigned int, unsigned int, unsigned int, unsigned int, void **));


risky_combination *risky_combination_new(int order, int comb[order], uint8_t **possible_genotypes_combinations, 
                                         int num_risky, int *risky_idx, void *aux_info, masks_info info);

risky_combination* risky_combination_copy(int order, int comb[order], uint8_t** possible_genotypes_combinations, 
                                          int num_risky, int* risky_idx, void *aux_info, risky_combination* risky);

void risky_combination_free(risky_combination *combination);


/* **************************
 *  Evaluation and ranking  *
 * **************************/

double test_model(int order, risky_combination *risky_comb, uint8_t **genotypes, 
                  uint8_t *fold_masks, enum evaluation_subset mode, int training_size[2], int testing_size[2], 
                  masks_info info, unsigned int *conf_matrix);

void confusion_matrix(int order, risky_combination *combination, uint8_t **genotypes, 
                      uint8_t *fold_masks, enum evaluation_subset mode, int training_size[2], int testing_size[2], 
                      masks_info info, unsigned int *matrix);

double evaluate_model(unsigned int *confusion_matrix, enum eval_function function);

int add_to_model_ranking(risky_combination *risky_comb, int max_ranking_size, struct heap *ranking_risky,
                         compare_risky_heap_func priority_func);


int compare_risky_heap_count_max(struct heap_node* a, struct heap_node* b);

int compare_risky_heap_count_min(struct heap_node* a, struct heap_node* b);

int compare_risky_heap_accuracy_max(struct heap_node* a, struct heap_node* b);

int compare_risky_heap_accuracy_min(struct heap_node* a, struct heap_node* b);

#endif
