#ifndef EPISTASIS_DATASET_H
#define EPISTASIS_DATASET_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <commons/file_utils.h>
#include <containers/array_list.h>

#include "hpg_variant_utils.h"

/* ***************************
 *  Whole dataset management *
 * ***************************/

uint8_t *epistasis_dataset_load(int *num_affected, int *num_unaffected, size_t *num_variants, size_t *file_len, size_t *genotypes_offset, char *filename);

int epistasis_dataset_close(uint8_t *contents, size_t file_len);


/* *********************************************
 *  Combinations of blocks, SNPs and genotypes *
 * *********************************************/

int get_block_stride(size_t block_operations, int order);

int get_next_block(int num_blocks, int order, int block_coordinates[order]);

void get_first_combination_in_block(int order, int init_coordinates[order], int block_coordinates[order], int stride);

int get_next_combination_in_block(int order, int comb[order], int block_coordinates[order], int stride, int num_variants);

uint8_t **get_genotype_combinations(int order, int *num_combinations);

uint8_t get_next_genotype_combination(int order, uint8_t comb[order]);


/* ***************************
 *        Input/Output       *
 * ***************************/

void print_combination(int comb[], unsigned long idx, int k);

void print_gt_combination(uint8_t comb[], unsigned long idx, int k);

#endif