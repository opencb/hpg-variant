#ifndef EPISTASIS_DATASET
#define EPISTASIS_DATASET

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <containers/array_list.h>


uint8_t *epistasis_dataset_process_records(vcf_record_t** variants, size_t num_variants, int* destination, int num_samples);


int get_block_stride(size_t block_operations, int order);

int get_next_block(int num_blocks, int order, int block_coordinates[order]);

int* get_first_combination_in_block(int order, int block_coordinates[order], int stride);

int get_next_combination_in_block(int order, int comb[order], int block_coordinates[order], int stride);

int get_next_genotype_combination(int order, int comb[order]);


void print_combination(int comb[], unsigned long idx, int k);

void print_gt_combination(uint8_t comb[], unsigned long idx, int k);

#endif