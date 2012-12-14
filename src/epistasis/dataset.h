#ifndef EPISTASIS_DATASET
#define EPISTASIS_DATASET

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <containers/array_list.h>

typedef struct {
    uint16_t num_affected;
    uint16_t num_unaffected;
    array_list_t *genotype_counts;
} epistasis_dataset;

// epistasis_dataset *epistasis_dataset_new();

// void epistasis_dataset_free(epistasis_dataset *dataset);


uint8_t *epistasis_dataset_process_records(vcf_record_t **variants, size_t num_variants, int num_samples);

void epistasis_dataset_add_entry(uint16_t *genotype_count, epistasis_dataset *dataset);


size_t epistasis_dataset_get_num_variants(epistasis_dataset *dataset);

uint8_t *epistasis_dataset_get_variant_counts(size_t index, epistasis_dataset *dataset);


#endif