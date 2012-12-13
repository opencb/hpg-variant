#ifndef EPISTASIS_DATASET
#define EPISTASIS_DATASET

#include <stdint.h>
#include <stdlib.h>
// #include <sys/types.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <containers/array_list.h>

typedef struct {
    uint16_t num_affected;
    uint16_t num_unaffected;
    array_list_t *genotype_counts;
} epistasis_dataset;

epistasis_dataset *epistasis_dataset_new();

void epistasis_dataset_free(epistasis_dataset *dataset);


void epistasis_dataset_process_records(vcf_record_t **variants, int num_variants, int num_samples, int *phenotypes, epistasis_dataset *dataset);

void epistasis_dataset_add_entry(uint16_t *genotype_count, epistasis_dataset *dataset);


uint16_t epistasis_dataset_get_num_variants(epistasis_dataset *dataset);

uint16_t *epistasis_dataset_get_variant_counts(size_t index, epistasis_dataset *dataset);


#endif