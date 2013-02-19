#ifndef EPISTASIS_DATASET_CREATOR
#define EPISTASIS_DATASET_CREATOR

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <omp.h>

#include <bioformats/family/family.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/ped/ped_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <commons/file_utils.h>
#include <containers/list.h>

#include "hpg_variant_utils.h"
#include "shared_options.h"

int create_dataset_from_vcf(shared_options_data_t* shared_options_data);


uint8_t *epistasis_dataset_process_records(vcf_record_t** variants, size_t num_variants, int* destination, 
                                           int num_samples, int num_threads);



static uint8_t *get_individual_phenotypes(vcf_file_t* vcf, ped_file_t* ped, int* num_affected, int* num_unaffected);

static int *group_individuals_by_phenotype(uint8_t *phenotypes, int num_affected, int num_unaffected);

static individual_t **sort_individuals(vcf_file_t *vcf, ped_file_t *ped);

static cp_hashtable* associate_samples_and_positions(vcf_file_t* file);

#endif
