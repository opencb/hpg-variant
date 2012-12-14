#ifndef EPISTASIS_DATASET_CREATOR
#define EPISTASIS_DATASET_CREATOR

#include <stdbool.h>
#include <stdlib.h>
#include <omp.h>

#include <bioformats/family/family.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/ped/ped_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <commons/file_utils.h>
#include <containers/list.h>

#include "dataset.h"
#include "hpg_variant_utils.h"
#include "shared_options.h"

int create_dataset_from_vcf(shared_options_data_t* shared_options_data);



static int flatten_phenotypes(vcf_file_t *vcf, ped_file_t *ped, int *num_affected, int *num_unaffected);

static uint8_t *get_individual_phenotypes(vcf_file_t *vcf, ped_file_t *ped);

static individual_t **sort_individuals(vcf_file_t *vcf, ped_file_t *ped);

static cp_hashtable* associate_samples_and_positions(vcf_file_t* file);

#endif
