#ifndef TRANSMISSION_DISEQUILIBRIUM_TEST_H
#define TRANSMISSION_DISEQUILIBRIUM_TEST_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cprops/hashtable.h>
#include <omp.h>

#include <family.h>
#include <ped_file.h>
#include <ped_file_structure.h>
#include <string_utils.h>
#include <vcf_batch.h>
#include <vcf_file.h>
#include <vcf_filters.h>
#include <vcf_util.h>

#include "global_options.h"
#include "gwas.h"
#include "hpg_variant_utils.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

typedef struct {
    // TODO define fields
} tdt_result_t;


int run_tdt_test(ped_file_t *ped_file, global_options_data_t *global_options_data, gwas_options_data_t *options_data);

// tdt_result_t *
int
tdt_test(ped_file_t *ped_file, list_item_t *variants, int num_variants, cp_hashtable *sample_ids);

#endif
