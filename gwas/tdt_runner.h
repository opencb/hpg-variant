#ifndef TDT_RUNNER_H
#define TDT_RUNNER_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cprops/hashtable.h>
#include <omp.h>

#include <bioformats/family.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/ped/ped_file_structure.h>
#include <bioformats/vcf/vcf_batch.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_util.h>
#include <commons/log.h>
#include <commons/string_utils.h>
#include <containers/list.h>

#include "global_options.h"
#include "gwas.h"
#include "hpg_variant_utils.h"
#include "tdt.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


int run_tdt_test(global_options_data_t *global_options_data, gwas_options_data_t *options_data);


static cp_hashtable *associate_samples_and_positions(vcf_file_t *file);

#endif
