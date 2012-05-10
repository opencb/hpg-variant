#ifndef TDT_RUNNER_H
#define TDT_RUNNER_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cprops/hashtable.h>
#include <omp.h>

#include <family.h>
#include <list.h>
#include <log.h>
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
#include "tdt.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


int run_tdt_test(global_options_data_t *global_options_data, gwas_options_data_t *options_data);


static cp_hashtable *associate_samples_and_positions(vcf_file_t *file);

#endif
