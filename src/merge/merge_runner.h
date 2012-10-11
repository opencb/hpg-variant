#ifndef MERGE_RUNNER_H
#define MERGE_RUNNER_H

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <cprops/hashtable.h>
#include <cprops/linked_list.h>
#include <khash.h>
#include <omp.h>

#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/array_list.h>
#include <containers/list.h>
#include <bioformats/features/region/region.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_reader.h>

#include "hpg_vcf_tools_utils.h"
#include "merge.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

int run_merge(shared_options_data_t *shared_options_data, merge_options_data_t *options_data);

static void compose_key_value(const char *chromosome, const long position, char *key);

#endif
