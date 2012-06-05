#ifndef SPLIT_RUNNER_H
#define SPLIT_RUNNER_H

#include <stdlib.h>

#include <cprops/hashtable.h>
#include <omp.h>

#include <commons/file_utils.h>
#include <containers/list.h>
#include <commons/log.h>
#include <bioformats/vcf/vcf_batch.h>
#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_file.h>

#include "hpg_vcf_tools_utils.h"
#include "split.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

int run_split(global_options_data_t *global_options_data, split_options_data_t *options_data);

static int initialize_output(cp_hashtable **output_files);

static void free_output(cp_hashtable *output_files);

static void free_file_key(char *key);

static void free_file_descriptor(FILE *fd);


#endif
