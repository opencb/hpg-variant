#ifndef MERGE_RUNNER_H
#define MERGE_RUNNER_H

#include <assert.h>
#include <limits.h>
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

#include "hpg_variant_utils.h"
#include "merge.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

KHASH_MAP_INIT_STR(pos, array_list_t*);

int run_merge(shared_options_data_t *shared_options_data, merge_options_data_t *options_data);

static int insert_position_read(char key[64], vcf_record_file_link *link, kh_pos_t* positions_read); 

static void calculate_merge_interval(vcf_record_t* current_record, char** max_chromosome_merged, long unsigned int* max_position_merged,
                                     char **chromosome_order, int num_chromosomes);

static void merge_interval(kh_pos_t* positions_read, char *max_chromosome_merged, unsigned long max_position_merged, 
                           char **chromosome_order, int num_chromosomes, vcf_file_t *files, 
                           shared_options_data_t *shared_options_data, merge_options_data_t *options_data, list_t *output_list);

static void merge_remaining_interval(kh_pos_t* positions_read, vcf_file_t *files, 
                                     shared_options_data_t *shared_options_data, merge_options_data_t *options_data, list_t *output_list);




static void compose_key_value(const char *chromosome, const long position, char *key);

#endif
