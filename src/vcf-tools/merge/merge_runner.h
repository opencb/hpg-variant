/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * hpg-variant is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-variant. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MERGE_RUNNER_H
#define MERGE_RUNNER_H

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

#include <omp.h>

#include <bioformats/features/region/region.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_reader.h>
#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/array_list.h>
#include <containers/khash.h>
#include <containers/list.h>
#include <containers/cprops/hashtable.h>

#include "hpg_variant_utils.h"
#include "merge.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

KHASH_MAP_INIT_STR(pos, array_list_t*);

int run_merge(shared_options_data_t *shared_options_data, merge_options_data_t *options_data);

static int insert_position_read(char key[64], vcf_record_file_link *link, kh_pos_t* positions_read); 

/**
 * @brief Given a record, checks that its position is before the last one marked for merging, and updates it in that case
 * @details 
 * 
 * @param current_record
 * @param max_chromosome_merged
 * @param max_position_merged
 * @param chromosome_order
 * @param num_chromosomes
 * @return Whether the max chromosome and position merged were updated
 */
static int calculate_merge_interval(vcf_record_t* current_record, char** max_chromosome_merged, long unsigned int* max_position_merged,
                                     char **chromosome_order, int num_chromosomes);

static int merge_interval(kh_pos_t* positions_read, char *max_chromosome_merged, unsigned long max_position_merged,
                           char **chromosome_order, int num_chromosomes, vcf_file_t **files, 
                           shared_options_data_t *shared_options_data, merge_options_data_t *options_data, list_t *output_list);

static int merge_remaining_interval(kh_pos_t* positions_read, vcf_file_t **files,
                                     shared_options_data_t *shared_options_data, merge_options_data_t *options_data, list_t *output_list);




static void compose_key_value(const char *chromosome, const long position, char *key);

static int record_cmp(const void *data1, const void *data2);

#endif
