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
                           char **chromosome_order, int num_chromosomes, vcf_file_t **files, 
                           shared_options_data_t *shared_options_data, merge_options_data_t *options_data, list_t *output_list);

static void merge_remaining_interval(kh_pos_t* positions_read, vcf_file_t **files, 
                                     shared_options_data_t *shared_options_data, merge_options_data_t *options_data, list_t *output_list);




static void compose_key_value(const char *chromosome, const long position, char *key);

#endif
