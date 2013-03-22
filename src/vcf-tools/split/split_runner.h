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

#ifndef SPLIT_RUNNER_H
#define SPLIT_RUNNER_H

#include <stdlib.h>

#include <omp.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_filters.h>
#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/list.h>
#include <containers/cprops/hashtable.h>

#include "hpg_variant_utils.h"
#include "split.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

int run_split(shared_options_data_t *shared_options_data, split_options_data_t *options_data);

static int initialize_output(cp_hashtable **output_files);

static void free_output(cp_hashtable *output_files);

static void free_file_key(char *key);

static void free_file_descriptor(FILE *fd);


#endif
