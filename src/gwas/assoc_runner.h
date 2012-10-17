/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
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

#ifndef ASSOC_RUNNER_H
#define ASSOC_RUNNER_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cprops/hashtable.h>
#include <omp.h>

#include <bioformats/family.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/ped/ped_file_structure.h>
// #include <bioformats/vcf/vcf_batch.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_reader.h>
#include <bioformats/vcf/vcf_util.h>
#include <commons/log.h>
#include <commons/string_utils.h>
#include <containers/list.h>
#include <math/fisher.h>

#include "assoc.h"
#include "assoc_basic_test.h"
#include "shared_options.h"
#include "gwas.h"
#include "hpg_variant_utils.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


int run_association_test(shared_options_data_t *global_options_data, gwas_options_data_t *options_data);


static cp_hashtable *associate_samples_and_positions(vcf_file_t *file);

#endif
