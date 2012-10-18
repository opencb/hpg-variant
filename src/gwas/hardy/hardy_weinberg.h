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

#ifndef TRANSMISSION_DISEQUILIBRIUM_TEST_H
#define TRANSMISSION_DISEQUILIBRIUM_TEST_H

#include <stdlib.h>
#include <string.h>

#include <cprops/hashtable.h>
#include <cprops/linked_list.h>
#include <gsl/gsl_cdf.h>
#include <omp.h>

#include <bioformats/ped/ped_file.h>
#include <bioformats/ped/ped_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_util.h>
#include <commons/log.h>
#include <containers/list.h>

#include "checks_family.h"

int hardy_weinberg_test(vcf_record_t **variants, int num_variants, individual_t **individuals, int num_individuals, 
                        cp_hashtable *sample_ids, list_t *output_list);

individual_t **get_founders_from_families(family_t **families, int num_families, int *num_individuals);

#endif