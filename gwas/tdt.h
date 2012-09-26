/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (CGI-CIPF)
 * Copyright (c) 2012 Ignacio Medina (CGI-CIPF)
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

typedef struct {
    char *chromosome;
    char *reference;
    char *alternate;
    
    unsigned long int position;
    
    int t1;
    int t2;
    double odds_ratio;
    double chi_square;
    double p_value;
} tdt_result_t;

int tdt_test(vcf_record_t **variants, int num_variants, family_t **families, int num_families, cp_hashtable *sample_ids, list_t *output_list);


//tdt_result_t *tdt_result_new(char *chromosome, unsigned long int position, char *reference, char *alternate, double t1, double t2, double chi_square);
tdt_result_t* tdt_result_new(char *chromosome, int chromosome_len, unsigned long int position, char *reference, int reference_len,
                             char *alternate, int alternate_len, double t1, double t2, double chi_square);

void tdt_result_free(tdt_result_t *result);

#endif
