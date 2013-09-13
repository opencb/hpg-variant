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

#include <gsl/gsl_cdf.h>
#include <omp.h>

#include <bioformats/family/checks_family.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/ped/ped_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_util.h>
#include <commons/log.h>
#include <commons/argtable/argtable2.h>
#include <commons/config/libconfig.h>
#include <containers/list.h>
#include <containers/khash.h>
#include <containers/cprops/hashtable.h>

#include "error.h"
#include "hpg_variant_utils.h"
#include "shared_options.h"

/**
 * Number of options applicable to the TDT tool.
 */
#define NUM_TDT_OPTIONS  29


typedef struct tdt_options { } tdt_options_t;

static tdt_options_t *new_tdt_cli_options(void);


/* **********************************************
 *                Options parsing               *
 * **********************************************/

/**
 * @brief Reads the configuration parameters of the tdt tool.
 * @param filename file the options data are read from
 * @param options_data local options values (host URL, species, num-threads...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 * 
 * Reads the basic configuration parameters of the tool. If the configuration
 * file can't be read, these parameters should be provided via the command-line
 * interface.
 */
int read_tdt_configuration(const char *filename, tdt_options_t *tdt_options, shared_options_t *shared_options);

/**
 * @brief Parses the tool options from the command-line.
 * @param argc Number of arguments from the command-line
 * @param argv List of arguments from the command line
 * @param[out] options_data Struct where the tool-specific options are stored in
 * @param[out] global_options_data Struct where the application options are stored in
 * 
 * Reads the arguments from the command-line, checking they correspond to an option for the 
 * tdt tool, and stores them in the local or global structure, depending on their scope.
 */
void **parse_tdt_options(int argc, char *argv[], tdt_options_t *tdt_options, shared_options_t *shared_options);

void **merge_tdt_options(tdt_options_t *tdt_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * @brief Checks semantic dependencies among the tool options.
 * @param global_options_data Application-wide options to check
 * @param options_data Tool-wide options to check
 * @return Zero (0) if the options are correct, non-zero otherwise
 * 
 * Checks that all dependencies among options are satisfied, i.e.: option A is mandatory, 
 * option B can't be provided at the same time as option C, and so on.
 */
int verify_tdt_options(tdt_options_t *tdt_options, shared_options_t *shared_options);


/* **********************************************
 *                Test execution                *
 * **********************************************/

typedef struct {
    char *chromosome;
    char *id;
    char *reference;
    char *alternate;
    
    unsigned long int position;
    
    int t1;
    int t2;
    double odds_ratio;
    double chi_square;
    double p_value;
} tdt_result_t;

int tdt_test(vcf_record_t **variants, int num_variants, family_t **families, int num_families, khash_t(ids) *sample_ids, list_t *output_list);

tdt_result_t* tdt_result_new(char *chromosome, int chromosome_len, unsigned long int position, char *id, int id_len, 
                             char *reference, int reference_len, char *alternate, int alternate_len, double t1, double t2, double chi_square);

void tdt_result_free(tdt_result_t *result);

#endif
