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

#ifndef ASSOC_H
#define ASSOC_H

#include <gsl/gsl_cdf.h>
#include <omp.h>

#include <bioformats/family/family.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/ped/ped_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_util.h>
#include <commons/log.h>
#include <commons/argtable/argtable2.h>
#include <commons/config/libconfig.h>
#include <containers/list.h>
#include <containers/cprops/hashtable.h>

#include "assoc_basic_test.h"
#include "assoc_fisher_test.h"
#include "error.h"
#include "hpg_variant_utils.h"
#include "shared_options.h"


/**
 * Number of options applicable to the assoc tool.
 */
#define NUM_ASSOC_OPTIONS  31

typedef struct assoc_options {
    struct arg_lit *chisq;
    struct arg_lit *fisher;
} assoc_options_t;

enum ASSOC_task { NONE, CHI_SQUARE, FISHER };

/**
 * @brief Values for the options of the assoc tool.
 * 
 * This struct contains the values for all the options of the assoc tool,
 * such as different parts of the web service URL or the parallelism 
 * parameters (number of threads, variants sent per request, and so on).
 */
typedef struct assoc_options_data {
    enum ASSOC_task task; /**< Task to perform */
} assoc_options_data_t;


static assoc_options_t *new_assoc_cli_options(void);

/**
 * @brief Initializes an assoc_options_data_t structure mandatory members.
 * @return A new assoc_options_data_t structure.
 * 
 * Initializes a new assoc_options_data_t structure mandatory members, which are the buffers for 
 * the URL parts, as well as its numerical ones.
 */
static assoc_options_data_t *new_assoc_options_data(assoc_options_t *options);

/**
 * @brief Free memory associated to a assoc_options_data_t structure.
 * @param options_data the structure to be freed
 * 
 * Free memory associated to a assoc_options_data_t structure, including its text buffers.
 */
static void free_assoc_options_data(assoc_options_data_t *options_data);


/* **********************************************
 *                Options parsing               *
 * **********************************************/

/**
 * @brief Reads the configuration parameters of the assoc tool.
 * @param filename file the options data are read from
 * @param options_data local options values (host URL, species, num-threads...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 * 
 * Reads the basic configuration parameters of the tool. If the configuration
 * file can't be read, these parameters should be provided via the command-line
 * interface.
 */
int read_assoc_configuration(const char *filename, assoc_options_t *assoc_options, shared_options_t *shared_options);

/**
 * @brief Parses the tool options from the command-line.
 * @param argc Number of arguments from the command-line
 * @param argv List of arguments from the command line
 * @param[out] options_data Struct where the tool-specific options are stored in
 * @param[out] global_options_data Struct where the application options are stored in
 * 
 * Reads the arguments from the command-line, checking they correspond to an option for the 
 * assoc tool, and stores them in the local or global structure, depending on their scope.
 */
void **parse_assoc_options(int argc, char *argv[], assoc_options_t *assoc_options, shared_options_t *shared_options);

void **merge_assoc_options(assoc_options_t *assoc_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * @brief Checks semantic dependencies among the tool options.
 * @param global_options_data Application-wide options to check
 * @param options_data Tool-wide options to check
 * @return Zero (0) if the options are correct, non-zero otherwise
 * 
 * Checks that all dependencies among options are satisfied, i.e.: option A is mandatory, 
 * option B can't be provided at the same time as option C, and so on.
 */
int verify_assoc_options(assoc_options_t *assoc_options, shared_options_t *shared_options);


/* **********************************************
 *                Test execution                *
 * **********************************************/

//void assoc_test(enum ASSOC_task test_type, vcf_record_t **variants, int num_variants, family_t **families, int num_families,
//                cp_hashtable *sample_ids, const void *opt_input, list_t *output_list);
void assoc_test(enum ASSOC_task test_type, vcf_record_t **variants, int num_variants, individual_t **samples, int num_samples,
                const void *opt_input, list_t *output_list);

void assoc_count_individual(individual_t *individual, vcf_record_t *record, int allele1, int allele2, 
                           int *affected1, int *affected2, int *unaffected1, int *unaffected2);

#endif
