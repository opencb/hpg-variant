/*
 * Copyright (c) 2012-2014 Cristina Yenyxe Gonzalez Garcia (EMBL-EBI)
 * Copyright (c) 2012-2014 Ignacio Medina (EMBL-EBI)
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

#ifndef VCF_TOOLS_AGGREGATE_H
#define	VCF_TOOLS_AGGREGATE_H

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

//#include <bioformats/db/db_utils.h>
//#include <bioformats/ped/ped_file.h>
//#include <bioformats/vcf/vcf_db.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
//#include <bioformats/vcf/vcf_aggregate.h>
//#include <bioformats/vcf/vcf_aggregate_report.h>
#include <commons/file_utils.h>
#include <commons/log.h>
#include <commons/config/libconfig.h>
//#include <commons/sqlite/sqlite3.h>
#include <containers/khash.h>

#include "error.h"
#include "shared_options.h"
#include "hpg_variant_utils.h"

#define NUM_AGGREGATE_OPTIONS  12
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


typedef struct aggregate_options {
    struct arg_lit *overwrite;              /**< Whether to overwrite the values already in the INFO column. */
//    struct arg_str *info_fields;            /**< Attributes of the new INFO fields generated */
//    struct arg_lit *genotypes_count;        /**< Whether to calculate the genotypes counts. */
} aggregate_options_t;

/**
 * @struct aggregate_options_data
 * 
 */
typedef struct aggregate_options_data {
//    char **info_fields;     /**< List of attributes of the new INFO fields generated */
//    int num_info_fields;    /**< Number of attributes of the new INFO fields generated */
    int overwrite;          /**< Whether to overwrite the values already in the INFO column. */
} aggregate_options_data_t;


static aggregate_options_t *new_aggregate_cli_options(void);

/**
 * Initialize a aggregate_options_data_t structure mandatory fields.
 */
static aggregate_options_data_t *new_aggregate_options_data(aggregate_options_t *options);

/**
 * Free memory associated to a aggregate_options_data_t structure.
 */
static void free_aggregate_options_data(aggregate_options_data_t *options_data);


/* ******************************
 *       Tool execution         *
 * ******************************/

int run_aggregate(shared_options_data_t *shared_options_data, aggregate_options_data_t *options_data);


/* ******************************
 *      Options parsing         *
 * ******************************/

/**
 * Read the basic configuration parameters of the tool. If the configuration
 * file can't be read, these parameters should be provided via the command-line
 * interface.
 * 
 * @param filename File the options data are read from
 * @param options_data Local options values (filtering, sorting...) 
 * 
 * @return If the configuration has been successfully read
 */
int read_aggregate_configuration(const char *filename, aggregate_options_t *options, shared_options_t *shared_options);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param shared_options_data
 */
void **parse_aggregate_options(int argc, char *argv[], aggregate_options_t *aggregate_options, shared_options_t *shared_options);

void **merge_aggregate_options(aggregate_options_t *aggregate_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * 
 * @param options
 */
int verify_aggregate_options(aggregate_options_t *aggregate_options, shared_options_t *shared_options);



#endif
