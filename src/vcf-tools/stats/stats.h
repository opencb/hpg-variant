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

#ifndef VCF_TOOLS_STATS_H
#define VCF_TOOLS_STATS_H

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include <bioformats/db/db_utils.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/vcf/vcf_db.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_stats.h>
#include <bioformats/vcf/vcf_stats_report.h>
#include <commons/file_utils.h>
#include <commons/log.h>
#include <commons/config/libconfig.h>
#include <commons/sqlite/sqlite3.h>
#include <containers/khash.h>

#include "error.h"
#include "shared_options.h"
#include "hpg_variant_utils.h"

#define NUM_STATS_OPTIONS  3
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


typedef struct stats_options {
    struct arg_lit *variant_stats;      /**< Whether to get stats about variants. */
    struct arg_lit *sample_stats;       /**< Whether to get stats about samples. */
    struct arg_lit *save_db;            /**< Whether to save stats to a database. */
    
    int num_options;
} stats_options_t;

/**
 * @struct stats_options_data
 * 
 */
typedef struct stats_options_data {
    int variant_stats;  /**< Whether to get stats about variants. */
    int sample_stats;   /**< Whether to get stats about samples. */
    int save_db;        /**< Whether to save stats to a database. */
} stats_options_data_t;


static stats_options_t *new_stats_cli_options(void);

/**
 * Initialize a stats_options_data_t structure mandatory fields.
 */
static stats_options_data_t *new_stats_options_data(stats_options_t *options);

/**
 * Free memory associated to a stats_options_data_t structure.
 */
static void free_stats_options_data(stats_options_data_t *options_data);


/* ******************************
 *       Tool execution         *
 * ******************************/

int run_stats(shared_options_data_t *shared_options_data, stats_options_data_t *options_data);


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
int read_stats_configuration(const char *filename, stats_options_t *options, shared_options_t *shared_options);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param shared_options_data
 */
void **parse_stats_options(int argc, char *argv[], stats_options_t *stats_options, shared_options_t *shared_options);

void **merge_stats_options(stats_options_t *stats_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * 
 * @param options
 */
int verify_stats_options(stats_options_t *stats_options, shared_options_t *shared_options);



#endif
