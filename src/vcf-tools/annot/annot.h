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

#ifndef VCF_TOOLS_ANNOT_H
#define VCF_TOOLS_ANNOT_H

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

#define NUM_ANNOT_OPTIONS  3
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


typedef struct annot_options {
    struct arg_lit *variant_annot;      /**< Whether to get annot about variants. */
    struct arg_lit *sample_annot;       /**< Whether to get annot about samples. */
    struct arg_lit *save_db;            /**< Whether to save annot to a database. */
    
    int num_options;
} annot_options_t;

/**
 * @struct annot_options_data
 * 
 */
typedef struct annot_options_data {
    int variant_annot;  /**< Whether to get annot about variants. */
    int sample_annot;   /**< Whether to get annot about samples. */
    int save_db;        /**< Whether to save annot to a database. */
} annot_options_data_t;


typedef struct vcf_annot_sample{
    char *name;
    array_list_t *chromosomes;
} vcf_annot_sample_t;

typedef struct vcf_annot_chr{
    char *name;
    array_list_t *positions;
} vcf_annot_chr_t;

typedef struct vcf_annot_pos{
    unsigned int pos;
    short int dp;
} vcf_annot_pos_t;

static annot_options_t *new_annot_cli_options(void);

/**
 * Initialize a annot_options_data_t structure mandatory fields.
 */
static annot_options_data_t *new_annot_options_data(annot_options_t *options);

/**
 * Free memory associated to a annot_options_data_t structure.
 */
static void free_annot_options_data(annot_options_data_t *options_data);


/* ******************************
 *       Tool execution         *
 * ******************************/

int run_annot(shared_options_data_t *shared_options_data, annot_options_data_t *options_data);


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
int read_annot_configuration(const char *filename, annot_options_t *options, shared_options_t *shared_options);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param shared_options_data
 */
void **parse_annot_options(int argc, char *argv[], annot_options_t *annot_options, shared_options_t *shared_options);

void **merge_annot_options(annot_options_t *annot_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * 
 * @param options
 */
int verify_annot_options(annot_options_t *annot_options, shared_options_t *shared_options);



#endif
