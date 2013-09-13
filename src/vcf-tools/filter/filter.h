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

#ifndef VCF_TOOLS_FILTER_H
#define VCF_TOOLS_FILTER_H


#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_write.h>
#include <commons/file_utils.h>
#include <commons/log.h>
#include <commons/config/libconfig.h>
#include <containers/list.h>

#include "error.h"
#include "shared_options.h"
#include "hpg_variant_utils.h"

#define NUM_FILTER_OPTIONS  29

typedef struct filter_options {
    struct arg_lit *save_rejected;  /**< Flag that sets whether to write a file containing the rejected records */
} filter_options_t;

/**
 * @struct filter_options_data
 * 
 */
typedef struct filter_options_data {
    int save_rejected;      /**< Flag that sets whether to write a file containing the rejected records */
    filter_chain *chain;    /**< Chain of filters to apply to the VCF records. */
} filter_options_data_t;


static filter_options_t *new_filter_cli_options(void);

/**
 * Initialize a filter_options_data_t structure mandatory fields.
 */
static filter_options_data_t *new_filter_options_data(filter_options_t *options, shared_options_t *shared_options);

/**
 * Free memory associated to a filter_options_data_t structure.
 */
static void free_filter_options_data(filter_options_data_t *options_data);


/* ******************************
 *       Options parsing        *
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
int read_filter_configuration(const char *filename, filter_options_t *options, shared_options_t *shared_options);

/**
 * 
 * @param argc
 * @param argv
 * @param options
 * @param shared_options
 */
void **parse_filter_options(int argc, char *argv[], filter_options_t *filter_options, shared_options_t *shared_options);

void **merge_filter_options(filter_options_t *filter_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * 
 * @param options
 */
int verify_filter_options(filter_options_t *filter_options, shared_options_t *shared_options);


#endif
