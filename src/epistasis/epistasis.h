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

#ifndef EPISTASIS_H
#define EPISTASIS_H

/** 
 * @file epistasis.h
 * @brief Structures and functions associated to the options for the epistasis tool
 * 
 * This file defines the structures which store the options for the epistasis tool, and also 
 * functions to read their value from a configuration file or the command line.
 */ 

#include <libconfig.h>
#include <stdlib.h>

#include <argtable2.h>

#include <bioformats/vcf/vcf_filters.h>
#include <commons/log.h>

#include "error.h"
#include "shared_options.h"

/**
 * Number of options applicable to the epistasis tool.
 */
#define NUM_EPISTASIS_OPTIONS  2

typedef struct epistasis_options {
    int num_options;
    
    struct arg_lit *no_phenotypes; /**< Flag asking not to retrieve phenotypical information. */
    struct arg_str *excludes; /**< Comma-separated consequence types to exclude from the query. */
} epistasis_options_t;

/**
 * @brief Values for the options of the epistasis tool.
 * 
 * This struct contains the options specific to the epistasis tool, such as configuration of web 
 * services that retrieve phenotypical information, or the epistasiss that will not be retrieved.
 */
typedef struct epistasis_options_data {
    int no_phenotypes;  /**< Flag asking not to retrieve phenotypical information. */
    char *excludes;     /**< Comma-separated consequence types to exclude from the query. */
} epistasis_options_data_t;


static epistasis_options_t *new_epistasis_cli_options(void);

/**
 * @brief Initializes an epistasis_options_data_t structure mandatory members.
 * @return A new epistasis_options_data_t structure.
 * 
 * Initializes a new epistasis_options_data_t structure mandatory members, which are the buffers for 
 * the URL parts, as well as its numerical ones.
 */
static epistasis_options_data_t *new_epistasis_options_data(epistasis_options_t *options);

/**
 * @brief Free memory associated to a epistasis_options_data_t structure.
 * @param options_data the structure to be freed
 * 
 * Free memory associated to a epistasis_options_data_t structure, including its text buffers.
 */
static void free_epistasis_options_data(epistasis_options_data_t *options_data);


/* **********************************************
 *                Options parsing               *
 * **********************************************/

/**
 * @brief Reads the configuration parameters of the epistasis tool.
 * @param filename file the options data are read from
 * @param options_data local options values (host URL, species, num-threads...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 * 
 * Reads the basic configuration parameters of the tool. If the configuration
 * file can't be read, these parameters should be provided via the command-line
 * interface.
 */
int read_epistasis_configuration(const char *filename, epistasis_options_t *epistasis_options, shared_options_t *shared_options);

/**
 * @brief Parses the tool options from the command-line.
 * @param argc Number of arguments from the command-line
 * @param argv List of arguments from the command line
 * @param[out] options_data Struct where the tool-specific options are stored in
 * @param[out] global_options_data Struct where the application options are stored in
 * 
 * Reads the arguments from the command-line, checking they correspond to an option for the 
 * epistasis tool, and stores them in the local or global structure, depending on their scope.
 */
void **parse_epistasis_options(int argc, char *argv[], epistasis_options_t *epistasis_options, shared_options_t *shared_options);

void **merge_epistasis_options(epistasis_options_t *epistasis_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * @brief Checks semantic dependencies among the tool options.
 * @param global_options_data Application-wide options to check
 * @param options_data Tool-wide options to check
 * @return Zero (0) if the options are correct, non-zero otherwise
 * 
 * Checks that all dependencies among options are satisfied, i.e.: option A is mandatory, 
 * option B can't be provided at the same time as option C, and so on.
 */
int verify_epistasis_options(epistasis_options_t *epistasis_options, shared_options_t *shared_options);


#endif
