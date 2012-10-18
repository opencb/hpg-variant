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

#ifndef HPG_VARIANT_UTIL_H
#define HPG_VARIANT_UTIL_H

#include <math.h>
#include <stdlib.h>

#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/list.h>

#include "shared_options.h"

/* ***********************
 *     Initialization    *
 * ***********************/

char *find_configuration_file(int argc, char *argv[]);


/* **********************
 *    Job management    *
 * **********************/

FILE *new_job_status_file(char *path);

void update_job_status_file(int percentage, FILE *file);

void close_job_status_file(FILE *file);


/* ***********************
 *      Miscellaneous    *
 * ***********************/

void show_usage(char *tool, void **argtable, int num_arguments);

/**
 * @brief Creates output files depending on the tool options given as input
 * @param global_options_data Options that apply to decide whether to create the output files or not
 * @param passed_file File containing entries that pass filters tests
 * @param failed_file File containing entries that don't pass filters tests
 * @return 0 if no errors occurred, 1 otherwise
 * 
 * If an output filename has been provided as command-line argument, creates 2 files with that prefix. 
 * If no output filename has been provided but some filters are applied, uses the original filename as prefix.
 * If no output filename nor filters are provided, don't create output files.
 */
int get_output_files(shared_options_data_t *global_options_data, FILE **passed_file, FILE **failed_file);

/**
 * @brief Given a list of records, distributes them in chunks of similar size
 * @param records list of records to separate in chunks
 * @param max_chunk_size maximum size of a chunk
 * @param[out] num_chunks number of chunks created
 * @param[out] chunk_sizes size of each chunk
 * @return The indices of the beginning of each chunk
 * 
 * Given a list of records, defines another list of chunks whose elements point to records separated 
 * by a maximum distance of max_chunk_size. These records will mark the beginning of each chunk.
 */
int *create_chunks(int length, int max_chunk_size, int *num_chunks, int **chunk_sizes);


#endif
