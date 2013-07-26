/*
 * Copyright (c) 2012-2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
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

#ifndef HPG_VARIANT_UTIL_H
#define HPG_VARIANT_UTIL_H

#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <bioformats/family/family.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/vcf/vcf_file.h>
#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/khash.h>
#include <containers/list.h>

#include "shared_options.h"

#define HPG_VARIANT_VERSION     "0.99.3"

/* ***********************
 *     Initialization    *
 * ***********************/

char *find_configuration_file(int argc, char *argv[]);

array_list_t *get_configuration_search_paths(int argc, char *argv[]);
char *get_config_path_from_args(int argc, char *argv[]);
char *get_config_home_folder(char *config_dirpaths[], int num_dirpaths);
array_list_t *sort_config_paths_by_priority(char *config_arg_path, char *home_path);


/* **********************
 *     Config files     *
 * **********************/

char *retrieve_config_file(char *filename, array_list_t *paths_to_search);


/* **********************
 *    Job management    *
 * **********************/

FILE *new_job_status_file(char *path);

void update_job_status_file(int percentage, FILE *file);

void close_job_status_file(FILE *file);


/* ***********************
 *        Filtering      *
 * ***********************/

/**
 * @brief Creates files that will contain the output of filtering tools, using the VCF input filename as prefix.
 * @param shared_options Options for deciding the destination of the filter output
 * @param passed_file File containing entries that pass filters tests
 * @param failed_file File containing entries that don't pass filters tests
 * @return 0 if no errors occurred, 1 otherwise
 */
int get_filtering_output_files(shared_options_data_t *shared_options, FILE **passed_file, FILE **failed_file);

int write_filtering_output_files(array_list_t *passed_records, array_list_t *failed_records, FILE* passed_file, FILE* failed_file);

array_list_t *filter_records(filter_t** filters, int num_filters, individual_t **individuals, khash_t(ids) *sample_ids, 
                             array_list_t *input_records, array_list_t **failed_records);

void free_filtered_records(array_list_t *passed_records, array_list_t *failed_records, array_list_t *input_records);


/* ***********************
 *         Output        *
 * ***********************/

FILE *get_output_file(shared_options_data_t *shared_options_data, char *default_name, char **path);


/* ***********************
 *      Miscellaneous    *
 * ***********************/

void show_usage(char *tool, void **argtable, int num_arguments);

void show_version(char *tool);

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
