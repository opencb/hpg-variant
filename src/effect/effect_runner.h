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

#ifndef EFFECT_RUNNER_H
#define EFFECT_RUNNER_H

/** 
 * @file effect_runner.h
 * @brief Invocation of the effect web service
 * 
 * This file declares the functions that allow to invoke the web service to get the effect or consequences 
 * of a mutation. This process involves the composition of the URL, the invocation of the web service and 
 * the response management, creating a collection of output files.
 */ 

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <omp.h>

#include <bioformats/db/cellbase_connector.h>
#include <bioformats/features/variant/variant_effect.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_write.h>
#include <commons/file_utils.h>
#include <commons/http_utils.h>
#include <commons/log.h>
#include <commons/result.h>
#include <containers/array_list.h>
#include <containers/list.h>
#include <containers/cprops/hashtable.h>

#include "effect.h"
#include "error.h"
#include "hpg_variant_utils.h"

#define MAX_VARIANTS_PER_QUERY  1000
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

enum phenotype_source { SNP_PHENOTYPE = -1, MUTATION_PHENOTYPE = -2};

// Line buffers and their maximum size (one per thread)
extern char **effect_line, **snp_line, **mutation_line;
extern int *max_line_size, *snp_max_line_size, *mutation_max_line_size;


/**
 * @brief Performs the whole process of invocation of the effect web service and parsing of its output.
 * @param url URL of the web service
 * @param global_options_data options related to the whole application
 * @param options_data options related only to the effect tool
 * @return If the process was completed successfully
 * 
 * Producer-consumer implementation that allows to invoke the effect web service.
 * The producer reads a VCF file in batches of size given as an option in the options_data structure.
 * The consumer filters (optionally) the entries of each batch by certain conditions and invokes the 
 * web service for the ones that passed the filtering process. Finally, the response from the WS is 
 * processed and a file for the entries associated to each consequence type is created. There is 
 * also a summary file which contains the number of entries for each of these consequence types.
 */
int run_effect(char **urls, shared_options_data_t *global_options_data, effect_options_data_t *options_data);


/* **********************************************
 *              Response management             *
 * **********************************************/

/**
 * @brief Parses the response from the effect web service.
 * 
 * Reads the contents of the response from the effect web service
 */
static void parse_effect_response(int tid, char *output_directory, size_t output_directory_len, cp_hashtable *output_files, 
                                  list_t *output_list, cp_hashtable *summary_count, cp_hashtable *gene_list);

static void parse_snp_phenotype_response(int tid, list_t *output_list);

static void parse_mutation_phenotype_response(int tid, list_t *output_list);

/**
 * Writes a summary file containing the number of entries for each of the consequence types processed.
 */
void write_summary_file(cp_hashtable *summary_count, FILE *summary_file);

/**
 * Writes a file containing the list of genes with any variant taking place in them.
 */
void write_genes_with_variants_file(cp_hashtable *gene_list, char *output_directory);

/**
 * Writes an XML file containing the process input and output files, as well as some metadata about the 
 * process.
 */
void write_result_file(shared_options_data_t *global_options_data, effect_options_data_t *options_data, cp_hashtable *summary_count, char *output_directory);

/**
 * @param output_directory directory Where the files will be stored
 * @param output_directory_len length of the path to the output directory
 * @return Whether the output files descriptor where correctly initialized
 * 
 * Initialize the output files where the web service response will be written to.
 */
static int initialize_output_files(char *output_directory, size_t output_directory_len, cp_hashtable **output_files);

/**
 * @param num_threads the number of threads that parse the response
 * @param output_directory the directory where the output files will be written
 * 
 * Initialize the structures for storing the web service response and writing it to the output files.
 */
static void initialize_output_data_structures(shared_options_data_t *shared_options, list_t **output_list, cp_hashtable **summary_count, cp_hashtable **gene_list);

/**
 * 
 * Free the structures for storing the web service response.
 */
void free_output_data_structures(cp_hashtable *output_files, cp_hashtable *summary_count, cp_hashtable *gene_list);

/**
 * @param key the unique identifier in a pair (id, file descriptor)
 * 
 * Frees the unique identifier for storing a file descriptor in a map.
 */
static void free_file_key1(int *key);

/**
 * @param fd the file descriptor in a pair (id, file descriptor)
 * 
 * Frees a file descriptor stored in a map.
 */
static void free_file_descriptor(FILE *fd);

static void free_file_key2(char *key);

/**
 * @param count the number of entries in a pair (id, file descriptor)
 * 
 * Frees the pointer to the number of entries of a consequence type that have been read.
 */
static void free_summary_counter(int *count);


static int int_cmp(int *a, int *b);

#endif
