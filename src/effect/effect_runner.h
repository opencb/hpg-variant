/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (CGI-CIPF)
 * Copyright (c) 2012 Ignacio Medina (CGI-CIPF)
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

#include <cprops/hashtable.h>
#include <cprops/linked_list.h>
#include <curl/curl.h>
#include <omp.h>

// #include <bioformats/vcf/vcf_batch.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_write.h>
#include <commons/file_utils.h>
#include <commons/http_utils.h>
#include <commons/log.h>
#include <commons/result.h>
#include <commons/string_utils.h>
#include <containers/array_list.h>
#include <containers/list.h>

#include "effect.h"
#include "error.h"
#include "hpg_variant_utils.h"
#include "main.h"

#define CONSEQUENCE_TYPE_WS_NUM_PARAMS  3
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

enum phenotype_source { SNP_PHENOTYPE = -1, MUTATION_PHENOTYPE = -2};

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
 *              Web service request             *
 * **********************************************/

/**
 * @brief Given a list of arguments, compounds a URL to invoke a web service.
 * @param options_data options that define some request arguments
 * @return URL composed from the given arguments, or NULL if any of the mandatory arguments is NULL
 * 
 * Given a list of arguments, compounds a URL to invoke a web service.
 */
char *compose_effect_ws_request(const char *category, const char *method, shared_options_data_t *options_data);

/**
 * @brief Invokes the effect web service for a list of regions.
 * 
 * @param url URL to invoke the web service through
 * @param records VCF records whose variant effect will be predicted
 * @param num_regions number of regions
 * @param excludes consequence types to exclude from the response
 * @return Whether the request could be successfully serviced
 * 
 * Given a list of genome positions, invokes the web service that returns a list of effect or consequences 
 * of the mutations in them. A callback function in order to parse the response.
 */
int invoke_effect_ws(const char *url, vcf_record_t **records, int num_records, char *excludes);

int invoke_snp_phenotype_ws(const char *url, vcf_record_t **records, int num_records);

int invoke_mutation_phenotype_ws(const char *url, vcf_record_t **records, int num_records);


/* **********************************************
 *              Response management             *
 * **********************************************/

/**
 * @brief Parses the response from the effect web service.
 * @param contents buffer with the response text
 * @param size size of an element in the buffer
 * @param nmemb number of elements in the buffer
 * @param userdata[out] where to write the response to (non-used)
 * @return The number of bytes processed
 * 
 * Reads the contents of the response from the effect web service
 */
static size_t write_effect_ws_results(char *contents, size_t size, size_t nmemb, void *userdata);

static size_t write_snp_phenotype_ws_results(char *contents, size_t size, size_t nmemb, void *userdata);

static size_t write_mutation_phenotype_ws_results(char *contents, size_t size, size_t nmemb, void *userdata);

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
 * @param num_threads the number of threads that parse the response
 * @param output_directory the directory where the output files will be written
 * @return Whether all the structures were successfully initialized
 * 
 * Initialize the structures for storing the web service response and writing it to the output files.
 */
int initialize_ws_output(shared_options_data_t *global_options_data, effect_options_data_t *options_data);

/**
 * @param num_threads the number of threads that parsed the response
 * @return Whether all the structures were successfully freed
 * 
 * Free the structures for storing the web service response.
 */
int free_ws_output(int num_threads);

/**
 * @param key the unique identifier in a pair (id, file descriptor)
 * 
 * Frees the unique identifier for storing a file descriptor in a map.
 */
static void free_file_key1(int *key);
// static void free_file_key1(char *key);

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
