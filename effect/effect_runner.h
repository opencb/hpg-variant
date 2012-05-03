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

#include <cprops/hashtable.h>
#include <curl/curl.h>
#include <omp.h>

#include <file_utils.h>
#include <http_utils.h>
#include <list.h>
#include <log.h>
#include <string_utils.h>
#include <vcf_batch.h>
#include <vcf_file.h>
#include <vcf_filters.h>
#include <vcf_write.h>

#include "effect.h"
#include "error.h"
#include "main.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

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
int execute_effect_query(char *url, global_options_data_t *global_options_data, effect_options_data_t *options_data);


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
char *compose_effect_ws_request(effect_options_data_t *options_data);

/**
 * @brief Given a list of records, defines the chunks to be send to the web service
 * @param records list of records to separate in chunks
 * @param max_chunk_size maximum size of a chunk
 * @param[out] num_chunks number of chunks created
 * @return The list of chunks
 * 
 * Given a list of records, defines another list of chunks whose elements point to records separated 
 * by a maximum distance of max_chunk_size. These records will mark the beginning of each chunk.
 */
list_item_t **create_chunks(list_t *records, int max_chunk_size, int *num_chunks);

/**
 * @brief Invokes the effect web service for a list of regions.
 * 
 * @param url URL to invoke the web service through
 * @param first_item first item of the list of regions
 * @param num_regions number of regions
 * @return Whether the request could be successfully serviced
 * 
 * Given a list of regions, invokes the web service that returns a list of effect or consequences 
 * of the mutations in them. A callback function in order to parse the response.
 */
int invoke_effect_ws(const char *url, list_item_t *first_item, int num_regions);


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
 * Reads the contents of the response from the effect web service, 
 */
static size_t write_effect_ws_results(char *contents, size_t size, size_t nmemb, void *userdata);

/**
 * Writes a summary file containing the number of entries for each of the consequence types processed.
 */
void write_summary_file(void);

/**
 * @param num_threads the number of threads that parse the response
 * @param output_directory the directory where the output files will be written
 * @return Whether all the structures were successfully initialized
 * 
 * Initialize the structures for storing the web service response and writing it to the output files.
 */
int initialize_ws_output(int num_threads, char *output_directory);

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
// static void free_file_key1(int *key);
static void free_file_key1(char *key);

/**
 * @param fd the file descriptor in a pair (id, file descriptor)
 * 
 * Frees a file descriptor stored in a map.
 */
static void free_file_descriptor(FILE *fd);

/**
 * @param count the number of entries in a pair (id, file descriptor)
 * 
 * Frees the pointer to the number of entries of a consequence type that have been read.
 */
static void free_summary_counter(int *count);


static int int_cmp(int *a, int *b);

#endif
