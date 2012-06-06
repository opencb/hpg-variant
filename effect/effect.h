#ifndef EFFECT_H
#define EFFECT_H

/** 
 * @file effect.h
 * @brief Structures and functions associated to the options for the effect tool
 * 
 * This file defines the structures which store the options for the effect tool, and also 
 * functions to read their value from a configuration file or the command line.
 */ 

#include <getopt.h>
#include <libconfig.h>
#include <stdlib.h>

#include <bioformats/vcf/vcf_filters.h>
#include <commons/log.h>

#include "error.h"
#include "global_options.h"
#include "main.h"

/**
 * Number of options applicable to the effect tool.
 */
#define NUM_EFFECT_OPTIONS  9

static struct option effect_options[] = {
    // Filters applied to data
    {"alleles",                 required_argument, 0, 'a' },
    {"coverage",                required_argument, 0, 'c' },
    {"quality",                 required_argument, 0, 'q' },

    {"region-file",             required_argument, 0, 'f' },
    {"region",                  required_argument, 0, 'r' },

    {"snp",                     required_argument, 0, 's' },
    
    // Feature types to exclude
    { "exclude",                required_argument, 0, 'e' },
    
    // Multithreading options
    { "num-threads",            required_argument, 0, 'n' },
    { "variants-per-request",   required_argument, 0, 'p' },
    
    // Other options
    
    {NULL,                      0, 0, 0}
};


/**
 * @brief Values for the options of the effect tool.
 * 
 * This struct contains the values for all the options of the effect tool,
 * such as different parts of the web service URL or the parallelism 
 * parameters (number of threads, variants sent per request, and so on).
 */
typedef struct effect_options_data
{
	/*	uint32_t???	*/
	long int num_threads; /**< Number of threads that query the web service simultaneously. */
    long int max_batches; /**< Number of VCF records' batches that can be stored simultaneously. */
    long int batch_size; /**< Maximum size of a VCF records' batch. */
	long int variants_per_request; /**< Maximum number of variants sent in each web service query. */
	
// 	char *host_url; /**< URL of the host where the web service runs. */
// 	char *version; /**< Version of the WS to query. */
// 	char *species; /**< Species whose genome is taken as reference. */
	
    filter_chain *chain; /**< Chain of filters to apply to the VCF records, if that is the case. */
	char *excludes; /**< Consequence types to exclude from the query. */
} effect_options_data_t;


/**
 * @brief Initializes an effect_options_data_t structure mandatory members.
 * @return A new effect_options_data_t structure.
 * 
 * Initializes a new effect_options_data_t structure mandatory members, which are the buffers for 
 * the URL parts, as well as its numerical ones.
 */
static effect_options_data_t *init_options_data(void);

/**
 * @brief Free memory associated to a effect_options_data_t structure.
 * @param options_data the structure to be freed
 * 
 * Free memory associated to a effect_options_data_t structure, including its text buffers.
 */
static void free_options_data(effect_options_data_t *options_data);


/* **********************************************
 *                Options parsing               *
 * **********************************************/

/**
 * @brief Reads the configuration parameters of the effect tool.
 * @param filename file the options data are read from
 * @param options_data local options values (host URL, species, num-threads...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 * 
 * Reads the basic configuration parameters of the tool. If the configuration
 * file can't be read, these parameters should be provided via the command-line
 * interface.
 */
int read_effect_configuration(const char *filename, effect_options_data_t *options_data);

/**
 * @brief Parses the tool options from the command-line.
 * @param argc Number of arguments from the command-line
 * @param argv List of arguments from the command line
 * @param[out] options_data Struct where the tool-specific options are stored in
 * @param[out] global_options_data Struct where the application options are stored in
 * 
 * Reads the arguments from the command-line, checking they correspond to an option for the 
 * effect tool, and stores them in the local or global structure, depending on their scope.
 */
void parse_effect_options(int argc, char *argv[], effect_options_data_t *options_data, global_options_data_t *global_options_data);

/**
 * @brief Checks semantic dependencies among the tool options.
 * @param global_options_data Application-wide options to check
 * @param options_data Tool-wide options to check
 * @return Zero (0) if the options are correct, non-zero otherwise
 * 
 * Checks that all dependencies among options are satisfied, i.e.: option A is mandatory, 
 * option B can't be provided at the same time as option C, and so on.
 */
int verify_effect_options(global_options_data_t *global_options_data, effect_options_data_t *options_data);


#endif
