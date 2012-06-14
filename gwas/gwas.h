#ifndef GWAS_H
#define GWAS_H

/** 
 * @file gwas.h
 * @brief Structures and functions associated to the options for the gwas tool
 * 
 * This file defines the structures which store the options for the gwas tool, and also 
 * functions to read their value from a configuration file or the command line.
 */ 

#include <getopt.h>
#include <libconfig.h>
#include <stdlib.h>

#include <argtable2.h>

#include <bioformats/vcf/vcf_filters.h>
#include <commons/log.h>

#include "error.h"
#include "global_options.h"
#include "main.h"

/**
 * Number of options applicable to the gwas tool.
 */
#define NUM_GWAS_OPTIONS  10

static struct option gwas_options[] = {
    // Task options
    { "assoc",                  no_argument, 0, 'o' },
    { "fisher",                 no_argument, 0, 'i' },
//     { "copy-number",            no_argument, 0, 'c' },
//     { "loh",                    no_argument, 0, 'l' },
    { "tdt",                    no_argument, 0, 't' },
    
    // Filtering options
    {"alleles",                 required_argument, 0, 'a' },
    {"coverage",                required_argument, 0, 'c' },
    {"quality",                 required_argument, 0, 'q' },

    {"region-file",             required_argument, 0, 'f' },
    {"region",                  required_argument, 0, 'r' },

    {"snp",                     required_argument, 0, 's' },
    
    // Multithreading options
    { "num-threads",            required_argument, 0, 'n' },
    
    // Other options
    
    {NULL,                      0, 0, 0}
};

enum GWAS_task { NONE, TDT, ASSOCIATION_BASIC, FISHER };


/**
 * @brief Values for the options of the gwas tool.
 * 
 * This struct contains the values for all the options of the gwas tool,
 * such as different parts of the web service URL or the parallelism 
 * parameters (number of threads, variants sent per request, and so on).
 */
typedef struct gwas_options_data
{
    struct arg_lit *assoc;
    struct arg_lit *fisher;
    struct arg_lit *tdt;
    
    enum GWAS_task task; /**< Task to perform */
    
//     /*  uint32_t??? */
//     long int num_threads; /**< Number of threads that query the web service simultaneously. */
//     long int max_batches; /**< Number of VCF records' batches that can be stored simultaneously. */
//     long int batch_size; /**< Maximum size of a VCF records' batch. */
//     long int variants_per_request; /**< Maximum number of variants sent in each web service query. */
//     
//     filter_chain *chain; /**< Chain of filters to apply to the VCF records, if that is the case. */
} gwas_options_data_t;


/**
 * @brief Initializes an gwas_options_data_t structure mandatory members.
 * @return A new gwas_options_data_t structure.
 * 
 * Initializes a new gwas_options_data_t structure mandatory members, which are the buffers for 
 * the URL parts, as well as its numerical ones.
 */
static gwas_options_data_t *init_options_data(void);

/**
 * @brief Free memory associated to a gwas_options_data_t structure.
 * @param options_data the structure to be freed
 * 
 * Free memory associated to a gwas_options_data_t structure, including its text buffers.
 */
static void free_options_data(gwas_options_data_t *options_data);


/* **********************************************
 *                Options parsing               *
 * **********************************************/

/**
 * @brief Reads the configuration parameters of the gwas tool.
 * @param filename file the options data are read from
 * @param options_data local options values (host URL, species, num-threads...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 * 
 * Reads the basic configuration parameters of the tool. If the configuration
 * file can't be read, these parameters should be provided via the command-line
 * interface.
 */
int read_gwas_configuration(const char *filename, gwas_options_data_t *options_data);

/**
 * @brief Parses the tool options from the command-line.
 * @param argc Number of arguments from the command-line
 * @param argv List of arguments from the command line
 * @param[out] options_data Struct where the tool-specific options are stored in
 * @param[out] global_options_data Struct where the application options are stored in
 * 
 * Reads the arguments from the command-line, checking they correspond to an option for the 
 * gwas tool, and stores them in the local or global structure, depending on their scope.
 */
void parse_gwas_options(int argc, char *argv[], gwas_options_data_t *options_data, shared_options_data_t *global_options_data);

/**
 * @brief Checks semantic dependencies among the tool options.
 * @param global_options_data Application-wide options to check
 * @param options_data Tool-wide options to check
 * @return Zero (0) if the options are correct, non-zero otherwise
 * 
 * Checks that all dependencies among options are satisfied, i.e.: option A is mandatory, 
 * option B can't be provided at the same time as option C, and so on.
 */
int verify_gwas_options(shared_options_data_t *global_options_data, gwas_options_data_t *options_data);


#endif
