#ifndef GLOBAL_OPTIONS_H
#define GLOBAL_OPTIONS_H

/** 
 * @file global_options.h 
 * @brief Structures and functions associated to application-wide options
 * 
 * This file defines the structures which store the options for the whole application, and also 
 * functions to read their value from a configuration file or the command line.
 */ 

#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <libconfig.h>

#include <log.h>

/**
 * Number of options applicable to the whole application.
 */
#define NUM_GLOBAL_OPTIONS  3
    
static struct option global_options[] = {
    // File formats accepted (range available A-H)
    {"vcf-file",        required_argument, 0, 'A' },
//  {"bam-file",        required_argument, 0, 'B' },
    {"ped-file",      required_argument, 0, 'E' },
//    {"gff-file",      required_argument, 0, 'G' },
    
    // IO options (range available I-O)
    {"outdir",          required_argument, 0, 'N' },
    {"out",             required_argument, 0, 'O' },
    
    // Other options (range available P-Z)
    
    {NULL,          0, 0, 0}
};


/**
 * @brief Values for the application-wide options.
 * 
 * This struct contains the values for all the options of the application that can be applied to 
 * all its tools, such as filenames for different formats (VCF, GFF, BAM...) and the output files 
 * and folders.
 */
typedef struct global_options_data
{
    char *ped_filename; /**< PED file used as input. */
    char *vcf_filename; /**< VCF file used as input. */
    char *output_directory; /**< Directory where the output files will be stored. */
    char *output_filename; /**< Filename template for the main output file. */
} global_options_data_t;


/**
 * @brief Initializes an global_options_data_t structure mandatory members.
 * @return A new global_options_data_t structure.
 * 
 * Initializes the only mandatory member of a global_options_data_t, which is the output directory.
 */
global_options_data_t *init_global_options_data(void);

/**
 * @brief Free memory associated to a global_options_data_t structure.
 * @param options_data the structure to be freed
 * 
 * Free memory associated to a global_options_data_t structure, including its text buffers.
 */
void free_global_options_data(global_options_data_t *options_data);


/* **********************************************
 *                Options parsing               *
 * **********************************************/

/**
 * @brief Reads the configuration parameters of the application.
 * @param filename file the options data are read from
 * @param options_data options values (vcf filename, output directory...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 * 
 * Reads the basic configuration parameters of the application. If the configuration file can't be 
 * read, these parameters should be provided via the command-line interface.
 */
int read_global_configuration(const char *filename, global_options_data_t *options_data);

/**
 * @brief Given a list of options, parse the global ones it contains.
 * @param argc number of input token
 * @param argv input tokens
 * @param options_data structure for values of all the global options
 * @param start_index index of the next option to parse
 * @return The index where the parsing finished at
 * 
 * Given a list of options, parse the global ones it contains. Because it invokes
 * a nested getopt_long loop, the start_index argument sets the position of the next 
 * common option to parse.
 */
int parse_global_options(int argc, char *argv[], global_options_data_t *options_data, int start_index);

/**
 * @brief Checks semantic dependencies among the application options.
 * @param options_data Application-wide options to check
 * @return Zero (0) if the options are correct, non-zero otherwise
 * 
 * Checks that all dependencies among options are satisfied, i.e.: option A is mandatory, 
 * option B can't be provided at the same time as option C, and so on.
 */
int verify_global_options(global_options_data_t *options_data);

/**
 * @brief Merges application-wide options with the ones from a certain tool.
 * @param local_options options for a tool
 * @param num_local_options number of options for a tool
 * @return The global and local options merged in a unique structure
 * 
 * Given a set of options local to certain tool, merges them with the ones globally 
 * available to the whole application.
 */
struct option *merge_options(struct option local_options[], size_t num_local_options);

#endif
