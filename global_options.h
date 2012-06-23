#ifndef GLOBAL_OPTIONS_H
#define GLOBAL_OPTIONS_H

/** 
 * @file shared_options.h 
 * @brief Structures and functions associated to application-wide options
 * 
 * This file defines the structures which store the options for the whole application, and also 
 * functions to read their value from a configuration file or the command line.
 */ 

#include <stdlib.h>
#include <string.h>
//#include <getopt.h>

#include <argtable2.h>
#include <libconfig.h>

#include <bioformats/vcf/vcf_util.h>
#include <commons/log.h>

#include "error.h"

/**
 * Number of options applicable to the whole application.
 */
#define NUM_GLOBAL_OPTIONS	13
    
/*static struct option shared_options[] = {
    // File formats accepted (range available A-H)
    {"vcf-file",        required_argument, 0, 'A' },
//  {"bam-file",        required_argument, 0, 'B' },
    {"ped-file",        required_argument, 0, 'E' },
//    {"gff-file",      required_argument, 0, 'G' },
    
    // IO options (range available I-O)
    {"outdir",          required_argument, 0, 'N' },
    {"out",             required_argument, 0, 'O' },
    
    // Other options (range available P-Z)
    { "species",        required_argument, 0, 'S' },
    { "version",        required_argument, 0, 'V' },
    { "url",            required_argument, 0, 'U' },
    
    {NULL,          0, 0, 0}
};*/

typedef struct shared_options {
    struct arg_file *vcf_filename; /**< VCF file used as input. */
    struct arg_file *ped_filename; /**< PED file used as input. */
    struct arg_file *output_filename; /**< Filename template for the main output file. */
    struct arg_str *output_directory; /**< Directory where the output files will be stored. */

    struct arg_str *host_url; /**< URL of the host where the web service runs. */
    struct arg_str *version; /**< Version of the WS to query. */
    struct arg_str *species; /**< Species whose genome is taken as reference. */

    struct arg_int *max_batches; /**< Maximum number of batches stored at the same time. */
    struct arg_int *batch_size; /**< Maximum size of a batch. */
    struct arg_int *num_threads; /**< Number of threads when a task is perform in parallel. */
    struct arg_int *entries_per_thread; /**< Number of entries in a batch each thread processes. */

    struct arg_file *config_file; /**< Path to the configuration file */

    struct arg_lit *mmap_vcf_files; /**< Whether to use map VCF files to virtual memory or use the I/O API. */

    int num_options;
} shared_options_t;


/**
 * @brief Values for the application-wide options.
 * 
 * This struct contains the values for all the options of the application that can be applied to 
 * all its tools, such as filenames for different formats (VCF, GFF, BAM...) and the output files 
 * and folders.
 */
typedef struct shared_options_data
{
    char *ped_filename; /**< PED file used as input. */
    char *vcf_filename; /**< VCF file used as input. */
    char *output_directory; /**< Directory where the output files will be stored. */
    char *output_filename; /**< Filename template for the main output file. */
    
    char *host_url;
    char *version;
    char *species;

    int max_batches;
    int batch_size;
    int num_threads;
    int entries_per_thread;
} shared_options_data_t;

/**
 * @brief Initializes an global_options_t structure mandatory members.
 * @return A new global_options_t structure.
 * 
 * Initializes the only mandatory member of a global_options_t, which is the output directory.
 */
shared_options_t *new_shared_cli_options(void);

/**
 * @brief Initializes an global_options_data_t structure mandatory members.
 * @return A new global_options_data_t structure.
 * 
 * Initializes the only mandatory member of a global_options_data_t, which is the output directory.
 */
shared_options_data_t *new_shared_options_data(shared_options_t *options);

/**
 * @brief Free memory associated to a global_options_data_t structure.
 * @param options_data the structure to be freed
 * 
 * Free memory associated to a global_options_data_t structure, including its text buffers.
 */
void free_shared_options_data(shared_options_data_t *options_data);


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
int read_shared_configuration(const char *filename, shared_options_t *options);

#endif
