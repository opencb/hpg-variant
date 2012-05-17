#ifndef VCF_TOOLS_FILTER_H
#define VCF_TOOLS_FILTER_H


#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include <libconfig.h>
#include <omp.h>

#include <file_utils.h>
#include <list.h>
#include <log.h>
#include <vcf_batch.h>
#include <vcf_filters.h>
#include <vcf_file.h>
#include <vcf_write.h>

#include "error.h"
#include "global_options.h"

#define NUM_FILTER_OPTIONS  4


static struct option filter_options[] = {
    // Filters applied to data
    {"quality",     required_argument, 0, 'q' },

    {"region-file", required_argument, 0, 'f' },
    {"region",      required_argument, 0, 'r' },

    {"snp",         required_argument, 0, 's' },

    {NULL,          0, 0, 0}
};

/**
 * @struct filter_options_data
 * 
 */
typedef struct filter_options_data
{
	long int num_threads;   /**< Number of threads that to run simultaneously. */
    long int max_batches;   /**< Number of VCF records' batches that can be stored simultaneously. */
    long int batch_size;    /**< Maximum size of a VCF records' batch. */

    filter_chain *chain;    /**< Chain of filters to apply to the VCF records. */
} filter_options_data_t;


/**
 * Initialize a filter_options_data_t structure mandatory fields.
 */
static filter_options_data_t *new_filter_options_data();

/**
 * Free memory associated to a filter_options_data_t structure.
 */
static void free_filter_options_data(filter_options_data_t *options_data);


/* ******************************
 *       Tool execution         *
 * ******************************/

int run_filter(global_options_data_t *global_options_data, filter_options_data_t *options_data);


/* ******************************
 *       Options parsing        *
 * ******************************/

/**
 * Read the basic configuration parameters of the tool. If the configuration
 * file can't be read, these parameters should be provided via the command-line
 * interface.
 * 
 * @param filename File the options data are read from
 * @param options_data Local options values (filtering, sorting...) 
 * 
 * @return If the configuration has been successfully read
 */
int read_filter_configuration(const char *filename, filter_options_data_t *options_data);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param global_options_data
 */
void parse_filter_options(int argc, char *argv[], filter_options_data_t *options_data, global_options_data_t *global_options_data);

/**
 * 
 * @param options_data
 */
int verify_filter_options(global_options_data_t *global_options_data, filter_options_data_t *options_data);



#endif
