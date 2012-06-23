#ifndef VCF_TOOLS_FILTER_H
#define VCF_TOOLS_FILTER_H


#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include <libconfig.h>
#include <omp.h>

#include <bioformats/vcf/vcf_batch.h>
#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_write.h>
#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/list.h>

#include "error.h"
#include "global_options.h"

#define NUM_FILTER_OPTIONS  6


// static struct option filter_options[] = {
//     // Filters applied to data
//     {"alleles",     required_argument, 0, 'a' },
//     
//     {"coverage",    required_argument, 0, 'c' },
//     
//     {"quality",     required_argument, 0, 'q' },
// 
//     {"region-file", required_argument, 0, 'f' },
//     {"region",      required_argument, 0, 'r' },
// 
//     {"snp",         required_argument, 0, 's' },
// 
//     {NULL,          0, 0, 0}
// };

typedef struct filter_options {
    struct arg_int *num_alleles; /**< Filter by number of alleles. */
    struct arg_int *coverage; /**< Filter by coverage. */
    struct arg_int *quality; /**< Filter by quality. */
    struct arg_str *region; /**< Filter by region */
    struct arg_file *region_file; /**< Filter by region (using a GFF file) */
    struct arg_str *snp; /**< Filter by SNP */
    
    int num_options;
} filter_options_t;

/**
 * @struct filter_options_data
 * 
 */
typedef struct filter_options_data
{
//     long int num_threads;   /**< Number of threads that to run simultaneously. */
//     long int max_batches;   /**< Number of VCF records' batches that can be stored simultaneously. */
//     long int batch_size;    /**< Maximum size of a VCF records' batch. */
    filter_chain *chain;    /**< Chain of filters to apply to the VCF records. */
} filter_options_data_t;


static filter_options_t *new_filter_cli_options(void);

/**
 * Initialize a filter_options_data_t structure mandatory fields.
 */
static filter_options_data_t *new_filter_options_data(filter_options_t *options, shared_options_t *shared_options);

/**
 * Free memory associated to a filter_options_data_t structure.
 */
static void free_filter_options_data(filter_options_data_t *options_data);


/* ******************************
 *       Tool execution         *
 * ******************************/

int run_filter(shared_options_data_t *shared_options_data, filter_options_data_t *options_data);


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
int read_filter_configuration(const char *filename, filter_options_t *options, shared_options_t *shared_options);

/**
 * 
 * @param argc
 * @param argv
 * @param options
 * @param shared_options
 */
void **parse_filter_options(int argc, char *argv[], filter_options_t *filter_options, shared_options_t *shared_options);
// void parse_filter_options(int argc, char *argv[], filter_options_t *options, shared_options_t *shared_options);

void **merge_filter_options(filter_options_t *filter_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * 
 * @param options
 */
int verify_filter_options(filter_options_t *filter_options, shared_options_t *shared_options);
// int verify_filter_options(shared_options_t *shared_options, filter_options_t *options);



#endif
