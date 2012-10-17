#ifndef VCF_TOOLS_STATS_H
#define VCF_TOOLS_STATS_H

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include <cprops/linked_list.h>
#include <libconfig.h>
#include <omp.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_stats.h>
#include <commons/file_utils.h>
#include <commons/log.h>

#include "error.h"
#include "shared_options.h"
#include "hpg_variant_utils.h"

#define NUM_STATS_OPTIONS  0
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


typedef struct stats_options {
    int num_options;
} stats_options_t;

/**
 * @struct stats_options_data
 * 
 */
typedef struct stats_options_data {
    // TODO no options right now, will be included when more statistics are available
} stats_options_data_t;


static stats_options_t *new_stats_cli_options(void);

/**
 * Initialize a stats_options_data_t structure mandatory fields.
 */
static stats_options_data_t *new_stats_options_data(stats_options_t *options);

/**
 * Free memory associated to a stats_options_data_t structure.
 */
static void free_stats_options_data(stats_options_data_t *options_data);


/* ******************************
 *       Tool execution         *
 * ******************************/

int run_stats(shared_options_data_t *shared_options_data, stats_options_data_t *options_data);


/* ******************************
 *      Options parsing         *
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
int read_stats_configuration(const char *filename, stats_options_t *options, shared_options_t *shared_options);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param shared_options_data
 */
void **parse_stats_options(int argc, char *argv[], stats_options_t *stats_options, shared_options_t *shared_options);

void **merge_stats_options(stats_options_t *stats_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * 
 * @param options
 */
int verify_stats_options(stats_options_t *stats_options, shared_options_t *shared_options);



#endif
