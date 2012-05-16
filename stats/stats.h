#ifndef VCF_TOOLS_STATS_H
#define VCF_TOOLS_STATS_H

#include <getopt.h>
#include <stdlib.h>
#include <string.h>

#include <cprops/linked_list.h>
#include <libconfig.h>

#include <file_utils.h>
#include <log.h>
#include <vcf_batch.h>
#include <vcf_filters.h>
#include <vcf_file.h>

#include "error.h"
#include "global_options.h"

#define NUM_STATS_OPTIONS  0


static struct option stats_options[] = {

    {NULL,          0, 0, 0}
};

/**
 * @struct stats_options_data
 * 
 */
typedef struct stats_options_data {
    long int num_threads;   /**< Number of threads that to run simultaneously. */
    long int max_batches;   /**< Number of VCF records' batches that can be stored simultaneously. */
    long int batch_size;   /**< Maximum size of a VCF records' batch. */
} stats_options_data_t;


typedef struct {
    char *ref_allele;
    char **alternates;
    
    int num_alleles;
    int *alleles_count;
    int *genotypes_count;
    
    int missing_alleles;
    int missing_genotypes;
} variant_stats_t;

/**
 * Initialize a stats_options_data_t structure mandatory fields.
 */
static stats_options_data_t *new_stats_options_data();

/**
 * Free memory associated to a stats_options_data_t structure.
 */
static void free_stats_options_data(stats_options_data_t *options_data);


/**
 * Initialize a variant_stats_t structure mandatory fields.
 */
variant_stats_t *new_variant_stats(char *ref_allele);

/**
 * Free memory associated to a variant_stats_t structure.
 */
void free_variant_stats(variant_stats_t *variant_stats);



/* ******************************
 *       Tool execution         *
 * ******************************/


int run_stats(global_options_data_t *global_options_data, stats_options_data_t *options_data);

int get_variants_stats(list_item_t *variants, int num_variants, list_t *output_list);


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
int read_stats_configuration(const char *filename, stats_options_data_t *options_data);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param global_options_data
 */
void parse_stats_options(int argc, char *argv[], stats_options_data_t *options_data, global_options_data_t *global_options_data);

/**
 * 
 * @param options_data
 */
int verify_stats_options(global_options_data_t *global_options_data, stats_options_data_t *options_data);



#endif
