#ifndef VCF_TOOLS_SPLIT_H
#define VCF_TOOLS_SPLIT_H

#include <getopt.h>
#include <stdlib.h>
#include <string.h>

#include <cprops/linked_list.h>
#include <libconfig.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <commons/log.h>
#include <containers/list.h>

#include "error.h"
#include "global_options.h"

#define NUM_SPLIT_OPTIONS  1

enum Split_criterion { NONE, CHROMOSOME, GENE };

static struct option split_options[] = {
    {"criterion",    required_argument, 0, 'c' },
    
    {NULL,          0, 0, 0}
};

/**
 * @struct split_options_data
 * 
 */
typedef struct split_options_data {
    enum Split_criterion criterion;   /**< Criterion for splitting the file */

    long int max_batches;   /**< Number of VCF records' batches that can be stored simultaneously. */
    long int batch_size;   /**< Maximum size of a VCF records' batch. */
    long int num_threads;   /**< Number of threads that to run simultaneously. */
    long int variants_per_thread;   /**< Number of variants each thread will analyze. */
} split_options_data_t;


typedef struct {
    vcf_record_t *record;
    char *split_name;
} split_result_t;


/**
 * Initialize a split_options_data_t structure mandatory fields.
 */
static split_options_data_t *new_split_options_data();

/**
 * Free memory associated to a split_options_data_t structure.
 */
static void free_split_options_data(split_options_data_t *options_data);

/**
 * Initialize a variant_split_result_t structure mandatory fields.
 */
split_result_t *new_split_result(vcf_record_t *record, char *split_name);

/**
 * Free memory associated to a variant_split_result_t structure.
 */
void free_split_result(split_result_t* split_result);


/* ******************************
 *       Tool execution         *
 * ******************************/

int split_by_chromosome(list_item_t* variants, int num_variants, list_t* output_list);

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
int read_split_configuration(const char *filename, split_options_data_t *options_data);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param global_options_data
 */
void parse_split_options(int argc, char *argv[], split_options_data_t *options_data, global_options_data_t *global_options_data);

/**
 * 
 * @param options_data
 */
int verify_split_options(global_options_data_t *global_options_data, split_options_data_t *options_data);



#endif
