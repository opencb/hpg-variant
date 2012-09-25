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

typedef struct split_options {
    struct arg_str *criterion;   /**< Criterion for splitting the file */
    int num_options;
} split_options_t;

typedef struct split_options_data {
    enum Split_criterion criterion;   /**< Criterion for splitting the file */
} split_options_data_t;


typedef struct {
    vcf_record_t *record;
    char *split_name;
} split_result_t;


static split_options_t *new_split_cli_options(void);

/**
 * Initialize a split_options_data_t structure mandatory fields.
 */
static split_options_data_t *new_split_options_data(split_options_t *options);

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

int split_by_chromosome(vcf_record_t **variants, int num_variants, list_t* output_list);

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
int read_split_configuration(const char *filename, split_options_t *options_data, shared_options_t *shared_options);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param global_options_data
 */
void **parse_split_options(int argc, char *argv[], split_options_t *options_data, shared_options_t *shared_options_data);

void **merge_split_options(split_options_t *split_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * 
 * @param options_data
 */
int verify_split_options(split_options_t *options_data, shared_options_t *shared_options_data);



#endif
