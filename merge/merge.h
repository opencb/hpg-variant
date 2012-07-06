#ifndef VCF_TOOLS_MERGE_H
#define VCF_TOOLS_MERGE_H

#include <getopt.h>
#include <stdlib.h>
#include <string.h>

#include <cprops/linked_list.h>
#include <libconfig.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <commons/log.h>
#include <commons/string_utils.h>
#include <containers/list.h>

#include "error.h"
#include "global_options.h"

#define NUM_MERGE_OPTIONS  1

typedef struct merge_options {
    struct arg_str *input_files;   /**< List of files used as input */
    int num_options;
} merge_options_t;

typedef struct merge_options_data {
    char **input_files;   /**< List of files used as input */
    int num_files;   /**< Number of files used as input */
} merge_options_data_t;


// typedef struct {
//     vcf_record_t *record;
//     char *merge_name;
// } merge_result_t;


static merge_options_t *new_merge_cli_options(void);

/**
 * Initialize a merge_options_data_t structure mandatory fields.
 */
static merge_options_data_t *new_merge_options_data(merge_options_t *options);

/**
 * Free memory associated to a merge_options_data_t structure.
 */
static void free_merge_options_data(merge_options_data_t *options_data);

// /**
//  * Initialize a variant_merge_result_t structure mandatory fields.
//  */
// merge_result_t *new_merge_result(vcf_record_t *record, char *merge_name);
// 
// /**
//  * Free memory associated to a variant_merge_result_t structure.
//  */
// void free_merge_result(merge_result_t* merge_result);


/* ******************************
 *       Tool execution         *
 * ******************************/

int merge(vcf_record_t **variants, int num_variants, list_t* output_list);

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
int read_merge_configuration(const char *filename, merge_options_t *options_data, shared_options_t *shared_options);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param global_options_data
 */
void **parse_merge_options(int argc, char *argv[], merge_options_t *options_data, shared_options_t *shared_options_data);

void **merge_merge_options(merge_options_t *merge_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * 
 * @param options_data
 */
int verify_merge_options(merge_options_t *options_data, shared_options_t *shared_options_data);



#endif
