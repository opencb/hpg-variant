#ifndef VCF_TOOLS_H
#define VCF_TOOLS_H

#include <libconfig.h>

#include <file_utils.h>
#include <vcf_filters.h>

#include "error.h"
#include "global_options.h"

#define NUM_VCF_TOOLS_OPTIONS	3

static struct option vcf_tools_options[] = {
	// Filters applied to data
// 	  {"filter-quality",	required_argument, 0, 'q' },

	{"region-file",  required_argument, 0, 'f' },
	{"region",       required_argument, 0, 'r' },
    
	{"snp",          required_argument, 0, 's' },
	
	{NULL, 			0, 0, 0}
};

/**
 * @struct vcf_tools_options_data
 * 
 * @var vcf_tools_options_data::num_threads
 * Number of threads that to run simultaneously
 * 
 * @var vcf_tools_options_data::regions_descriptor
 * A string containing a list of regions or the name of the file where they must be read from
 * @var vcf_tools_options_data::use_region_file
 * Whether to get the regions directly from regions_descriptor or from a file
 */
typedef struct vcf_tools_options_data
{
// 	int use_region_file;
	long int num_threads;
    long int max_batches;
    long int batch_size;
	
// 	char *regions_descriptor;
    filter_chain *chain;
} vcf_tools_options_data_t;


/**
 * Initialize a vcf_tools_options_data_t structure mandatory fields.
 */
static vcf_tools_options_data_t *init_options_data();

/**
 * Free memory associated to a vcf_tools_options_data_t structure.
 */
static void free_options_data(vcf_tools_options_data_t *options_data);


/* **********************************************
 * 		Options parsing			*
 * **********************************************/

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
int read_vcf_tools_configuration(const char *filename, vcf_tools_options_data_t *options_data);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param global_options_data
 */
void parse_vcf_tools_options(int argc, char *argv[], vcf_tools_options_data_t *options_data, global_options_data_t *global_options_data);

/**
 * 
 * @param options_data
 */
int verify_vcf_tools_options(global_options_data_t *global_options_data, vcf_tools_options_data_t *options_data);



#endif
