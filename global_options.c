#include "global_options.h"

shared_options_t *new_shared_cli_options(void) {
    shared_options_t *options_data = (shared_options_t*) calloc (1, sizeof(shared_options_t));
    
    options_data->vcf_filename = arg_file1(NULL, "vcf-file", NULL, "VCF file used as input");
    options_data->ped_filename = arg_file0(NULL, "ped-file", NULL, "PED file used as input");
    options_data->output_filename = arg_file0(NULL, "out", NULL, "Filename prefix for main output files");
    options_data->output_directory = arg_str0(NULL, "outdir", NULL, "Directory where the output files will be stored");
    
    options_data->host_url = arg_str0(NULL, "url", NULL, "URL of the host where remote web services run");
    options_data->version = arg_str0(NULL, "version", NULL, "Version of the web service to query");
    options_data->species = arg_str0(NULL, "species", NULL, "Species whose genome is taken as reference");
    
    options_data->max_batches = arg_int0(NULL, "num-batches", NULL, "Maximum number of batches stored at the same time");
    options_data->batch_size = arg_int0(NULL, "batch-size", NULL, "Maximum size of a batch");
    options_data->num_threads = arg_int0(NULL, "num-threads", NULL, "Number of threads when a task is perform in parallel");
    options_data->entries_per_thread = arg_int0(NULL, "entries-per-thread", NULL, "Number of entries in a batch each thread processes");
    
    options_data->num_alleles = arg_int0(NULL, "alleles", NULL, "Filter: by number of alleles");
    options_data->coverage = arg_int0(NULL, "coverage", NULL, "Filter: by minimum coverage");
    options_data->quality = arg_int0(NULL, "quality", NULL, "Filter: by minimum quality");
    options_data->region = arg_str0(NULL, "region", NULL, "Filter: by a list of regions (chr1:start1-end1,chr2:start2-end2...)");
    options_data->region_file = arg_file0(NULL, "region-file", NULL, "Filter: by a list of regions (read from a GFF file)");
    options_data->snp = arg_str0(NULL, "snp", NULL, "Filter: by being a SNP or not");
    
    options_data->config_file = arg_file0(NULL, "config", NULL, "File that contains the parameters for configuring the application");
    
    options_data->num_options = NUM_GLOBAL_OPTIONS;
    
    return options_data;
}

shared_options_data_t* new_shared_options_data(shared_options_t* options) {
    shared_options_data_t *options_data = (shared_options_data_t*) calloc (1, sizeof(shared_options_data_t));
    
    options_data->vcf_filename = strdup(*(options->vcf_filename->filename));
    options_data->ped_filename = strdup(*(options->ped_filename->filename));
    options_data->output_filename = strdup(*(options->output_filename->filename));
    options_data->output_directory = strdup(*(options->output_directory->sval));
    
    options_data->host_url = strdup(*(options->host_url->sval));
    options_data->version = strdup(*(options->version->sval));
    options_data->species = strdup(*(options->species->sval));
    
    options_data->max_batches = *(options->max_batches->ival);
    options_data->batch_size = *(options->batch_size->ival);
    options_data->num_threads = *(options->num_threads->ival);
    options_data->entries_per_thread = *(options->entries_per_thread->ival);
    
    options_data->num_alleles = *(options->num_alleles->ival);
    options_data->coverage = *(options->coverage->ival);
    options_data->quality = *(options->quality->ival);
    options_data->region = options->region->count > 0 ? strdup(*(options->region->sval))
                                                      : ( options->region_file->count > 0 ? strdup(*(options->region_file->filename)) : NULL );
    options_data->snp = strdup(*(options->snp->sval));
    
    return options_data;
}


void free_shared_options_data(shared_options_data_t *options_data) {
    if (options_data->vcf_filename)     { free(options_data->vcf_filename); }
    if (options_data->output_directory) { free(options_data->output_directory); }
    if (options_data->output_filename)  { free(options_data->output_filename); }
    free(options_data);
}

int read_global_configuration(const char *filename, shared_options_t *options_data) {
    if (filename == NULL || options_data == NULL) {
        return -1;
    }
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    
    // Read output directory
    const char *tmp_string;
    ret_code = config_lookup_string(config, "global.outdir", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Output folder not found in configuration file, must be set via command-line");
    } else {
        *(options_data->output_directory->sval) = strdup(tmp_string);
        LOG_INFO_F("Output folder = %s (%zu chars)\n",
                   *(options_data->output_directory->sval), 
                   strlen(*(options_data->output_directory->sval)));
    }
    
    config_destroy(config);
    free(config);
    
    return ret_code;
}
