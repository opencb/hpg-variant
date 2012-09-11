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
    options_data->batch_lines = arg_int0(NULL, "batch-lines", NULL, "Maximum number of lines in a batch");
    options_data->batch_bytes = arg_int0(NULL, "batch-bytes", NULL, "Maximum number of bytes in a batch");
    options_data->num_threads = arg_int0(NULL, "num-threads", NULL, "Number of threads when a task is perform in parallel");
    options_data->entries_per_thread = arg_int0(NULL, "entries-per-thread", NULL, "Number of entries in a batch each thread processes");
    
    options_data->num_alleles = arg_int0(NULL, "alleles", NULL, "Filter: by number of alleles");
    options_data->coverage = arg_int0(NULL, "coverage", NULL, "Filter: by minimum coverage");
    options_data->quality = arg_int0(NULL, "quality", NULL, "Filter: by minimum quality");
    options_data->region = arg_str0(NULL, "region", NULL, "Filter: by a list of regions (chr1:start1-end1,chr2:start2-end2...)");
    options_data->region_file = arg_file0(NULL, "region-file", NULL, "Filter: by a list of regions (read from a GFF file)");
    options_data->snp = arg_str0(NULL, "snp", NULL, "Filter: by being a SNP or not");
    
    options_data->config_file = arg_file0(NULL, "config", NULL, "File that contains the parameters for configuring the application");
    
    options_data->mmap_vcf_files = arg_lit0(NULL, "mmap-vcf", "Whether to use map VCF files to virtual memory or use the I/O API");
    
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
    options_data->batch_lines = *(options->batch_lines->ival);
    options_data->batch_bytes = *(options->batch_bytes->ival);
    options_data->num_threads = *(options->num_threads->ival);
    options_data->entries_per_thread = *(options->entries_per_thread->ival);
    
    filter_t *filter;
    if (options->num_alleles->count > 0) {
        filter = create_num_alleles_filter(*(options->num_alleles->ival));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("number of alleles filter = %d\n", *(options->num_alleles->ival));
    }
    if (options->coverage->count > 0) {
        filter = create_coverage_filter(*(options->coverage->ival));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("minimum coverage filter = %d\n", *(options->coverage->ival));
    }
    if (options->quality->count > 0) {
        filter = create_quality_filter(*(options->quality->ival));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("minimum quality filter = %d\n", *(options->quality->ival));
    }
    if (options->snp->count > 0) {
        filter = create_snp_filter(strdup(*(options->snp->sval)));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("snp filter to %s SNPs\n", *(options->snp->sval));
    }
    if (options->region->count > 0) {
        filter = create_region_exact_filter(strdup(*(options->region->sval)), 0, options_data->host_url, options_data->species, options_data->version);
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("regions = %s\n", *(options->region->sval));
    } 
    if (options->region_file->count > 0) {
        filter = create_region_exact_filter(strdup(*(options->region_file->filename)), 1, options_data->host_url, options_data->species, options_data->version);
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("regions file = %s\n", *(options->region->sval));
    }
    
    // If not previously defined, set the value present in the command-line
    if (!mmap_vcf) {
        mmap_vcf = options->mmap_vcf_files->count;
    }
    
    return options_data;
}


void free_shared_options_data(shared_options_data_t *options_data) {
    if (options_data->vcf_filename)     { free(options_data->vcf_filename); }
    // TODO ped filename freed in ped_close
//     if (options_data->ped_filename)     { free(options_data->ped_filename); }
    if (options_data->output_directory) { free(options_data->output_directory); }
    if (options_data->output_filename)  { free(options_data->output_filename); }
    if (options_data->host_url)         { free(options_data->host_url); }
    if (options_data->version)          { free(options_data->version); }
    if (options_data->species)          { free(options_data->species); }
    free(options_data);
}

int read_global_configuration(const char *filename, shared_options_t *options_data) {
    if (filename == NULL || options_data == NULL) {
        return -1;
    }
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    
    const char *tmp_string;
    
    // Read output directory
    ret_code = config_lookup_string(config, "global.outdir", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Output folder not found in configuration file, must be set via command-line");
    } else {
        *(options_data->output_directory->sval) = strdup(tmp_string);
        LOG_DEBUG_F("Output folder = %s (%zu chars)\n",
                   *(options_data->output_directory->sval), 
                   strlen(*(options_data->output_directory->sval)));
    }
    
    // Read whether to mmap VCF files
    ret_code = config_lookup_bool(config, "global.mmap-vcf", &mmap_vcf);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("I/O strategy for VCF files not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("VCF files mapped to virtual memory = %d\n", mmap_vcf);
    }
    
    config_destroy(config);
    free(config);
    
    return ret_code;
}
