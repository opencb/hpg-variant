#include "effect.h"


int read_effect_configuration(const char *filename, effect_options_t *options_data, shared_options_t *global_options_data) {
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
    
    // Read number of threads that will make request to the web service
    ret_code = config_lookup_int(config, "effect.num-threads", global_options_data->num_threads->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of threads not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("num-threads = %ld\n", *(global_options_data->num_threads->ival));
    }

    // Read maximum number of batches that can be stored at certain moment
    ret_code = config_lookup_int(config, "effect.max-batches", global_options_data->max_batches->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Maximum number of batches not found in configuration file, must be set via command-line");
    } else {
        LOG_INFO_F("max-batches = %ld\n", *(global_options_data->max_batches->ival));
    }
    
    // Read size of every batch read
    ret_code = config_lookup_int(config, "effect.batch-size", global_options_data->batch_size->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Batch size not found in configuration file, must be set via command-line");
    } else {
        LOG_INFO_F("batch-size = %ld\n", *(global_options_data->batch_size->ival));
    }
    
    // Read number of variants per request to the web service
    ret_code = config_lookup_int(config, "effect.variants-per-request", global_options_data->entries_per_thread->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Variants per request not found in configuration file, must be set via command-line");
    } else {
        LOG_INFO_F("entries-per-thread = %ld\n", *(global_options_data->entries_per_thread->ival));
    }

    // Read host URL
    ret_code = config_lookup_string(config, "effect.url", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Web services URL not found in configuration file, must be set via command-line");
    } else {
//         free(global_options_data->host_url);
//         global_options_data->host_url = strdup(tmp_string);
        *(global_options_data->host_url->sval) = strdup(tmp_string);
        LOG_INFO_F("web services host URL = %s (%zu chars)\n",
                   *(global_options_data->host_url->sval), strlen(*(global_options_data->host_url->sval)));
    }
    
    // Read species
    ret_code = config_lookup_string(config, "effect.species", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Species not found in configuration file, must be set via command-line");
    } else {
//         free(global_options_data->species);
//         global_options_data->species = strdup(tmp_string);
        *(global_options_data->species->sval) = strdup(tmp_string);
        LOG_INFO_F("species = %s (%zu chars)\n",
                   *(global_options_data->species->sval), strlen(*(global_options_data->species->sval)));
    }
    
    // Read version
    ret_code = config_lookup_string(config, "effect.version", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Version not found in configuration file, must be set via command-line");
    } else {
//         free(global_options_data->version);
//         global_options_data->version = strdup(tmp_string);
        *(global_options_data->version->sval) = strdup(tmp_string);
        LOG_INFO_F("version = %s (%zu chars)\n",
                   *(global_options_data->version->sval), strlen(*(global_options_data->version->sval)));
    }

    config_destroy(config);
    free(config);

    return 0;
}

void **parse_effect_options(int argc, char *argv[], effect_options_t *options_data, shared_options_t *global_options_data) {
    struct arg_end *end = arg_end(options_data->num_options + global_options_data->num_options);
    void **argtable = merge_effect_options(options_data, global_options_data, end);
    
    int num_errors = arg_parse(argc, argv, argtable);
    if (num_errors > 0) {
        arg_print_errors(stdout, end, "hpg-variant");
    }
    
    return argtable;
}

void **merge_effect_options(effect_options_t *effect_options, shared_options_t *shared_options, struct arg_end *arg_end) {
    size_t opts_size = effect_options->num_options + shared_options->num_options + 1;
    void **tool_options = malloc (opts_size * sizeof(void*));
    tool_options[0] = shared_options->vcf_filename;
    tool_options[1] = shared_options->ped_filename;
    tool_options[2] = shared_options->output_filename;
    tool_options[3] = shared_options->output_directory;
    
    tool_options[4] = shared_options->host_url;
    tool_options[5] = shared_options->version;
    tool_options[6] = shared_options->species;
    
    tool_options[7] = shared_options->max_batches;
    tool_options[8] = shared_options->batch_size;
    tool_options[9] = shared_options->num_threads;
    tool_options[10] = shared_options->entries_per_thread;
    
    tool_options[11] = shared_options->num_alleles;
    tool_options[12] = shared_options->coverage;
    tool_options[13] = shared_options->quality;
    tool_options[14] = shared_options->region;
    tool_options[15] = shared_options->region_file;
    tool_options[16] = shared_options->snp;
    
    tool_options[17] = shared_options->config_file;
    
    tool_options[18] = effect_options->excludes;
               
    tool_options[19] = arg_end;
    
    return tool_options;
}


int verify_effect_options(shared_options_t *global_options_data, effect_options_t *options_data) {
    // Check whether the input VCF file is defined
    if (global_options_data->vcf_filename->count == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    // Check whether the host URL is defined
    if (global_options_data->host_url->sval == NULL || strlen(global_options_data->host_url->sval) == 0) {
        LOG_ERROR("Please specify the host URL to the web service.\n");
        return EFFECT_HOST_URL_NOT_SPECIFIED;
    }

    // Check whether the version is defined
    if (global_options_data->version->sval == NULL || strlen(global_options_data->version->sval) == 0) {
        LOG_ERROR("Please specify the version.\n");
        return EFFECT_VERSION_NOT_SPECIFIED;
    }

    // Check whether the species is defined
    if (global_options_data->species->sval == NULL || strlen(global_options_data->species->sval) == 0) {
        LOG_ERROR("Please specify the species to take as reference.\n");
        return EFFECT_SPECIES_NOT_SPECIFIED;
    }

    return 0;
}
