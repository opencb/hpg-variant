#include "stats.h"


int read_stats_configuration(const char *filename, stats_options_data_t *options_data) {
    config_t *config = (config_t*) malloc (sizeof(config_t));
    
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("config file error: %s\n", config_error_text(config));
        return ret_code;
    }
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "stats.max-batches", &(options_data->max_batches));
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("max-batches config error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    LOG_INFO_F("max-batches = %ld\n", options_data->max_batches);
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "stats.batch-size", &(options_data->batch_size));
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("batch-size config error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    LOG_INFO_F("batch-size = %ld\n", options_data->batch_size);
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "stats.num-threads", &(options_data->num_threads));
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("num-threads config error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    LOG_INFO_F("num-threads = %ld\n", options_data->num_threads);
    
    // Read number of variants each thread will handle
    ret_code = config_lookup_int(config, "stats.variants-per-threads", &(options_data->variants_per_thread));
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("variants-per-threads config error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    LOG_INFO_F("variants-per-threads = %ld\n", options_data->variants_per_thread);
    
    config_destroy(config);
    free(config);

    return 0;
}


void parse_stats_options(int argc, char *argv[], stats_options_data_t *options_data, global_options_data_t *global_options_data) {
    const struct option *options = merge_options(stats_options, NUM_STATS_OPTIONS);

    char *tmp_string_field;
    int tmp_int_field;
    
    int c;
    // Last option read, for when the global options parser is invoked
    int previous_opt_index = optind;
    
    int debug = 1;
    while ((c = getopt_long (argc, argv, "A:N:O:", options, &optind)) != -1) {
        LOG_DEBUG_F("<main> c = %c, opt_idx = %d\n", c, optind);
        switch (c) {
            case 'A':
            case 'N':
            case 'O':
                optind = parse_global_options(argc, argv, global_options_data, previous_opt_index);
                break;
            case '?':
            default:
                LOG_WARN("Option unknown\n");
                break;
        }
        
        previous_opt_index = optind;
    }
}


int verify_stats_options(global_options_data_t *global_options_data, stats_options_data_t *options_data)
{
    // Check whether the input VCF file is defined
    if (global_options_data->vcf_filename == NULL || strlen(global_options_data->vcf_filename) == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    return 0;
}

