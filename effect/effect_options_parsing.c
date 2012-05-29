#include "effect.h"


int read_effect_configuration(const char *filename, effect_options_data_t *options_data) {
    if (filename == NULL || options_data == NULL) {
        return -1;
    }
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }

    // Read number of threads that will make request to the web service
    ret_code = config_lookup_int(config, "effect.num-threads", &(options_data->num_threads));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of threads not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("num-threads = %ld\n", options_data->num_threads);
    }

    // Read maximum number of batches that can be stored at certain moment
    ret_code = config_lookup_int(config, "effect.max-batches", &(options_data->max_batches));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Maximum number of batches not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("max-batches = %ld\n", options_data->max_batches);
    }
    
    // Read size of every batch read
    ret_code = config_lookup_int(config, "effect.batch-size", &(options_data->batch_size));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Batch size not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("batch-size = %ld\n", options_data->batch_size);
    }
    
    // Read number of variants per request to the web service
    ret_code = config_lookup_int(config, "effect.variants-per-request", &(options_data->variants_per_request));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Variants per request not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("variants-per-request = %ld\n", options_data->variants_per_request);
    }

    config_destroy(config);
    free(config);

    return 0;
}


void parse_effect_options(int argc, char *argv[], effect_options_data_t *options_data, global_options_data_t *global_options_data) {
    struct option *options = merge_options(effect_options, NUM_EFFECT_OPTIONS);

    char *tmp_string_field;
    int tmp_int_field;
    filter_t *filter;

    int c;
    // Last option read, for when the global options parser is invoked
    int previous_opt_index = optind;

    while ((c = getopt_long (argc, argv, "A:N:O:S:U:V:Z:e:f:n:p:r:", options, &optind)) != -1) {
        LOG_DEBUG_F("<main> c = %c, opt_idx = %d\n", c, optind);
        switch (c) {
            case 'A':
            case 'N':
            case 'O':
            case 'S':
            case 'V':
            case 'U':
            case 'Z':   // config, not processed
                optind = parse_global_options(argc, argv, global_options_data, previous_opt_index);
                break;
            case 'e':
                options_data->excludes = (char*) calloc(strlen(optarg)+1, sizeof(char));
                strcat(options_data->excludes, optarg);
                LOG_INFO_F("exclude = %s\n", options_data->excludes);
                break;
            case 'f':
                tmp_string_field = strdup(optarg);
                filter = create_region_exact_filter(tmp_string_field, 1, global_options_data->host_url, global_options_data->species, global_options_data->version);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                LOG_INFO_F("regions file = %s\n", optarg);
                free(tmp_string_field);
                break;
            case 'n':
                tmp_int_field = atoi(optarg);
                if (!tmp_int_field) {
                    LOG_WARN_F("The requested number of threads is not valid (%s)\n", optarg);
                } else {
                    options_data->num_threads = tmp_int_field;
                }
                LOG_INFO_F("num-threads = %ld\n", options_data->num_threads);
                break;
            case 'p':
                tmp_int_field = atoi(optarg);
                if (!tmp_int_field) {
                    LOG_WARN_F("The requested number of variants per request is not valid (%s)\n", optarg);
                } else {
                    options_data->variants_per_request = tmp_int_field;
                }
                LOG_INFO_F("variants-per-request = %ld\n", options_data->variants_per_request);
                break;
            case 'r':
                tmp_string_field = strdup(optarg);
                filter = create_region_exact_filter(tmp_string_field, 0, global_options_data->host_url, global_options_data->species, global_options_data->version);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                LOG_INFO_F("regions = %s\n", optarg);
                free(tmp_string_field);
                break;
            case '?':
            default:
                LOG_WARN("Option unknown\n");
                break;
        }
        
        previous_opt_index = optind;
    }

    free(options);
}


int verify_effect_options(global_options_data_t *global_options_data, effect_options_data_t *options_data) {
    // Check whether the input VCF file is defined
    if (global_options_data->vcf_filename == NULL || strlen(global_options_data->vcf_filename) == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    // Check whether the host URL is defined
    if (global_options_data->host_url == NULL || strlen(global_options_data->host_url) == 0) {
        LOG_ERROR("Please specify the host URL to the web service.\n");
        return EFFECT_HOST_URL_NOT_SPECIFIED;
    }

    // Check whether the version is defined
    if (global_options_data->version == NULL || strlen(global_options_data->version) == 0) {
        LOG_ERROR("Please specify the version.\n");
        return EFFECT_VERSION_NOT_SPECIFIED;
    }

    // Check whether the species is defined
    if (global_options_data->species == NULL || strlen(global_options_data->species) == 0) {
        LOG_ERROR("Please specify the species to take as reference.\n");
        return EFFECT_SPECIES_NOT_SPECIFIED;
    }

    return 0;
}
