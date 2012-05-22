#include "filter.h"


int read_filter_configuration(const char *filename, filter_options_data_t *options_data) {
    config_t *config = (config_t*) malloc (sizeof(config_t));
    
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("config file error: %s\n", config_error_text(config));
        return ret_code;
    }
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "filter.num-threads", &(options_data->num_threads));
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("num-threads config error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    LOG_INFO_F("num-threads = %ld\n", options_data->num_threads);
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "filter.max-batches", &(options_data->max_batches));
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("max-batches config error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    LOG_INFO_F("max-batches = %ld\n", options_data->max_batches);
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "filter.batch-size", &(options_data->batch_size));
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("batch-size config error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    LOG_INFO_F("batch-size = %ld\n", options_data->batch_size);
    
    return 0;
}


void parse_filter_options(int argc, char *argv[], filter_options_data_t *options_data, global_options_data_t *global_options_data) {
    const struct option *options = merge_options(filter_options, NUM_FILTER_OPTIONS);

    char *tmp_string_field;
    int tmp_int_field;
    filter_t *filter;

    int c;
    // Last option read, for when the global options parser is invoked
    int previous_opt_index = optind;

    while ((c = getopt_long (argc, argv, "A:N:O:c:f:q:r:s:", options, &optind)) != -1) {
        LOG_DEBUG_F("<main> c = %c, opt_idx = %d\n", c, optind);
        switch (c) {
            case 'A':
            case 'N':
            case 'O':
                optind = parse_global_options(argc, argv, global_options_data, previous_opt_index);
                break;
            case 'c':
                if (!is_numeric(optarg)) {
                    LOG_WARN("Coverage filter argument must be a numeric value");
                } else {
                    tmp_int_field = atoi(optarg);
                    filter = create_coverage_filter(tmp_int_field);
                    options_data->chain = add_to_filter_chain(filter, options_data->chain);
                    LOG_INFO_F("coverage filter, minimum is = %d\n", tmp_int_field);
                }
                break;
            case 'f':
                tmp_string_field = (char*) calloc(strlen(optarg)+1, sizeof(char));
                strncat(tmp_string_field, optarg, strlen(optarg));
                filter = create_region_filter(tmp_string_field, 1);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                LOG_INFO_F("regions file = %s\n", optarg);
                break;
            case 'q':
                if (!is_numeric(optarg)) {
                    LOG_WARN("Quality filter argument must be a numeric value");
                } else {
                    tmp_int_field = atoi(optarg);
                    filter = create_quality_filter(tmp_int_field);
                    options_data->chain = add_to_filter_chain(filter, options_data->chain);
                    LOG_INFO_F("quality filter, minimum is = %d\n", tmp_int_field);
                }
                break;
            case 'r':
                tmp_string_field = (char*) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(tmp_string_field, optarg);
                filter = create_region_filter(tmp_string_field, 0);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                LOG_INFO_F("regions = %s\n", optarg);
            case 's':
                filter = create_snp_filter(optarg);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                LOG_INFO_F("snp filter to %s SNPs\n", (optarg == NULL)? "include" : optarg);
                break;
            case '?':
            default:
                LOG_WARN("Option unknown\n");
                break;
        }
        
        previous_opt_index = optind;
    }
}


int verify_filter_options(global_options_data_t *global_options_data, filter_options_data_t *options_data)
{
    // Check whether the input VCF file is defined
    if (global_options_data->vcf_filename == NULL || strlen(global_options_data->vcf_filename) == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    // Check whether a filter or more has been specified
    if (cp_heap_count(options_data->chain) == 0) {
        LOG_ERROR("Please specify at least one filter\n");
        return EMPTY_LIST_OF_FILTERS;
    }

    return 0;
}

