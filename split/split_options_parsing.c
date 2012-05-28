#include "split.h"


int read_split_configuration(const char *filename, split_options_data_t *options_data) {
    config_t *config = (config_t*) malloc (sizeof(config_t));
    
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("config file error: %s\n", config_error_text(config));
        return ret_code;
    }
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "split.max-batches", &(options_data->max_batches));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Maximum number of batches not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("max-batches = %ld\n", options_data->max_batches);
    }
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "split.batch-size", &(options_data->batch_size));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Batch size not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("batch-size = %ld\n", options_data->batch_size);
    }
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "split.num-threads", &(options_data->num_threads));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of threads not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("num-threads = %ld\n", options_data->num_threads);
    }
    
    // Read number of variants each thread will handle
    ret_code = config_lookup_int(config, "split.variants-per-threads", &(options_data->variants_per_thread));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of variants per thread not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("variants-per-threads = %ld\n", options_data->variants_per_thread);
    }
    
    config_destroy(config);
    free(config);

    return 0;
}


void parse_split_options(int argc, char *argv[], split_options_data_t *options_data, global_options_data_t *global_options_data) {
    const struct option *options = merge_options(split_options, NUM_SPLIT_OPTIONS);

    char *tmp_string_field;
    int tmp_int_field;
    
    int c;
    // Last option read, for when the global options parser is invoked
    int previous_opt_index = optind;
    
    while ((c = getopt_long (argc, argv, "A:N:O:", options, &optind)) != -1) {
        LOG_DEBUG_F("<main> c = %c, opt_idx = %d\n", c, optind);
        switch (c) {
            case 'A':
            case 'N':
            case 'O':
                optind = parse_global_options(argc, argv, global_options_data, previous_opt_index);
                break;
            case 'c':
                tmp_string_field = (char*) calloc(strlen(optarg)+1, sizeof(char));
                strcpy(tmp_string_field, optarg);
                if (!strcasecmp("chromosome", tmp_string_field)) {
                    options_data->criterion = CHROMOSOME;
                } else if (!strcasecmp("gene", tmp_string_field)) {
                    options_data->criterion = GENE;
                    LOG_FATAL("Gene criterion not implemented yet!");
                }
                LOG_INFO_F("Splitting by %s\n", tmp_string_field);
                break;
            case '?':
            default:
                LOG_WARN("Option unknown\n");
                break;
        }
        
        previous_opt_index = optind;
    }
}


int verify_split_options(global_options_data_t *global_options_data, split_options_data_t *options_data) {
    // Check whether the input VCF file is defined
    if (global_options_data->vcf_filename == NULL || strlen(global_options_data->vcf_filename) == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    if (options_data->criterion == NONE) {
        LOG_ERROR("Please specify a splittering criterion.\n");
        return NONE_CRITERION_SPECIFIED;
    }
    
    return 0;
}

