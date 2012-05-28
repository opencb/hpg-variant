#include "gwas.h"


int read_gwas_configuration(const char *filename, gwas_options_data_t *options_data) {
    if (filename == NULL || options_data == NULL) {
        return -1;
    }
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE) {
        fprintf(stderr, "config file error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }

    // Read number of threads that will make request to the web service
    ret_code = config_lookup_int(config, "gwas.num-threads", &(options_data->num_threads));
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("num-threads config error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    LOG_INFO_F("num-threads = %ld\n", options_data->num_threads);

    // Read maximum number of batches that can be stored at certain moment
    ret_code = config_lookup_int(config, "gwas.max-batches", &(options_data->max_batches));
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("max-batches config error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    LOG_INFO_F("max-batches = %ld\n", options_data->max_batches);
    
    // Read size of every batch read
    ret_code = config_lookup_int(config, "gwas.batch-size", &(options_data->batch_size));
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("batch-size config error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    LOG_INFO_F("batch-size = %ld\n", options_data->batch_size);
    
    config_destroy(config);
    free(config);

    return 0;
}


void parse_gwas_options(int argc, char *argv[], gwas_options_data_t *options_data, global_options_data_t *global_options_data) {
    struct option *options = merge_options(gwas_options, NUM_GWAS_OPTIONS);

    char *tmp_string_field;
    int tmp_int_field;
    filter_t *filter;

    int c;
    // Last option read, for when the global options parser is invoked
    int previous_opt_index = optind;

    while ((c = getopt_long (argc, argv, "A:N:O:f:n:r:t", options, &optind)) != -1)
    {
        LOG_DEBUG_F("<main> c = %c, opt_idx = %d\n", c, optind);
        switch (c)
        {
            case 'A':
            case 'N':
            case 'O':
                optind = parse_global_options(argc, argv, global_options_data, previous_opt_index);
                break;
            case 'f':
                tmp_string_field = (char*) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(tmp_string_field, optarg);
                filter = create_region_exact_filter(tmp_string_field, 1);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                LOG_INFO_F("regions file = %s\n", optarg);
                break;
            case 'n':
                tmp_int_field = atoi(optarg);
                if (!tmp_int_field)
                {
                    LOG_WARN_F("The requested number of threads is not valid (%s)\n", optarg);
                } else
                {
                    options_data->num_threads = tmp_int_field;
                }
                LOG_INFO_F("num-threads = %ld\n", options_data->num_threads);
                break;
            case 'r':
                tmp_string_field = (char*) calloc(strlen(optarg)+1, sizeof(char));
                strcat(tmp_string_field, optarg);
                filter = create_region_exact_filter(tmp_string_field, 0);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                LOG_INFO_F("regions = %s\n", optarg);
                free(tmp_string_field);
                break;
            case 't':
                if (options_data->task == NONE) {
                    options_data->task = TDT;
                    LOG_INFO("task = TDT\n");
                } else {
                    LOG_ERROR("Task already selected\n");
                }
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


int verify_gwas_options(global_options_data_t *global_options_data, gwas_options_data_t *options_data)
{
    // Check whether the input VCF file is defined
    if (global_options_data->vcf_filename == NULL || strlen(global_options_data->vcf_filename) == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    // Check whether the input PED file is defined
    if (global_options_data->ped_filename == NULL || strlen(global_options_data->ped_filename) == 0) {
        LOG_ERROR("Please specify the input PED file.\n");
        return PED_FILE_NOT_SPECIFIED;
    }

    // Check whether the task to perform is defined
    if (options_data->task == NONE) {
        LOG_ERROR("Please specify the task to perform.\n");
        return GWAS_TASK_NOT_SPECIFIED;
    }

    
    return 0;
}
