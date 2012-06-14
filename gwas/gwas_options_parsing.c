#include "gwas.h"


int read_gwas_configuration(const char *filename, gwas_options_data_t *options_data) {
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
    ret_code = config_lookup_int(config, "gwas.num-threads", &(options_data->num_threads));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of threads not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("num-threads = %ld\n", options_data->num_threads);
    }

    // Read maximum number of batches that can be stored at certain moment
    ret_code = config_lookup_int(config, "gwas.max-batches", &(options_data->max_batches));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Maximum number of batches not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("max-batches = %ld\n", options_data->max_batches);
    }
    
    // Read size of every batch read
    ret_code = config_lookup_int(config, "gwas.batch-size", &(options_data->batch_size));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Batch size not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("batch-size = %ld\n", options_data->batch_size);
    }
    
    // Read number of variants per request of test execution
    ret_code = config_lookup_int(config, "gwas.variants-per-request", &(options_data->variants_per_request));
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Variants per request not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("variants-per-request = %ld\n", options_data->variants_per_request);
    }

    config_destroy(config);
    free(config);

    return 0;
}


void parse_gwas_options(int argc, char *argv[], gwas_options_data_t *options_data, shared_options_data_t *global_options_data) {
    struct option *options = merge_options(gwas_options, NUM_GWAS_OPTIONS);

    char *tmp_string_field;
    int tmp_int_field;
    filter_t *filter;

    int c;
    // Last option read, for when the global options parser is invoked
    int previous_opt_index = optind;

    while ((c = getopt_long (argc, argv, "A:N:O:S:U:V:a:c:f:n:o:q:r:s:t", options, &optind)) != -1)
    {
        LOG_DEBUG_F("<main> c = %c, opt_idx = %d\n", c, optind);
        switch (c)
        {
            case 'A':
            case 'N':
            case 'O':
            case 'S':
            case 'V':
            case 'U':
                optind = parse_global_options(argc, argv, global_options_data, previous_opt_index);
                break;
            case 'a':
                if (!is_numeric(optarg)) {
                    LOG_WARN("The argument of the filter by number of alleles must be a numeric value");
                } else {
                    tmp_int_field = atoi(optarg);
                    filter = create_num_alleles_filter(tmp_int_field);
                    options_data->chain = add_to_filter_chain(filter, options_data->chain);
                    LOG_INFO_F("number of alleles filter = %d\n", tmp_int_field);
                }
                break;
            case 'c':
                if (!is_numeric(optarg)) {
                    LOG_WARN("The argument of the coverage filter must be a numeric value");
                } else {
                    tmp_int_field = atoi(optarg);
                    filter = create_coverage_filter(tmp_int_field);
                    options_data->chain = add_to_filter_chain(filter, options_data->chain);
                    LOG_INFO_F("coverage filter, minimum is = %d\n", tmp_int_field);
                }
                break;
            case 'f':
                tmp_string_field = strdup(optarg);
                filter = create_region_exact_filter(tmp_string_field, 1, global_options_data->host_url, global_options_data->species, global_options_data->version);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                LOG_INFO_F("regions file = %s\n", optarg);
                free(tmp_string_field);
                break;
            case 'i':
                if (options_data->task == NONE) {
                    options_data->task = FISHER;
                    LOG_INFO("task = Fisher's exact test\n");
                } else {
                    LOG_ERROR("Task already selected\n");
                }
                break;
            case 'n':
                if (!is_numeric(optarg)) {
                    LOG_WARN_F("The requested number of threads is not valid (%s)\n", optarg);
                } else {
                    options_data->num_threads = atoi(optarg);
                    LOG_INFO_F("num-threads = %ld\n", options_data->num_threads);
                }
                break;
            case 'o':
                if (options_data->task == NONE) {
                    options_data->task = ASSOCIATION_BASIC;
                    LOG_INFO("task = Basic case/control association\n");
                } else {
                    LOG_ERROR("Task already selected\n");
                }
                break;
            case 'q':
                if (!is_numeric(optarg)) {
                    LOG_WARN("The argument of the quality filter must be a numeric value");
                } else {
                    tmp_int_field = atoi(optarg);
                    filter = create_quality_filter(tmp_int_field);
                    options_data->chain = add_to_filter_chain(filter, options_data->chain);
                    LOG_INFO_F("quality filter, minimum is = %d\n", tmp_int_field);
                }
                break;
            case 'r':
                tmp_string_field = strdup(optarg);
                filter = create_region_exact_filter(tmp_string_field, 0, global_options_data->host_url, global_options_data->species, global_options_data->version);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                LOG_INFO_F("regions = %s\n", optarg);
                free(tmp_string_field);
                break;
            case 's':
                filter = create_snp_filter(optarg);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                LOG_INFO_F("snp filter to %s SNPs\n", (optarg == NULL)? "include" : optarg);
                break;
            case 't':
                if (options_data->task == NONE) {
                    options_data->task = TDT;
                    LOG_INFO("task = TDT\n");
                } else {
                    LOG_ERROR("Task already chosen\n");
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


int verify_gwas_options(shared_options_data_t *global_options_data, gwas_options_data_t *options_data)
{
    // Check whether the input VCF file is defined
    if (global_options_data->vcf_filename == NULL || strlen(global_options_data->vcf_filename) == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    // Check whether the task to perform is defined
    if (options_data->task == NONE) {
        LOG_ERROR("Please specify the task to perform.\n");
        return GWAS_TASK_NOT_SPECIFIED;
    }

    // Check whether the input PED file is defined
    if (global_options_data->ped_filename == NULL || strlen(global_options_data->ped_filename) == 0) {
        LOG_ERROR("Please specify the input PED file.\n");
        return PED_FILE_NOT_SPECIFIED;
    }
    
    return 0;
}
