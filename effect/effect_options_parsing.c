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
        global_options_data->host_url->sval = strdup(tmp_string);
        LOG_INFO_F("web services host URL = %s (%zu chars)\n",
                   global_options_data->host_url->sval, strlen(global_options_data->host_url->sval));
    }
    
    // Read species
    ret_code = config_lookup_string(config, "effect.species", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Species not found in configuration file, must be set via command-line");
    } else {
//         free(global_options_data->species);
//         global_options_data->species = strdup(tmp_string);
        global_options_data->species->sval = strdup(tmp_string);
        LOG_INFO_F("species = %s (%zu chars)\n",
                   global_options_data->species->sval, strlen(global_options_data->species->sval));
    }
    
    // Read version
    ret_code = config_lookup_string(config, "effect.version", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Version not found in configuration file, must be set via command-line");
    } else {
//         free(global_options_data->version);
//         global_options_data->version = strdup(tmp_string);
        global_options_data->version->sval = strdup(tmp_string);
        LOG_INFO_F("version = %s (%zu chars)\n",
                   global_options_data->version->sval, strlen(global_options_data->version->sval));
    }

    config_destroy(config);
    free(config);

    return 0;
}


void *parse_effect_options(int argc, char *argv[], effect_options_t *options_data, shared_options_t *global_options_data) {
    struct arg_end *arg_end = arg_end(options_data->num_options + global_options_data->num_options);
    void *argtable = merge_effect_options(options_data, global_options_data, arg_end);
    
    int num_errors = arg_parse(argc, argv, argtable);
    if (num_errors > 0) {
        arg_print_errors(stdout, arg_end, "hpg-variant");
    }
    
    return argtable;
}
//     struct option *options = merge_options(effect_options, NUM_EFFECT_OPTIONS);
// 
//     char *tmp_string_field;
//     int tmp_int_field;
//     filter_t *filter;
// 
//     int c;
//     // Last option read, for when the global options parser is invoked
//     int previous_opt_index = optind;
// 
//     while ((c = getopt_long (argc, argv, "A:N:O:S:U:V:Z:a:c:e:f:n:p:q:r:s:", options, &optind)) != -1) {
//         LOG_INFO_F("<main> c = %c, opt_idx = %d\n", c, optind);
//         switch (c) {
//             case 'A':
//             case 'N':
//             case 'O':
//             case 'S':
//             case 'V':
//             case 'U':
//             case 'Z':   // config, not processed
//                 optind = parse_global_options(argc, argv, global_options_data, previous_opt_index);
//                 break;
//             case 'a':
//                 if (!is_numeric(optarg)) {
//                     LOG_WARN("The argument of the filter by number of alleles must be a numeric value");
//                 } else {
//                     tmp_int_field = atoi(optarg);
//                     filter = create_num_alleles_filter(tmp_int_field);
//                     options_data->chain = add_to_filter_chain(filter, options_data->chain);
//                     LOG_INFO_F("number of alleles filter = %d\n", tmp_int_field);
//                 }
//                 break;
//             case 'c':
//                 if (!is_numeric(optarg)) {
//                     LOG_WARN("The argument of the coverage filter must be a numeric value");
//                 } else {
//                     tmp_int_field = atoi(optarg);
//                     filter = create_coverage_filter(tmp_int_field);
//                     options_data->chain = add_to_filter_chain(filter, options_data->chain);
//                     LOG_INFO_F("coverage filter, minimum is = %d\n", tmp_int_field);
//                 }
//                 break;
//             case 'e':
//                 options_data->excludes = strdup(optarg);
//                 LOG_INFO_F("exclude = %s\n", options_data->excludes);
//                 break;
//             case 'f':
//                 tmp_string_field = strdup(optarg);
//                 filter = create_region_exact_filter(tmp_string_field, 1, global_options_data->host_url, global_options_data->species, global_options_data->version);
//                 options_data->chain = add_to_filter_chain(filter, options_data->chain);
//                 LOG_INFO_F("regions file = %s\n", optarg);
//                 free(tmp_string_field);
//                 break;
//             case 'n':
//                 if (!is_numeric(optarg)) {
//                     LOG_WARN_F("The requested number of threads is not valid (%s)\n", optarg);
//                 } else {
//                     options_data->num_threads = atoi(optarg);
//                     LOG_INFO_F("num-threads = %ld\n", options_data->num_threads);
//                 }
//                 break;
//             case 'p':
//                 if (!is_numeric(optarg)) {
//                     LOG_WARN_F("The requested number of variants per request is not valid (%s)\n", optarg);
//                 } else {
//                     options_data->variants_per_request = atoi(optarg);
//                     LOG_INFO_F("variants-per-request = %ld\n", options_data->variants_per_request);
//                 }
//                 break;
//             case 'q':
//                 if (!is_numeric(optarg)) {
//                     LOG_WARN("The argument of the quality filter must be a numeric value");
//                 } else {
//                     tmp_int_field = atoi(optarg);
//                     filter = create_quality_filter(tmp_int_field);
//                     options_data->chain = add_to_filter_chain(filter, options_data->chain);
//                     LOG_INFO_F("quality filter, minimum is = %d\n", tmp_int_field);
//                 }
//                 break;
//             case 'r':
//                 tmp_string_field = strdup(optarg);
//                 filter = create_region_exact_filter(tmp_string_field, 0, global_options_data->host_url, global_options_data->species, global_options_data->version);
//                 options_data->chain = add_to_filter_chain(filter, options_data->chain);
//                 LOG_INFO_F("regions = %s\n", optarg);
//                 free(tmp_string_field);
//                 break;
//             case 's':
//                 filter = create_snp_filter(optarg);
//                 options_data->chain = add_to_filter_chain(filter, options_data->chain);
//                 LOG_INFO_F("snp filter to %s SNPs\n", (optarg == NULL)? "include" : optarg);
//                 break;
//             case '?':
//             default:
//                 LOG_WARN("Option unknown\n");
//                 break;
//         }
//         
//         printf("previous opt index = %d\n", optind);
//         previous_opt_index = optind;
//         
//         printf("argc = %d\targv = %d\toptions = %d\toptind = %d\n", argc, argv != NULL, options != NULL, optind);
//     }
// 
//     free(options);
// }

void* merge_effect_options(effect_options_t *options_data, shared_options_t *global_options_data, struct arg_end *arg_end) {
    size_t opts_size = options_data->num_options + global_options_data->num_options + 1;
    void *tool_options = malloc (opts_size * sizeof(void*));
    
    // Add global options
    for (int i = 0; i < global_options_data->num_options; i++) {
        tool_options[i] = global_options_data[i];
    }
    
    // Add local options
    for (int i = 0; i < options_data->num_options; i++) {
        tool_options[global_options_data->num_options + i] = options_data[i];
    }

    tool_options[global_options_data->num_options + options_data->num_options] = arg_end;
    
    return tool_options;
}


int verify_effect_options(shared_options_t *global_options_data, effect_options_t *options_data) {
    // Check whether the input VCF file is defined
    if (global_options_data->vcf_filename->count == 0) {// || strlen(global_options_data->vcf_filename->filename) == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    // Check whether the host URL is defined
    if (global_options_data->host_url->count == 0) {// || strlen(global_options_data->host_url->sval) == 0) {
        LOG_ERROR("Please specify the host URL to the web service.\n");
        return EFFECT_HOST_URL_NOT_SPECIFIED;
    }

    // Check whether the version is defined
    if (global_options_data->version->count == 0) {// || strlen(global_options_data->version->sval) == 0) {
        LOG_ERROR("Please specify the version.\n");
        return EFFECT_VERSION_NOT_SPECIFIED;
    }

    // Check whether the species is defined
    if (global_options_data->species->count == 0) {// || strlen(global_options_data->species->sval) == 0) {
        LOG_ERROR("Please specify the species to take as reference.\n");
        return EFFECT_SPECIES_NOT_SPECIFIED;
    }

    return 0;
}
