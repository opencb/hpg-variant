#include "stats.h"


int read_stats_configuration(const char *filename, stats_options_t *options, shared_options_t *shared_options) {
    if (filename == NULL || options == NULL || shared_options == NULL) {
        return -1;
    }
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    printf("** 0.5\n");
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("config file error: %s\n", config_error_text(config));
        return ret_code;
    }
    
    // Read number of threads that will make request to the web service
    ret_code = config_lookup_int(config, "stats.num-threads", shared_options->num_threads->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of threads not found in config file, must be set via command-line");
    } else {
        LOG_INFO_F("num-threads = %ld\n", *(shared_options->num_threads->ival));
    }

    // Read maximum number of batches that can be stored at certain moment
    ret_code = config_lookup_int(config, "stats.max-batches", shared_options->max_batches->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Maximum number of batches not found in configuration file, must be set via command-line");
    } else {
        LOG_INFO_F("max-batches = %ld\n", *(shared_options->max_batches->ival));
    }
    
    // Read size of every batch read
    ret_code = config_lookup_int(config, "stats.batch-size", shared_options->batch_size->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Batch size not found in configuration file, must be set via command-line");
    } else {
        LOG_INFO_F("batch-size = %ld\n", *(shared_options->batch_size->ival));
    }
    
    // Read number of variants per request to the web service
    ret_code = config_lookup_int(config, "stats.entries-per-thread", shared_options->entries_per_thread->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Variants per request not found in configuration file, must be set via command-line");
    } else {
        LOG_INFO_F("entries-per-thread = %ld\n", *(shared_options->entries_per_thread->ival));
    }
    
    config_destroy(config);
    free(config);

    return 0;
}

void **parse_stats_options(int argc, char *argv[], stats_options_t *stats_options, shared_options_t *shared_options) {
    struct arg_end *end = arg_end(stats_options->num_options + shared_options->num_options);
    void **argtable = merge_stats_options(stats_options, shared_options, end);
    
    int num_errors = arg_parse(argc, argv, argtable);
    if (num_errors > 0) {
        arg_print_errors(stdout, end, "hpg-variant");
    }
    
    return argtable;
}

void **merge_stats_options(stats_options_t *stats_options, shared_options_t *shared_options, struct arg_end *arg_end) {
    size_t opts_size = stats_options->num_options + shared_options->num_options + 1;
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
    
    tool_options[11] = shared_options->config_file;
    tool_options[12] = shared_options->mmap_vcf_files;
    
    tool_options[13] = arg_end;
    
    return tool_options;
}


int verify_stats_options(stats_options_t *stats_options, shared_options_t *shared_options) {
    // Check whether the input VCF file is defined
    if (shared_options->vcf_filename->count == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    return 0;
}


// void parse_stats_options(int argc, char *argv[], stats_options_data_t *options_data, global_options_data_t *global_options_data) {
//     const struct option *options = merge_options(stats_options, NUM_STATS_OPTIONS);
// 
//     char *tmp_string_field;
//     int tmp_int_field;
//     
//     int c;
//     // Last option read, for when the global options parser is invoked
//     int previous_opt_index = optind;
//     
//     while ((c = getopt_long (argc, argv, "A:N:O:", options, &optind)) != -1) {
//         LOG_DEBUG_F("<main> c = %c, opt_idx = %d\n", c, optind);
//         switch (c) {
//             case 'A':
//             case 'N':
//             case 'O':
//                 optind = parse_global_options(argc, argv, global_options_data, previous_opt_index);
//                 break;
//             case '?':
//             default:
//                 LOG_WARN("Option unknown\n");
//                 break;
//         }
//         
//         previous_opt_index = optind;
//     }
// }
// 
// 
// int verify_stats_options(global_options_data_t *global_options_data, stats_options_data_t *options_data)
// {
//     // Check whether the input VCF file is defined
//     if (global_options_data->vcf_filename == NULL || strlen(global_options_data->vcf_filename) == 0) {
//         LOG_ERROR("Please specify the input VCF file.\n");
//         return VCF_FILE_NOT_SPECIFIED;
//     }
//     
//     return 0;
// }

