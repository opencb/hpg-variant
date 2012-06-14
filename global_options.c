#include "global_options.h"

shared_options_t *new_global_options(void) {
    shared_options_t *options_data = (shared_options_t*) calloc (1, sizeof(shared_options_t));
    
    options_data->vcf_filename = arg_file1(NULL, "vcf-file", NULL, "VCF file used as input");
    options_data->ped_filename = arg_file0(NULL, "ped-file", NULL, "PED file used as input");
    options_data->output_filename = arg_file0(NULL, "out", NULL, "Filename prefix for main output files");
    options_data->output_directory = arg_str0(NULL, "outdir", NULL, "Directory where the output files will be stored");
    
    options_data->host_url = arg_str0(NULL, "url", NULL, "URL of the host where remote web services run");
    options_data->version = arg_str0(NULL, "version", NULL, "Version of the web service to query");
    options_data->species = arg_str0(NULL, "species", NULL, "Species whose genome is taken as reference");
    
    options_data->max_batches = arg_int0(NULL, "num-batches", NULL, "Maximum number of batches stored at the same time");
    options_data->batch_size = arg_int0(NULL, "batch-size", NULL, "Maximum size of a batch");
    options_data->num_threads = arg_int0(NULL, "num-threads", NULL, "Number of threads when a task is perform in parallel");
    options_data->entries_per_thread = arg_int0(NULL, "entries-per-thread", NULL, "Number of entries in a batch each thread processes");
    
    options_data->num_alleles = arg_int0(NULL, "alleles", NULL, "Filter: by number of alleles");
    options_data->coverage = arg_int0(NULL, "coverage", NULL, "Filter: by minimum coverage");
    options_data->quality = arg_int0(NULL, "quality", NULL, "Filter: by minimum quality");
    options_data->region = arg_str0(NULL, "region", NULL, "Filter: by a list of regions (chr1:start1-end1,chr2:start2-end2...)");
    options_data->region_file = arg_file0(NULL, "region-file", NULL, "Filter: by a list of regions (read from a GFF file)");
    options_data->snp = arg_str0(NULL, "snp", NULL, "Filter: by being a SNP or not");
    
    options_data->num_options = NUM_GLOBAL_OPTIONS;
//     options_data->output_directory = (char*) calloc (6, sizeof(char));
//     strncat(options_data->output_directory, "/tmp/", 5);
//     options_data->ped_filename = NULL;
//     options_data->vcf_filename = NULL;
//     options_data->output_filename = NULL;
    
    return options_data;
}

shared_options_data_t* new_global_options_data(shared_options_t* options) {
    shared_options_data_t *options_data = (shared_options_data_t*) calloc (1, sizeof(shared_options_data_t));
    
    options_data->vcf_filename = *(options->vcf_filename->filename);
    options_data->ped_filename = *(options->ped_filename->filename);
    options_data->output_filename = *(options->output_filename->filename);
    options_data->output_directory = *(options->output_directory->sval);
    
    options_data->host_url = *(options->host_url->sval);
    options_data->version = *(options->version->sval);
    options_data->species = *(options->species->sval);
    
    options_data->max_batches = *(options->max_batches->ival);
    options_data->batch_size = *(options->batch_size->ival);
    options_data->num_threads = *(options->num_threads->ival);
    options_data->entries_per_thread = *(options->entries_per_thread->ival);
    
    options_data->num_alleles = *(options->num_alleles->ival);
    options_data->coverage = *(options->coverage->ival);
    options_data->quality = *(options->quality->ival);
    options_data->region = options->region->count > 0 ? *(options->region->sval) 
                                                      : ( options->region_file->count > 0 ? *(options->region_file->filename) : NULL );
    options_data->snp = *(options->snp->sval);
    
    return options_data;
}


void free_global_options_data(shared_options_data_t *options_data) {
//     if (options_data->ped_filename)     { free((void*) options_data->ped_filename); }
    if (options_data->vcf_filename)     { free((void*) options_data->vcf_filename); }
    if (options_data->output_directory) { free((void*) options_data->output_directory); }
    if (options_data->output_filename)  { free((void*) options_data->output_filename); }
    free(options_data);
}

int read_global_configuration(const char *filename, shared_options_t *options_data) {
    if (filename == NULL || options_data == NULL) {
        return -1;
    }
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    
    // Read output directory
    const char *tmp_string;
    ret_code = config_lookup_string(config, "global.outdir", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Output folder not found in configuration file, must be set via command-line");
    } else {
//         free(options_data->output_directory);
//         options_data->output_directory = strdup(tmp_string);
        options_data->output_directory->sval = strdup(tmp_string);
        LOG_INFO_F("Output folder = %s (%zu chars)\n",
                   *(options_data->output_directory->sval), 
                   strlen(*(options_data->output_directory->sval)));
    }
    
//     // Read host URL
//     ret_code = config_lookup_string(config, "global.url", &tmp_string);
//     if (ret_code == CONFIG_FALSE) {
//         LOG_WARN("Web services URL not found in config file, must be set via command-line");
//     } else {
//         free(options_data->host_url);
//         options_data->host_url = strdup(tmp_string);
//         LOG_INFO_F("web services host URL = %s (%zu chars)\n",
//                    options_data->host_url, strlen(options_data->host_url));
//     }
//     
//     // Read species
//     ret_code = config_lookup_string(config, "global.species", &tmp_string);
//     if (ret_code == CONFIG_FALSE) {
//         LOG_WARN("Species not found in config file, must be set via command-line");
//     } else {
//         free(options_data->species);
//         options_data->species = strdup(tmp_string);
//         LOG_INFO_F("species = %s (%zu chars)\n",
//                    options_data->species, strlen(options_data->species));
//     }
//     
//     // Read version
//     ret_code = config_lookup_string(config, "global.version", &tmp_string);
//     if (ret_code == CONFIG_FALSE) {
//         LOG_WARN("Version not found in config file, must be set via command-line");
//     } else {
//         free(options_data->version);
//         options_data->version = strdup(tmp_string);
//         LOG_INFO_F("version = %s (%zu chars)\n",
//                    options_data->version, strlen(options_data->version));
//     }
    
    config_destroy(config);
    free(config);
    
    return ret_code;
}

// int parse_global_options(int argc, char *argv[], shared_options_t *options_data, int start_index) {
//     int c;
//     int finished = 0;
//     int last_opt_index = optind = start_index;
//     
//     LOG_DEBUG("Parsing global options...\n");
//     while (!finished && (c = getopt_long (argc, argv, "A:E:N:O:S:V:U:Z:", global_options, &optind)) != 1) {
//         switch (c) {
//             case 'A':
//                 options_data->vcf_filename = strdup(optarg);
//                 last_opt_index = optind;
//                 LOG_INFO_F("input vcf filename = %s\n", options_data->vcf_filename);
//                 break;
//             case 'E':
//                 options_data->ped_filename = strdup(optarg);
//                 last_opt_index = optind;
//                 LOG_INFO_F("input ped filename = %s\n", options_data->ped_filename);
//                 break;
//             case 'N':
//                 options_data->output_directory = strdup(optarg);
//                 last_opt_index = optind;
//                 LOG_INFO_F("output directory = %s\n", options_data->output_directory);
//                 break;
//             case 'O':
//                 options_data->output_filename = strdup(optarg);
//                 last_opt_index = optind;
//                 LOG_INFO_F("output filename template = %s\n", options_data->output_filename);
//                 break;
//             case 'S':
//                 options_data->species = strdup(optarg);
//                 last_opt_index = optind;
//                 LOG_INFO_F("species = %s\n", options_data->species);
//                 break;
//             case 'V':
//                 options_data->version = strdup(optarg);
//                 last_opt_index = optind;
//                 LOG_INFO_F("version = %s\n", options_data->version);
//                 break;
//             case 'U':
//                 options_data->host_url = strdup(optarg);
//                 last_opt_index = optind;
//                 LOG_INFO_F("web services URL = %s\n", options_data->host_url);
//                 break;
//             case 'Z':
//                 last_opt_index = optind;
//                 LOG_DEBUG("--config found but not processed");
//                 break;
//             case '?':
//                 LOG_WARN("undefined common option\n");
//             default:
//                 finished = 1;
//                 LOG_INFO("global opts - default\n");
//         }
//     }
//     
//     LOG_INFO_F("last opt index = %d\n", last_opt_index);
//     return last_opt_index;
// }

int verify_global_options(shared_options_t *options_data) {
    // Nothing to verify yet
    return 0;
}

// void* merge_options(effect_options_data_t *options_data, global_options_data_t *global_options_data, size_t num_local_options) {
//     size_t opts_size = num_local_options + NUM_GLOBAL_OPTIONS + 1;
//     void *tool_options = malloc (opts_size * sizeof(void*));
//     
//     // Add global options
//     for (int i = 0; i < NUM_GLOBAL_OPTIONS; i++) {
//         tool_options[i] = global_options_data[i];
//     }
//     
//     // Add local options
//     for (int i = 0; i < num_local_options; i++) {
//         tool_options[NUM_GLOBAL_OPTIONS + i] = options_data[i];
//     }
// 
//     return tool_options;
// }


// struct option *merge_options(struct option local_options[], size_t num_local_options) {
//     size_t opts_size = num_local_options + NUM_GLOBAL_OPTIONS + 1;
//     struct option *tool_options = (struct option*) malloc (opts_size * sizeof(struct option));
//     
//     // Add global options
//     for (int i = 0; i < NUM_GLOBAL_OPTIONS; i++) {
//         tool_options[i] = global_options[i];
//     }
//     
//     // Add local options
//     for (int i = 0; i < num_local_options; i++) {
//         tool_options[NUM_GLOBAL_OPTIONS + i] = local_options[i];
//     }
//     
//     return tool_options;
// }
