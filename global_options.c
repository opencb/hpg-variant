#include "global_options.h"

global_options_data_t *new_global_options_data(void)
{
    global_options_data_t *options_data = (global_options_data_t*) calloc (1, sizeof(global_options_data_t));
    
//     options_data->output_directory = (char*) calloc (6, sizeof(char));
//     strncat(options_data->output_directory, "/tmp/", 5);
//     options_data->ped_filename = NULL;
//     options_data->vcf_filename = NULL;
//     options_data->output_filename = NULL;
    
    return options_data;
}

void free_global_options_data(global_options_data_t *options_data)
{
//     if (options_data->ped_filename)     { free((void*) options_data->ped_filename); }
    if (options_data->vcf_filename)     { free((void*) options_data->vcf_filename); }
    if (options_data->output_directory) { free((void*) options_data->output_directory); }
    if (options_data->output_filename)  { free((void*) options_data->output_filename); }
    free(options_data);
}

int read_global_configuration(const char *filename, global_options_data_t *options_data)
{
    if (filename == NULL || options_data == NULL) {
        return -1;
    }
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE)
    {
        LOG_ERROR_F("config file error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }
    
    // Read output directory
    const char *tmp_string;
    ret_code = config_lookup_string(config, "global.outdir", &tmp_string);
    if (ret_code == CONFIG_FALSE)
    {
        LOG_WARN("Output directory not found in config file, must be set via command-line");
    } else {
        free(options_data->output_directory);
        options_data->output_directory = strdup(tmp_string);
        LOG_INFO_F("output directory = %s (%zu chars)\n", options_data->output_directory, 
                    strlen(options_data->output_directory));
    }
    
    // Read host URL
    ret_code = config_lookup_string(config, "global.url", &tmp_string);
    if (ret_code == CONFIG_FALSE)
    {
        LOG_WARN("Web services URL not found in config file, must be set via command-line");
    } else {
        free(options_data->host_url);
        options_data->host_url = strdup(tmp_string);
        LOG_INFO_F("web services host URL = %s (%zu chars)\n",
                   options_data->host_url, strlen(options_data->host_url));
    }
    
    // Read species
    ret_code = config_lookup_string(config, "global.species", &tmp_string);
    if (ret_code == CONFIG_FALSE)
    {
        LOG_WARN("Species not found in config file, must be set via command-line");
    } else {
        free(options_data->species);
        options_data->species = strdup(tmp_string);
        LOG_INFO_F("species = %s (%zu chars)\n",
                   options_data->species, strlen(options_data->species));
    }
    
    // Read version
    ret_code = config_lookup_string(config, "global.version", &tmp_string);
    if (ret_code == CONFIG_FALSE)
    {
        LOG_WARN("Version not found in config file, must be set via command-line");
    } else {
        free(options_data->version);
        options_data->version = strdup(tmp_string);
        LOG_INFO_F("version = %s (%zu chars)\n",
                   options_data->version, strlen(options_data->version));
    }
    
    config_destroy(config);
    free(config);
    
    return ret_code;
}

int parse_global_options(int argc, char *argv[], global_options_data_t *options_data, int start_index)
{
    int c;
    int finished = 0;
    int last_opt_index = optind = start_index;
    
    LOG_DEBUG("Parsing global options...\n");
    while (!finished && (c = getopt_long (argc, argv, "A:E:N:O:S:V:U:Z:", global_options, &optind)) != 1)
    {
        switch (c)
        {
            case 'A':
                options_data->vcf_filename = strdup(optarg);
                last_opt_index = optind;
                LOG_INFO_F("input vcf filename = %s\n", options_data->vcf_filename);
                break;
            case 'E':
                options_data->ped_filename = strdup(optarg);
                last_opt_index = optind;
                LOG_INFO_F("input ped filename = %s\n", options_data->ped_filename);
                break;
            case 'N':
                options_data->output_directory = strdup(optarg);
                last_opt_index = optind;
                LOG_INFO_F("output directory = %s\n", options_data->output_directory);
                break;
            case 'O':
                options_data->output_filename = strdup(optarg);
                last_opt_index = optind;
                LOG_INFO_F("output filename template = %s\n", options_data->output_filename);
                break;
            case 'S':
                options_data->species = strdup(optarg);
                last_opt_index = optind;
                LOG_INFO_F("species = %s\n", options_data->species);
                break;
            case 'V':
                options_data->version = strdup(optarg);
                last_opt_index = optind;
                LOG_INFO_F("version = %s\n", options_data->version);
                break;
            case 'U':
                options_data->host_url = strdup(optarg);
                last_opt_index = optind;
                LOG_INFO_F("web services URL = %s\n", options_data->host_url);
                break;
            case 'Z':
                last_opt_index = optind;
                LOG_DEBUG("--config found but not processed");
                break;
            case '?':
                LOG_WARN("undefined common option\n");
            default:
                finished = 1;
                LOG_INFO("global opts - default\n");
        }
    }
    
    return last_opt_index;
}

int verify_global_options(global_options_data_t *options_data)
{
    // Nothing to verify yet
    return 0;
}

struct option *merge_options(struct option local_options[], size_t num_local_options)
{
    size_t opts_size = num_local_options + NUM_GLOBAL_OPTIONS;
    struct option *tool_options = (struct option*) malloc (opts_size * sizeof(struct option));
    
    // Add global options
    for (int i = 0; i < NUM_GLOBAL_OPTIONS; i++)
    {
        tool_options[i] = global_options[i];
    }
    
    // Add local options
    for (int i = 0; i < num_local_options; i++)
    {
        tool_options[NUM_GLOBAL_OPTIONS + i] = local_options[i];
    }
    
    return tool_options;
}
