#include "global_options.h"

global_options_data_t *init_global_options_data(void)
{
    global_options_data_t *options_data = (global_options_data_t*) malloc (sizeof(global_options_data_t));
    
    options_data->output_directory = (char*) calloc (5, sizeof(char));
    strncat(options_data->output_directory, "/tmp", 4);
    
    return options_data;
}

void free_global_options_data(global_options_data_t *options_data)
{
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
    
    config_t *config = (config_t*) malloc (sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE)
    {
        fprintf(stderr, "config file error: %s\n", config_error_text(config));
        return ret_code;
    }
    
    // Read output directory
    const char *tmp_string;
    ret_code = config_lookup_string(config, "global.outdir", &tmp_string);
    if (ret_code == CONFIG_FALSE)
    {
        fprintf(stderr, "output directory config error: %s\n", config_error_text(config));
        return ret_code;
    } else {
        free(options_data->output_directory);
        options_data->output_directory = (char*) calloc (strlen(tmp_string)+1, sizeof(char));
        strncat(options_data->output_directory, tmp_string, strlen(tmp_string));
        LOG_INFO_F("output directory = %s (%zu chars)\n", options_data->output_directory, 
                    strlen(options_data->output_directory));
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
    while (!finished && (c = getopt_long (argc, argv, ":A:O:", global_options, &optind)) != 1)
    {
        switch (c)
        {
            case 'A':
                options_data->vcf_filename = (char*) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(options_data->vcf_filename, optarg);
                last_opt_index = optind;
                LOG_INFO_F("input vcf filename = %s\n", options_data->vcf_filename);
                break;
            case 'N':
                options_data->output_directory = (char*) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(options_data->output_directory, optarg);
                last_opt_index = optind;
                LOG_INFO_F("output directory = %s\n", options_data->output_directory);
                break;
            case 'O':
                options_data->output_filename = (char*) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(options_data->output_filename, optarg);
                last_opt_index = optind;
                LOG_INFO_F("output filename template = %s\n", options_data->output_filename);
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
