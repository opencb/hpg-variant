#include "filter.h"

int vcf_tool_filter(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    global_options_data_t *global_options_data = new_global_options_data();
    filter_options_data_t *options_data = new_filter_options_data();


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_global_configuration(configuration_file, global_options_data);
    config_errors &= read_filter_configuration(configuration_file, options_data);
    LOG_INFO_F("Config read with errors = %d\n", config_errors);

    if (config_errors) {
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    parse_filter_options(argc, argv, options_data, global_options_data);

    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_global_opts, check_vcf_tools_opts;

    check_global_opts = verify_global_options(global_options_data);
    if (check_global_opts > 0)
    {
        return check_global_opts;
    }

    check_vcf_tools_opts = verify_filter_options(global_options_data, options_data);
    if (check_vcf_tools_opts > 0)
    {
        return check_vcf_tools_opts;
    }

    // Step 4: Perform the requested task
    int result = run_filter(global_options_data, options_data);

    free_filter_options_data(options_data);
    free_global_options_data(global_options_data);

    return 0;
}


filter_options_data_t *new_filter_options_data()
{
    filter_options_data_t *options_data = (filter_options_data_t*) malloc (sizeof(filter_options_data_t));

    options_data->num_threads = 1;
    options_data->max_batches = 10;
    options_data->batch_size = 20000;

    return options_data;
}

void free_filter_options_data(filter_options_data_t *options_data)
{
    free(options_data);
}

