#include "stats.h"


int vcf_tool_stats(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options();
    stats_options_t *options = new_stats_cli_options();


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_stats_configuration(configuration_file, options, shared_options);
    LOG_INFO_F("Config read with errors = %d\n", config_errors);

    if (config_errors) {
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    void **argtable = parse_stats_options(argc, argv, options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_vcf_tools_opts = verify_stats_options(options, shared_options);
    if (check_vcf_tools_opts > 0) {
        return check_vcf_tools_opts;
    }

    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    stats_options_data_t *options_data = new_stats_options_data(options);

    // Step 5: Perform the requested task
    int result = run_stats(shared_options_data, options_data);

    free_stats_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, options->num_options + shared_options->num_options - 7);

    return 0;
}

stats_options_t *new_stats_cli_options() {
    stats_options_t *options = (stats_options_t*) malloc (sizeof(stats_options_t));
    options->num_options = NUM_STATS_OPTIONS;
    return options;
}

stats_options_data_t *new_stats_options_data(stats_options_t *options) {
    stats_options_data_t *options_data = (stats_options_data_t*) malloc (sizeof(stats_options_data_t));
    return options_data;
}

void free_stats_options_data(stats_options_data_t *options_data) {
    free(options_data);
}

