#include "split.h"
#include "split_runner.h"


int vcf_tool_split(int argc, char *argv[]) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    global_options_data_t *global_options_data = new_global_options_data();
    split_options_data_t *options_data = new_split_options_data();


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_global_configuration("hpg-vcf-tools.cfg", global_options_data);
    config_errors &= read_split_configuration("hpg-vcf-tools.cfg", options_data);
    LOG_INFO_F("Config read with errors = %d\n", config_errors);

    if (config_errors) {
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    parse_split_options(argc, argv, options_data, global_options_data);
    
    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_global_opts, check_vcf_tools_opts;

    check_global_opts = verify_global_options(global_options_data);
    if (check_global_opts > 0)
    {
        return check_global_opts;
    }

    check_vcf_tools_opts = verify_split_options(global_options_data, options_data);
    if (check_vcf_tools_opts > 0)
    {
        return check_vcf_tools_opts;
    }

    // Step 4: Perform the requested task
    int result = run_split(global_options_data, options_data);

    free_split_options_data(options_data);
    free_global_options_data(global_options_data);

    return 0;
}


split_options_data_t *new_split_options_data() {
    split_options_data_t *options_data = (split_options_data_t*) malloc (sizeof(split_options_data_t));

    options_data->criterion = NONE;
    options_data->num_threads = 1;
    options_data->max_batches = 10;
    options_data->batch_size = 8000;
    options_data->variants_per_thread = 2000;

    return options_data;
}

void free_split_options_data(split_options_data_t *options_data) {
    free(options_data);
}

