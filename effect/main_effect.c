#include "effect.h"
#include "effect_runner.h"

int effect(int argc, char *argv[], const char *configuration_file) {
    LOG_DEBUG_F("effect called with %d args\n", argc);

    /* ******************************
     * 	    Modifiable options	    *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options();
    effect_options_t *effect_options = new_effect_cli_options();


    /* ******************************
     * 	    Execution steps	        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_global_configuration(configuration_file, shared_options);
    config_errors &= read_effect_configuration(configuration_file, effect_options, shared_options);
    LOG_INFO_F("Config read with errors = %d\n", config_errors);
    
    if (config_errors) {
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    void **argtable = parse_effect_options(argc, argv, effect_options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory options that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_effect_opts = verify_effect_options(effect_options, shared_options);
    if (check_effect_opts > 0) {
        return check_effect_opts;
    }
    
    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    effect_options_data_t *effect_options_data = new_effect_options_data(effect_options);

    // Step 5: Create the web service request with all the parameters
    char *url = compose_effect_ws_request(shared_options_data);

    // Step 6: Execute request and manage its response (as CURL request callback function)
    int result = run_effect(url, shared_options_data, effect_options_data);

    free(url);
    free_effect_options_data(effect_options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, effect_options->num_options + shared_options->num_options);

    return 0;
}

effect_options_t *new_effect_cli_options(void) {
    effect_options_t *options = (effect_options_t*) malloc (sizeof(effect_options_t));
    options->num_options = NUM_EFFECT_OPTIONS;
    options->excludes = arg_str0(NULL, "exclude", NULL, "Consequence types to exclude from the query");
    return options;
}

effect_options_data_t *new_effect_options_data(effect_options_t *options) {
    effect_options_data_t *options_data = (effect_options_data_t*) malloc (sizeof(effect_options_data_t));
    options_data->excludes = strdup(*(options->excludes->sval));
    return options_data;
}

void free_effect_options_data(effect_options_data_t *options_data) {
    if (options_data->excludes) { free(options_data->excludes); }
    free(options_data);
}
