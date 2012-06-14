#include "effect.h"
#include "effect_runner.h"

int effect(int argc, char *argv[], const char *configuration_file) {
    LOG_DEBUG_F("effect called with %d args\n", argc);

    /* ******************************
     * 	    Modifiable options	    *
     * ******************************/

    shared_options_t *shared_options = new_global_options();
    effect_options_t *effect_options = init_options();


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
    void *argtable = parse_effect_options(argc, argv, effect_options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory options that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_global_opts, check_effect_opts;

    check_global_opts = verify_global_options(shared_options);
    if (check_global_opts > 0)
    {
        return check_global_opts;
    }

    check_effect_opts = verify_effect_options(shared_options, effect_options);
    if (check_effect_opts > 0)
    {
        return check_effect_opts;
    }
    
    // TODO create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_global_options_data(shared_options);
    effect_options_data_t *effect_options_data = init_options_data(effect_options);

    // Step 4: Create the web service request with all the parameters
    char *url = compose_effect_ws_request(shared_options_data);

    // Step 5: Execute request and manage its response (as CURL request callback function)
    int result = run_effect(url, shared_options_data, effect_options_data);

    free(url);
    free_options_data(effect_options_data);
    free_global_options_data(shared_options_data);
    arg_freetable(argtable, effect_options->num_options + shared_options->num_options);

    return 0;
}

effect_options_t *init_options(void) {
    effect_options_t *options = (effect_options_t*) malloc (sizeof(effect_options_t));
    options->excludes = arg_str0(NULL, "exclude", NULL, "Consequence types to exclude from the query");
    return options;
}

effect_options_data_t *init_options_data(effect_options_t *options) {
    effect_options_data_t *options_data = (effect_options_data_t*) malloc (sizeof(effect_options_data_t));
    
    options_data->excludes = options->excludes->sval;
//     options_data->chain = NULL;
//     options_data->excludes = NULL;
//     options_data->max_batches = 10;
//     options_data->batch_size = 2000;
//     options_data->num_threads = 4;
//     options_data->variants_per_request = 1000;

    return options_data;
}

void free_options_data(effect_options_data_t *options_data) {
//     if (options_data->chain)    { cp_heap_destroy(options_data->chain); }
    if (options_data->excludes) { free((void*) options_data->excludes); }
    free(options_data);
}
