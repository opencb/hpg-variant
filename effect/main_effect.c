#include "effect.h"
#include "effect_runner.h"

int effect(int argc, char *argv[])
{
	LOG_DEBUG_F("effect called with %d args\n", argc);
	
	/* ******************************
	 * 	    Modifiable options	    *
	 * ******************************/
	
	global_options_data_t *global_options_data = init_global_options_data();
	effect_options_data_t *options_data = init_options_data();
	
	
	/* ******************************
	 * 	    Execution steps	        *
	 * ******************************/
	
	// Step 1: read options from configuration file
	int config_read = read_global_configuration("hpg-variant.cfg", global_options_data);
	config_read &= read_effect_configuration("hpg-variant.cfg", options_data);
	LOG_INFO_F("Config read successfully = %d\n", config_read);
	
	// Step 2: parse command-line options
	parse_effect_options(argc, argv, options_data, global_options_data);
	
	// Step 3: check that all options are set with valid values
	// Mandatory options that couldn't be read from the config file must be set via command-line
	// If not, return error code!
	int check_global_opts, check_effect_opts;
	
	check_global_opts = verify_global_options(global_options_data);
	if (check_global_opts > 0)
	{
		return check_global_opts;
	}
	
	check_effect_opts = verify_effect_options(global_options_data, options_data);
	if (check_effect_opts > 0)
	{
		return check_effect_opts;
	}
	
	// Step 4: Create the web service request with all the parameters
	char *url = compose_effect_ws_request(options_data);
	
    // Step 5: Execute request and manage its response (as CURL request callback function)
    int result = execute_effect_query(url, global_options_data, options_data);

    free(url);
	free_options_data(options_data);
	free_global_options_data(global_options_data);
	
	return 0;
}

effect_options_data_t *init_options_data(void)
{
	effect_options_data_t *options_data = (effect_options_data_t*) malloc (sizeof(effect_options_data_t));
	
	options_data->host_url = (char*) malloc (256 * sizeof(char));
	options_data->version = (char*) malloc (32 * sizeof(char));
	options_data->species = (char*) malloc (32 * sizeof(char));
    options_data->max_batches = 10;
    options_data->batch_size = 2000;
    options_data->num_threads = 4;
    options_data->variants_per_request = 1000;
	
	return options_data;
}

void free_options_data(effect_options_data_t *options_data)
{
	if (options_data->host_url) { free((void*) options_data->host_url); }
	if (options_data->version)  { free((void*) options_data->version); }
	if (options_data->species)  { free((void*) options_data->species); }
	if (options_data->chain)    { cp_heap_destroy(options_data->chain); }
	if (options_data->excludes) { free((void*) options_data->excludes); }
	free(options_data);
}
