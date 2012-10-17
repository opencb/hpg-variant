/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hpg-variant is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-variant. If not, see <http://www.gnu.org/licenses/>.
 */

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
    const int num_urls = 3;
    char **urls = malloc (num_urls * sizeof(char*));
    urls[0] = compose_effect_ws_request("genomic/variant", "consequence_type", shared_options_data);
    urls[1] = compose_effect_ws_request("feature/snp", "phenotype", shared_options_data);
    urls[2] = compose_effect_ws_request("genomic/variant", "mutation_phenotype", shared_options_data);

    LOG_DEBUG_F("URL #1 = '%s'\nURL #2 = '%s'\nURL #3 = '%s'\n", urls[0], urls[1], urls[2]);
    
    // Step 6: Execute request and manage its response (as CURL request callback function)
    int result = run_effect(urls, shared_options_data, effect_options_data);

    // Step 7: Free memory
    for (int i = 0; i < num_urls; i++) {
        free(urls[i]);
    }
    free(urls);
    
    free_effect_options_data(effect_options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, effect_options->num_options + shared_options->num_options);

    return 0;
}

effect_options_t *new_effect_cli_options(void) {
    effect_options_t *options = (effect_options_t*) malloc (sizeof(effect_options_t));
    options->num_options = NUM_EFFECT_OPTIONS;
    options->no_phenotypes = arg_lit0(NULL, "no-phenotypes", "Flag asking not to retrieve phenotypical information");
    options->excludes = arg_str0(NULL, "exclude", NULL, "Consequence types to exclude from the query");
    return options;
}

effect_options_data_t *new_effect_options_data(effect_options_t *options) {
    effect_options_data_t *options_data = (effect_options_data_t*) malloc (sizeof(effect_options_data_t));
    options_data->no_phenotypes = options->no_phenotypes->count;
    options_data->excludes = strdup(*(options->excludes->sval));
    return options_data;
}

void free_effect_options_data(effect_options_data_t *options_data) {
    if (options_data->excludes) { free(options_data->excludes); }
    free(options_data);
}
