/*
 * Copyright (c) 2012-2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
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

int main(int argc, char *argv[]) {

    /* ******************************
     * 	    Modifiable options	    *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(0);
    effect_options_t *effect_options = new_effect_cli_options();
    void **argtable;

    // If no arguments or only -h / --help are provided, show usage
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_effect_options(effect_options, shared_options, arg_end(NUM_EFFECT_OPTIONS));
        show_usage(argv[0], argtable);
        arg_freetable(argtable, NUM_EFFECT_OPTIONS);
        return 0;
    } else if (!strcmp(argv[1], "--version")) {
        show_version("Effect");
        return 0;
    }

    /* ******************************
     * 	    Execution steps	    *
     * ******************************/

    init_log_custom(LOG_DEFAULT_LEVEL, 1, "hpg-var-effect.log", "w");
    
    array_list_t *config_search_paths = get_configuration_search_paths(argc, argv);
    const char *configuration_file = retrieve_config_file("hpg-variant.conf", config_search_paths);
    
    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_effect_configuration(configuration_file, effect_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    argtable = parse_effect_options(argc, argv, effect_options, shared_options);
    
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

    init_log_custom(shared_options_data->log_level, 1, "hpg-var-effect.log", "w");
    
    // Step 5: Create the web service request with all the parameters
    const int num_urls = 3;
    char **urls = malloc (num_urls * sizeof(char*));
    urls[0] = compose_cellbase_ws_request(shared_options_data->host_url, shared_options_data->version, shared_options_data->species, 
                                          "genomic/variant", "consequence_type");
    urls[1] = compose_cellbase_ws_request(shared_options_data->host_url, shared_options_data->version, shared_options_data->species, 
                                          "feature/snp", "phenotype");
    urls[2] = compose_cellbase_ws_request(shared_options_data->host_url, shared_options_data->version, shared_options_data->species, 
                                          "genomic/variant", "mutation_phenotype");

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
    arg_freetable(argtable, NUM_EFFECT_OPTIONS);
    array_list_free(config_search_paths, free);
    free(configuration_file);

    stop_log();
    
    return 0;
}

effect_options_t *new_effect_cli_options(void) {
    effect_options_t *options = (effect_options_t*) malloc (sizeof(effect_options_t));
    options->no_phenotypes = arg_lit0(NULL, "no-phenotypes", "Flag asking not to retrieve phenotypical information");
    options->excludes = arg_str0(NULL, "exclude", NULL, "Consequence types to exclude from the query (comma-separated)");
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
