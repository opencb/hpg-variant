/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
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

#include "epistasis.h"
#include "dataset_creator.h"

int main(int argc, char *argv[]) {

    /* ******************************
     * 	    Modifiable options	    *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options();
    epistasis_options_t *epistasis_options = new_epistasis_cli_options();
    void **argtable;

    // If no arguments or only -h / --help are provided, show usage
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_epistasis_options(epistasis_options, shared_options, arg_end(epistasis_options->num_options + shared_options->num_options));
        show_usage(argv[0], argtable, epistasis_options->num_options + shared_options->num_options);
        arg_freetable(argtable, epistasis_options->num_options + shared_options->num_options);
        return 0;
    }

    /* ******************************
     * 	    Execution steps	        *
     * ******************************/

    init_log_custom(2, 1, "hpg-var-epistasis.log", "w");
    
//    const char *configuration_file = find_configuration_file(argc, argv);
    array_list_t *config_search_paths = get_configuration_search_paths(argc, argv);
    const char *configuration_file = retrieve_config_file("hpg-variant.conf", config_search_paths);
    
    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_epistasis_configuration(configuration_file, epistasis_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    argtable = parse_epistasis_options(argc, argv, epistasis_options, shared_options);
    
    // Step 3: check that all options are set with valid values
    // Mandatory options that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_epistasis_opts = verify_epistasis_options(epistasis_options, shared_options);
    if (check_epistasis_opts > 0) {
        return check_epistasis_opts;
    }
    
    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    epistasis_options_data_t *epistasis_options_data = new_epistasis_options_data(epistasis_options);

    // Step 5: Execute request and manage its response (as CURL request callback function)
    int result = create_dataset_from_vcf(shared_options_data);
    
    free_epistasis_options_data(epistasis_options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, epistasis_options->num_options + shared_options->num_options - 1);

    stop_log();
    
    return 0;
}

epistasis_options_t *new_epistasis_cli_options(void) {
    epistasis_options_t *options = (epistasis_options_t*) malloc (sizeof(epistasis_options_t));
    options->num_options = NUM_EPISTASIS_OPTIONS;
    options->no_phenotypes = arg_lit0(NULL, "no-phenotypes", "Flag asking not to retrieve phenotypical information");
    options->excludes = arg_str0(NULL, "exclude", NULL, "Consequence types to exclude from the query (comma-separated)");
    return options;
}

epistasis_options_data_t *new_epistasis_options_data(epistasis_options_t *options) {
    epistasis_options_data_t *options_data = (epistasis_options_data_t*) malloc (sizeof(epistasis_options_data_t));
    options_data->no_phenotypes = options->no_phenotypes->count;
    options_data->excludes = strdup(*(options->excludes->sval));
    return options_data;
}

void free_epistasis_options_data(epistasis_options_data_t *options_data) {
    if (options_data->excludes) { free(options_data->excludes); }
    free(options_data);
}
