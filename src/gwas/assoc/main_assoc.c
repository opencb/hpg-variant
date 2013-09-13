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

#include "assoc.h"
#include "assoc_runner.h"

int association(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(1);
    assoc_options_t *assoc_options = new_assoc_cli_options();

    // If no arguments or only --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "--help")) {
        argtable = merge_assoc_options(assoc_options, shared_options, arg_end(NUM_ASSOC_OPTIONS));
        show_usage("hpg-var-gwas assoc", argtable);
        arg_freetable(argtable, NUM_ASSOC_OPTIONS);
        return 0;
    }

    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_assoc_configuration(configuration_file, assoc_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }
    
    // Step 2: parse command-line options
    argtable = parse_assoc_options(argc, argv, assoc_options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory options that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_assoc_opts = verify_assoc_options(assoc_options, shared_options);
    if (check_assoc_opts > 0) {
        return check_assoc_opts;
    }
    
    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    assoc_options_data_t *options_data = new_assoc_options_data(assoc_options);

    init_log_custom(shared_options_data->log_level, 1, "hpg-var-gwas.log", "w");
    
    // Step 5: Perform the GWAS test
    run_association_test(shared_options_data, options_data);
    
    free_assoc_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, NUM_ASSOC_OPTIONS);
    free(configuration_file);

    return 0;
}

assoc_options_t *new_assoc_cli_options(void) {
    assoc_options_t *options = (assoc_options_t*) malloc (sizeof(assoc_options_t));
    options->chisq = arg_lit0(NULL, "chisq", "Chi-square association test");
    options->fisher = arg_lit0(NULL, "fisher", "Fisher's exact test");
    return options;
}

assoc_options_data_t *new_assoc_options_data(assoc_options_t *options) {
    assoc_options_data_t *options_data = (assoc_options_data_t*) calloc (1, sizeof(assoc_options_data_t));
    if (options->chisq->count > 0) {
        options_data->task = CHI_SQUARE;
    } else if (options->fisher->count > 0) {
        options_data->task = FISHER;
    } else {
        options_data->task = NONE;
    }
    return options_data;
}

void free_assoc_options_data(assoc_options_data_t *options_data) {
    free(options_data);
}
