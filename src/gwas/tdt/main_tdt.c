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

#include "tdt.h"
#include "tdt_runner.h"

int tdt(int argc, char *argv[], const char *configuration_file) {
    
    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(1);
    tdt_options_t *tdt_options = new_tdt_cli_options();

    // If no arguments or only --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "--help")) {
        argtable = merge_tdt_options(tdt_options, shared_options, arg_end(NUM_TDT_OPTIONS));
        show_usage("hpg-var-gwas tdt", argtable);
        arg_freetable(argtable, NUM_TDT_OPTIONS);
        return 0;
    }

    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_tdt_configuration(configuration_file, tdt_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }
    
    // Step 2: parse command-line options
    // If no arguments or only --help are provided, show usage
    argtable = parse_tdt_options(argc, argv, tdt_options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory options that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_tdt_opts = verify_tdt_options(tdt_options, shared_options);
    if (check_tdt_opts > 0) {
        return check_tdt_opts;
    }
    
    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
//     tdt_options_data_t *options_data = new_tdt_options_data(tdt_options);

    init_log_custom(shared_options_data->log_level, 1, "hpg-var-gwas.log", "w");
    
    // Step 5: Perform the operations related to the selected GWAS sub-tool
//     run_tdt_test(shared_options_data, options_data);
    run_tdt_test(shared_options_data);
    
//     free_tdt_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, NUM_TDT_OPTIONS);
    free(configuration_file);

    return 0;
}

tdt_options_t *new_tdt_cli_options(void) {
    tdt_options_t *options = (tdt_options_t*) malloc (sizeof(tdt_options_t));
    return options;
}
