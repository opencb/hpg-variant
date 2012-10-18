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

#include "gwas.h"
#include "assoc_runner.h"
#include "tdt_runner.h"

int genomic_analysis(int argc, char *argv[], const char *configuration_file) {
    LOG_DEBUG_F("gwas called with %d args\n", argc);

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options();
    gwas_options_t *gwas_options = new_gwas_cli_options();


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_gwas_configuration(configuration_file, gwas_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }
    
    // Step 2: parse command-line options
    // If no arguments or only --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "--help")) {
        argtable = merge_gwas_options(gwas_options, shared_options, arg_end(gwas_options->num_options + shared_options->num_options));
        show_usage("gwas", argtable, gwas_options->num_options + shared_options->num_options);
        arg_freetable(argtable, gwas_options->num_options + shared_options->num_options);
        return 0;
    } else {
        argtable = parse_gwas_options(argc, argv, gwas_options, shared_options);
    }

    // Step 3: check that all options are set with valid values
    // Mandatory options that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_gwas_opts = verify_gwas_options(gwas_options, shared_options);
    if (check_gwas_opts > 0) {
        return check_gwas_opts;
    }
    
    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    gwas_options_data_t *options_data = new_gwas_options_data(gwas_options);

    // Step 5: Perform the operations related to the selected GWAS sub-tool
    switch (options_data->task) {
        case TDT:
            run_tdt_test(shared_options_data, options_data);
        break;
        case ASSOCIATION_BASIC:
        case FISHER:
            run_association_test(shared_options_data, options_data);
        break;
    }
    
    free_gwas_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, gwas_options->num_options + shared_options->num_options);

    return 0;
}

gwas_options_t *new_gwas_cli_options(void) {
    gwas_options_t *options = (gwas_options_t*) malloc (sizeof(gwas_options_t));
    options->num_options = NUM_GWAS_OPTIONS;
    options->assoc = arg_lit0(NULL, "assoc", "Basic case/control association test");
    options->fisher = arg_lit0(NULL, "fisher", "Fisher's exact test");
    options->tdt = arg_lit0(NULL, "tdt", "Transmission disequilibrium test");
    return options;
}

gwas_options_data_t *new_gwas_options_data(gwas_options_t *options) {
    gwas_options_data_t *options_data = (gwas_options_data_t*) calloc (1, sizeof(gwas_options_data_t));
    if (options->assoc->count > 0) {
        options_data->task = ASSOCIATION_BASIC;
    } else if (options->fisher->count > 0) {
        options_data->task = FISHER;
    } else if (options->tdt->count > 0) {
        options_data->task = TDT;
    }
    return options_data;
}

void free_gwas_options_data(gwas_options_data_t *options_data) {
    free(options_data);
}
