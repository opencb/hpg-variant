/*
 * Copyright (c) 2012-2014 Cristina Yenyxe Gonzalez Garcia (EMBL-EBI)
 * Copyright (c) 2012-2014 Ignacio Medina (EMBL-EBI)
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

#include "aggregate.h"


int vcf_tool_aggregate(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(0);
    aggregate_options_t *aggregate_options = new_aggregate_cli_options();

    // If no arguments or only --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_aggregate_options(aggregate_options, shared_options, arg_end(NUM_AGGREGATE_OPTIONS));
        show_usage("hpg-var-vcf aggregate", argtable);
        arg_freetable(argtable, NUM_AGGREGATE_OPTIONS);
        return 0;
    }


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_aggregate_configuration(configuration_file, aggregate_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }
    
    // Step 2: parse command-line options
    argtable = parse_aggregate_options(argc, argv, aggregate_options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_vcf_tools_opts = verify_aggregate_options(aggregate_options, shared_options);
    if (check_vcf_tools_opts > 0) {
        return check_vcf_tools_opts;
    }

    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    aggregate_options_data_t *options_data = new_aggregate_options_data(aggregate_options);

    init_log_custom(shared_options_data->log_level, 1, "hpg-var-vcf.log", "w");

    // Step 5: Perform the requested task
    int result = run_aggregate(shared_options_data, options_data);

    free_aggregate_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, NUM_AGGREGATE_OPTIONS);

    return result;
}

aggregate_options_t *new_aggregate_cli_options() {
    aggregate_options_t *options = (aggregate_options_t*) malloc (sizeof(aggregate_options_t));
    options->overwrite = arg_lit0(NULL, "overwrite", "Overwrite fields in the INFO column (AC, AF, AN...)");
    return options;
}

aggregate_options_data_t *new_aggregate_options_data(aggregate_options_t *options) {
    aggregate_options_data_t *options_data = (aggregate_options_data_t*) malloc (sizeof(aggregate_options_data_t));
    options_data->overwrite = options->overwrite->count;
    return options_data;
}

void free_aggregate_options_data(aggregate_options_data_t *options_data) {
    free(options_data);
}

