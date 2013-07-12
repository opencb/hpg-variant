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

#include "annot.h"


int vcf_tool_annot(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(0);
    annot_options_t *annot_options = new_annot_cli_options();

    // If no arguments or only --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_annot_options(annot_options, shared_options, arg_end(annot_options->num_options + shared_options->num_options));
        show_usage("hpg-var-vcf annot", argtable, annot_options->num_options + shared_options->num_options);
        arg_freetable(argtable, 15);
        return 0;
    }


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_annot_configuration(configuration_file, annot_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }
    
    // Step 2: parse command-line options
    argtable = parse_annot_options(argc, argv, annot_options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_vcf_tools_opts = verify_annot_options(annot_options, shared_options);
    if (check_vcf_tools_opts > 0) {
        return check_vcf_tools_opts;
    }

    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    annot_options_data_t *options_data = new_annot_options_data(annot_options);

    // Step 5: Perform the requested task
    int result = run_annot(shared_options_data, options_data);

    free_annot_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, 15);

    return 0;

}

annot_options_t *new_annot_cli_options() {
    annot_options_t *options = (annot_options_t*) malloc (sizeof(annot_options_t));
    options->sample_annot = arg_lit0(NULL, "samples", "Get statistics about samples");
    options->variant_annot = arg_lit0(NULL, "variants", "Get statistics about variants, both per variant and per file (default)");
    options->save_db = arg_lit0(NULL, "db", "Save statistics to SQLite3 database file");
    options->num_options = NUM_ANNOT_OPTIONS;
    return options;
}

annot_options_data_t *new_annot_options_data(annot_options_t *options) {
    annot_options_data_t *options_data = (annot_options_data_t*) malloc (sizeof(annot_options_data_t));
    options_data->sample_annot = options->sample_annot->count;
    options_data->variant_annot = options->variant_annot->count;
    options_data->save_db = options->save_db->count;
    return options_data;
}

void free_annot_options_data(annot_options_data_t *options_data) {
    free(options_data);
}

