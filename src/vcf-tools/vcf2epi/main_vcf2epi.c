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

#include "vcf2epi.h"
#include "dataset_creator.h"

int vcf_tool_vcf2epi(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     * 	    Modifiable options	    *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(1);
    vcf2epi_options_t *vcf2epi_options = new_vcf2epi_cli_options();

    // If no arguments or only -h / --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_vcf2epi_options(vcf2epi_options, shared_options, arg_end(vcf2epi_options->num_options + shared_options->num_options));
        show_usage("hpg-vcf-var vcf2epi", argtable, vcf2epi_options->num_options + shared_options->num_options);
        arg_freetable(argtable, vcf2epi_options->num_options + shared_options->num_options);
        return 0;
    }

    
    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_vcf2epi_configuration(configuration_file, vcf2epi_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    argtable = parse_vcf2epi_options(argc, argv, vcf2epi_options, shared_options);
    
    // Step 3: check that all options are set with valid values
    // Mandatory options that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_vcf2epi_opts = verify_vcf2epi_options(vcf2epi_options, shared_options);
    if (check_vcf2epi_opts > 0) {
        return check_vcf2epi_opts;
    }
    
    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    vcf2epi_options_data_t *vcf2epi_options_data = new_vcf2epi_options_data(vcf2epi_options);

    // Step 5: Execute request and manage its response (as CURL request callback function)
    int result = create_dataset_from_vcf(shared_options_data);
    
    free_vcf2epi_options_data(vcf2epi_options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, vcf2epi_options->num_options + shared_options->num_options - 1);

    return 0;
}

vcf2epi_options_t *new_vcf2epi_cli_options(void) {
    vcf2epi_options_t *options = (vcf2epi_options_t*) malloc (sizeof(vcf2epi_options_t));
    options->num_options = NUM_EPISTASIS_OPTIONS;
    return options;
}

vcf2epi_options_data_t *new_vcf2epi_options_data(vcf2epi_options_t *options) {
    vcf2epi_options_data_t *options_data = (vcf2epi_options_data_t*) malloc (sizeof(vcf2epi_options_data_t));
    return options_data;
}

void free_vcf2epi_options_data(vcf2epi_options_data_t *options_data) {
    free(options_data);
}
