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

#include "split.h"
#include "split_runner.h"
#include "commons/string_utils.h"


int vcf_tool_split(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(0);
    split_options_t *split_options = new_split_cli_options();

    // If no arguments or only --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_split_options(split_options, shared_options, arg_end(split_options->num_options + shared_options->num_options));
        show_usage("hpg-var-vcf split", argtable, split_options->num_options + shared_options->num_options);
        arg_freetable(argtable, 12);
        return 0;
    }


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_split_configuration(configuration_file, split_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }
    
    // Step 2: parse command-line options
    argtable = parse_split_options(argc, argv, split_options, shared_options);
    
    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_vcf_tools_opts = verify_split_options(split_options, shared_options);
    if (check_vcf_tools_opts > 0) {
        return check_vcf_tools_opts;
    }

    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    split_options_data_t *options_data = new_split_options_data(split_options);

    init_log_custom(shared_options_data->log_level, 1, "hpg-var-vcf.log", "w");

    // Step 5: Perform the requested task
    int result = run_split(shared_options_data, options_data);

    free_split_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, 12);

    return 0;
}

split_options_t *new_split_cli_options() {
    split_options_t *options = (split_options_t*) malloc (sizeof(split_options_t));
    options->num_options = NUM_SPLIT_OPTIONS;
    options->criterion = arg_str1(NULL, "criterion", NULL, "Criterion for splitting the file");
    options->intervals = arg_str0(NULL, "intervals", NULL, "Values of intervals for splitting the file");
    return options;
}

split_options_data_t *new_split_options_data(split_options_t *options) {
    split_options_data_t *options_data = (split_options_data_t*) malloc (sizeof(split_options_data_t));

    if (!strcasecmp("chromosome", *(options->criterion->sval))) {
        options_data->criterion = SPLIT_CHROMOSOME;
        
    } else if (!strcasecmp("coverage", *(options->criterion->sval))) {
        options_data->criterion = SPLIT_COVERAGE;
	char* intervals_str = strdup(*(options->intervals->sval));
	char** tokens = split(intervals_str, ",", &(options_data->num_intervals));
	
	options_data->intervals = malloc(options_data->num_intervals * sizeof(long));
	for (int i = 0; i < options_data->num_intervals; i++) {
            options_data->intervals[i] = atol(tokens[i]);
            free(tokens[i]);
	}
	
        free(tokens);
        free(intervals_str);
    }

    return options_data;
}

void free_split_options_data(split_options_data_t *options_data) {
    free(options_data);
}

