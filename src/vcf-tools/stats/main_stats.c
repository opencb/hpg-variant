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

#include "stats.h"


int vcf_tool_stats(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(0);
    stats_options_t *stats_options = new_stats_cli_options();

    // If no arguments or only --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_stats_options(stats_options, shared_options, arg_end(NUM_STATS_OPTIONS));
        show_usage("hpg-var-vcf stats", argtable);
        arg_freetable(argtable, NUM_STATS_OPTIONS);
        return 0;
    }


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_stats_configuration(configuration_file, stats_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }
    
    // Step 2: parse command-line options
    argtable = parse_stats_options(argc, argv, stats_options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_vcf_tools_opts = verify_stats_options(stats_options, shared_options);
    if (check_vcf_tools_opts > 0) {
        return check_vcf_tools_opts;
    }

    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    stats_options_data_t *options_data = new_stats_options_data(stats_options);

    init_log_custom(shared_options_data->log_level, 1, "hpg-var-vcf.log", "w");

    // Step 5: Perform the requested task
    int result = run_stats(shared_options_data, options_data);

    free_stats_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, NUM_STATS_OPTIONS);

    return result;
}

stats_options_t *new_stats_cli_options() {
    stats_options_t *options = (stats_options_t*) malloc (sizeof(stats_options_t));
    options->sample_stats = arg_lit0(NULL, "samples", "Get statistics about samples");
    options->variant_stats = arg_lit0(NULL, "variants", "Get statistics about variants, both per variant and per file (default)");
    options->save_db = arg_lit0(NULL, "db", "Save statistics to SQLite3 database file");
    options->variable = arg_str0(NULL, "variable", NULL, "Name for the variable field");
    options->variable_groups = arg_str0(NULL, "variable-group", NULL, "Sequence of variable groups");
    options->phenotype = arg_str0(NULL, "phenotype",NULL, "Affected,Unaffected phenotype values");
    
    return options;
}

stats_options_data_t *new_stats_options_data(stats_options_t *options) {
    stats_options_data_t *options_data = (stats_options_data_t*) malloc (sizeof(stats_options_data_t));
    options_data->sample_stats = options->sample_stats->count;
    options_data->variant_stats = options->variant_stats->count;
    options_data->save_db = options->save_db->count;

    options_data->variable = options->variable->count? strdup(*(options->variable->sval)) :NULL;
    options_data->variable_groups = options->variable_groups->count? strdup(*(options->variable_groups->sval)) :NULL;
    options_data->phenotype = options->phenotype->count? strdup(*options->phenotype->sval) :NULL;
    return options_data;
}

void free_stats_options_data(stats_options_data_t *options_data) {
    if(options_data->variable) free(options_data->variable);
    if(options_data->variable_groups) free(options_data->variable_groups);
    if(options_data->phenotype) free(options_data->phenotype);
    free(options_data);
}

