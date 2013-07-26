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

#include "merge.h"
#include "merge_runner.h"


int vcf_tool_merge(int argc, char *argv[], const char *configuration_file, array_list_t *config_search_paths) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(0);
    merge_options_t *merge_options = new_merge_cli_options();

    // If no arguments or only --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_merge_options(merge_options, shared_options, arg_end(merge_options->num_options + shared_options->num_options));
        show_usage("hpg-var-vcf merge", argtable, merge_options->num_options + shared_options->num_options);
        arg_freetable(argtable, 19);
        return 0;
    }


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_merge_configuration(configuration_file, merge_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }
    
    // Step 2: parse command-line options
    // If no arguments or only --help are provided, show usage
    argtable = parse_merge_options(argc, argv, merge_options, shared_options);
    
    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_vcf_tools_opts = verify_merge_options(merge_options, shared_options);
    if (check_vcf_tools_opts > 0) {
        return check_vcf_tools_opts;
    }

    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    merge_options_data_t *options_data = new_merge_options_data(merge_options, config_search_paths);

    init_log_custom(shared_options_data->log_level, 1, "hpg-var-vcf.log", "w");

    // Step 5: Perform the requested task
    int result = run_merge(shared_options_data, options_data);

    free_merge_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, 19);

    return 0;
}

merge_options_t *new_merge_cli_options() {
    merge_options_t *options = (merge_options_t*) malloc (sizeof(merge_options_t));
    options->num_options = NUM_MERGE_OPTIONS;
    options->input_files = arg_str1(NULL, "vcf-list", NULL, "List of comma-separated input VCF files");
    options->missing_mode = arg_str0(NULL, "missing-mode", NULL, "How to fill missing genotypes (missing = ./., reference = 0/0)");
    options->info_fields = arg_str0(NULL, "info-fields", NULL, "Information to generate in the new INFO column");
    options->strict_reference = arg_lit0(NULL, "strict-ref", "Whether to reject variants whose reference allele is not the same in all files");
    options->copy_filter = arg_lit0(NULL, "copy-filter", "Whether to copy the FILTER column from the original files into the samples");
    options->copy_info = arg_lit0(NULL, "copy-info", "Whether to copy the INFO column from the original files into the samples");
    return options;
}

merge_options_data_t *new_merge_options_data(merge_options_t *options, array_list_t *config_search_paths) {
    merge_options_data_t *options_data = (merge_options_data_t*) calloc (1, sizeof(merge_options_data_t));
    if (!strcasecmp(*(options->missing_mode->sval), "missing")) {
        options_data->missing_mode = MISSING;
    } else if (!strcasecmp(*(options->missing_mode->sval), "reference")) {
        options_data->missing_mode = REFERENCE;
    }
    options_data->input_files = split(*(options->input_files->sval), ",", &(options_data->num_files));

    if (options->info_fields->count > 0) {
    	options_data->info_fields = split(*(options->info_fields->sval), ",", &(options_data->num_info_fields));
    } else {
    	options_data->num_info_fields = 0;
    }
    options_data->strict_reference = options->strict_reference->count;
    options_data->copy_filter = options->copy_filter->count;
    options_data->copy_info = options->copy_info->count;
    options_data->config_search_paths = config_search_paths;
    return options_data;
}

void free_merge_options_data(merge_options_data_t *options_data) {
    for (int i = 0; i < options_data->num_files; i++) {
        free(options_data->input_files[i]);
    }
    free(options_data->input_files);
    
    for (int i = 0; i < options_data->num_info_fields; i++) {
        free(options_data->info_fields[i]);
    }
    free(options_data->info_fields);
    
    free(options_data);
}

