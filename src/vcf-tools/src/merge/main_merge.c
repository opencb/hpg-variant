#include "merge.h"
#include "merge_runner.h"


int vcf_tool_merge(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options();
    merge_options_t *options = new_merge_cli_options();


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_merge_configuration(configuration_file, options, shared_options);
    LOG_INFO_F("Config read with errors = %d\n", config_errors);

    if (config_errors) {
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    void **argtable = parse_merge_options(argc, argv, options, shared_options);
    
    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_vcf_tools_opts = verify_merge_options(options, shared_options);
    if (check_vcf_tools_opts > 0) {
        return check_vcf_tools_opts;
    }

    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    merge_options_data_t *options_data = new_merge_options_data(options);

    // Step 5: Perform the requested task
    int result = run_merge(shared_options_data, options_data);

    free_merge_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, options->num_options + shared_options->num_options + 1 - 2);

    return 0;
}

merge_options_t *new_merge_cli_options() {
    merge_options_t *options = (merge_options_t*) malloc (sizeof(merge_options_t));
    options->num_options = NUM_MERGE_OPTIONS;
    options->input_files = arg_str1(NULL, "vcf-list", NULL, "List of comma-separated input VCF files");
    options->missing_mode = arg_str0(NULL, "missing-mode", NULL, "How to fill missing genotypes (missing = ./., reference = 0/0)");
    options->info_fields = arg_str1(NULL, "info-fields", NULL, "Information to generate in the new INFO column");
    options->copy_filter = arg_lit0(NULL, "copy-filter", "Whether to copy the FILTER column from the original files into the samples");
    options->copy_info = arg_lit0(NULL, "copy-info", "Whether to copy the INFO column from the original files into the samples");
    return options;
}

merge_options_data_t *new_merge_options_data(merge_options_t *options) {
    merge_options_data_t *options_data = (merge_options_data_t*) calloc (1, sizeof(merge_options_data_t));
    if (!strcasecmp(*(options->missing_mode->sval), "missing")) {
        options_data->missing_mode = MISSING;
    } else if (!strcasecmp(*(options->missing_mode->sval), "reference")) {
        options_data->missing_mode = REFERENCE;
    }
    options_data->input_files = split(*(options->input_files->sval), ",", &(options_data->num_files));
    options_data->info_fields = split(*(options->info_fields->sval), ",", &(options_data->num_info_fields));
    
    options_data->copy_filter = options->copy_filter->count;
    options_data->copy_info = options->copy_info->count;
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

