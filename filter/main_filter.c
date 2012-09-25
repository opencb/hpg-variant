#include "filter.h"

int vcf_tool_filter(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options();
    filter_options_t *options = new_filter_cli_options();


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_filter_configuration(configuration_file, options, shared_options);
    LOG_INFO_F("Config read with errors = %d\n", config_errors);

    if (config_errors) {
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    void **argtable = parse_filter_options(argc, argv, options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_vcf_tools_opts = verify_filter_options(options, shared_options);
    if (check_vcf_tools_opts > 0)
    {
        return check_vcf_tools_opts;
    }

    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    filter_options_data_t *options_data = new_filter_options_data(options, shared_options);

    // Step 5: Perform the requested task
    int result = run_filter(shared_options_data, options_data);

    free_filter_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, options->num_options + shared_options->num_options);

    return 0;
}

filter_options_t *new_filter_cli_options(void) {
    filter_options_t *options = (filter_options_t*) malloc (sizeof(filter_options_t));
    options->num_options = NUM_FILTER_OPTIONS;
    
    options->num_alleles = arg_int0(NULL, "alleles", NULL, "Filter: by number of alleles");
    options->coverage = arg_int0(NULL, "coverage", NULL, "Filter: by minimum coverage");
    options->quality = arg_int0(NULL, "quality", NULL, "Filter: by minimum quality");
    options->region = arg_str0(NULL, "region", NULL, "Filter: by a list of regions (chr1:start1-end1,chr2:start2-end2...)");
    options->region_file = arg_file0(NULL, "region-file", NULL, "Filter: by a list of regions (read from a GFF file)");
    options->snp = arg_str0(NULL, "snp", NULL, "Filter: by being a SNP or not");
    
//     options->excludes = arg_str0(NULL, "exclude", NULL, "Consequence types to exclude from the query");
    return options;
}


filter_options_data_t *new_filter_options_data(filter_options_t *options, shared_options_t *shared_options) {
    filter_options_data_t *options_data = (filter_options_data_t*) calloc (1, sizeof(filter_options_data_t));
    
    filter_t *filter;
    if (options->num_alleles->count > 0) {
        filter = num_alleles_filter_new(*(options->num_alleles->ival));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_INFO_F("number of alleles filter = %d\n", *(options->num_alleles->ival));
    }
    if (options->coverage->count > 0) {
        filter = coverage_filter_new(*(options->coverage->ival));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_INFO_F("minimum coverage filter = %d\n", *(options->coverage->ival));
    }
    if (options->quality->count > 0) {
        filter = quality_filter_new(*(options->quality->ival));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_INFO_F("minimum quality filter = %d\n", *(options->quality->ival));
    }
    if (options->snp->count > 0) {
        if (!strcmp(*(options->snp->sval), "exclude")) {
            filter = snp_filter_new(0);
        } else {
            filter = snp_filter_new(1);
        }
//         filter = snp_filter_new(strdup(*(options->snp->sval)));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_INFO_F("snp filter to %s SNPs\n", *(options->snp->sval));
    }
    if (options->region->count > 0) {
        filter = region_exact_filter_new(strdup(*(options->region->sval)), 0,
                                            *(shared_options->host_url->sval), *(shared_options->species->sval), *(shared_options->version->sval));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_INFO_F("regions = %s\n", *(options->region->sval));
    } 
    if (options->region_file->count > 0) {
        filter = region_exact_filter_new(strdup(*(options->region->sval)), 1, 
                                            *(shared_options->host_url->sval), *(shared_options->species->sval), *(shared_options->version->sval));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_INFO_F("regions file = %s\n", *(options->region->sval));
    }
//     options_data->num_threads = 1;
//     options_data->max_batches = 10;
//     options_data->batch_size = 20000;

    return options_data;
}

void free_filter_options_data(filter_options_data_t *options_data) {
    if (options_data->chain) {
        free_filter_chain(options_data->chain);
    };
    free(options_data);
}

