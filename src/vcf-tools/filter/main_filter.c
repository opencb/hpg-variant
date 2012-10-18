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

#include "filter.h"

int vcf_tool_filter(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options();
    filter_options_t *filter_options = new_filter_cli_options();

    // If no arguments or only --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_filter_options(filter_options, shared_options, arg_end(filter_options->num_options + shared_options->num_options));
        show_usage("hpg-var-vcf", argtable, filter_options->num_options + shared_options->num_options);
        arg_freetable(argtable, filter_options->num_options + shared_options->num_options);
        return 0;
    }

    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_filter_configuration(configuration_file, filter_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }
    
    // Step 2: parse command-line options
    argtable = parse_filter_options(argc, argv, filter_options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_vcf_tools_opts = verify_filter_options(filter_options, shared_options);
    if (check_vcf_tools_opts > 0)
    {
        return check_vcf_tools_opts;
    }

    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    filter_options_data_t *options_data = new_filter_options_data(filter_options, shared_options);

    // Step 5: Perform the requested task
    int result = run_filter(shared_options_data, options_data);

    free_filter_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, filter_options->num_options + shared_options->num_options);

    return 0;
}

filter_options_t *new_filter_cli_options(void) {
    filter_options_t *options = (filter_options_t*) malloc (sizeof(filter_options_t));
    options->num_options = NUM_FILTER_OPTIONS;
    
    options->num_alleles = arg_int0(NULL, "alleles", NULL, "Filter: by number of alleles");
    options->coverage = arg_int0(NULL, "coverage", NULL, "Filter: by minimum coverage");
    options->quality = arg_int0(NULL, "quality", NULL, "Filter: by minimum quality");
    options->maf = arg_dbl0(NULL, "maf", NULL, "Filter: by maximum MAF (minimum allele frequency)");
    options->region = arg_str0(NULL, "region", NULL, "Filter: by a list of regions (chr1:start1-end1,chr2:start2-end2...)");
    options->region_file = arg_file0(NULL, "region-file", NULL, "Filter: by a list of regions (read from a GFF file)");
    options->snp = arg_str0(NULL, "snp", NULL, "Filter: by being a SNP or not");
    
    options->save_rejected = arg_lit0(NULL, "save-rejected", "Write a file containing the rejected records");
    return options;
}


filter_options_data_t *new_filter_options_data(filter_options_t *options, shared_options_t *shared_options) {
    filter_options_data_t *options_data = (filter_options_data_t*) calloc (1, sizeof(filter_options_data_t));
    
    filter_t *filter;
    if (options->num_alleles->count > 0) {
        filter = num_alleles_filter_new(*(options->num_alleles->ival));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("number of alleles filter = %d\n", *(options->num_alleles->ival));
    }
    if (options->coverage->count > 0) {
        filter = coverage_filter_new(*(options->coverage->ival));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("minimum coverage filter = %d\n", *(options->coverage->ival));
    }
    if (options->quality->count > 0) {
        filter = quality_filter_new(*(options->quality->ival));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("minimum quality filter = %d\n", *(options->quality->ival));
    }
    if (options->maf->count > 0) {
        filter = maf_filter_new(*(options->maf->dval));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("maximum MAF = %.3f\n", ((maf_filter_args*)filter->args)->max_maf);
    }
    if (options->snp->count > 0) {
        if (!strcmp(*(options->snp->sval), "exclude")) {
            filter = snp_filter_new(0);
        } else {
            filter = snp_filter_new(1);
        }
//         filter = snp_filter_new(strdup(*(options->snp->sval)));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("snp filter to %s SNPs\n", *(options->snp->sval));
    }
    if (options->region->count > 0) {
        filter = region_exact_filter_new(strdup(*(options->region->sval)), 0,
                                            *(shared_options->host_url->sval), *(shared_options->species->sval), *(shared_options->version->sval));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("regions = %s\n", *(options->region->sval));
    } 
    if (options->region_file->count > 0) {
        filter = region_exact_filter_new(strdup(*(options->region->sval)), 1, 
                                            *(shared_options->host_url->sval), *(shared_options->species->sval), *(shared_options->version->sval));
        options_data->chain = add_to_filter_chain(filter, options_data->chain);
        LOG_DEBUG_F("regions file = %s\n", *(options->region->sval));
    }
    
    options_data->save_rejected = (options->save_rejected->count > 0);

    return options_data;
}

void free_filter_options_data(filter_options_data_t *options_data) {
    if (options_data->chain) {
        free_filter_chain(options_data->chain);
    };
    free(options_data);
}

