/*
 * Copyright (c) 2013 Alejandro Alemán Ramos (ICM-CIPF)
 * Copyright (c) 2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2013 Ignacio Medina (ICM-CIPF)
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
        argtable = merge_annot_options(annot_options, shared_options, arg_end(NUM_ANNOT_OPTIONS));
        show_usage("hpg-var-vcf annot", argtable);
        arg_freetable(argtable, NUM_ANNOT_OPTIONS);
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

    // Step 5: Create the web service request with all the parameters
    const int num_urls = 2;
    char **urls = malloc (num_urls * sizeof(char*));
    urls[0] = compose_cellbase_ws_request(shared_options_data->host_url, shared_options_data->version, shared_options_data->species, 
                                          "genomic/variant", "consequence_type");
    urls[1] = compose_cellbase_ws_request(shared_options_data->host_url, shared_options_data->version, shared_options_data->species, 
                                          "genomic/position", "snp");
    
    LOG_DEBUG_F("URL #1 = '%s'\n", urls[0]);
    LOG_DEBUG_F("URL #2 = '%s'\n", urls[1]);
    
    // Step 6: Perform the requested task
    int result = run_annot(urls, shared_options_data, options_data);

    // Step 7: Free memory
    for (int i = 0; i < num_urls; i++) {
        free(urls[i]);
    }
    free(urls);
    
    free_annot_options_data(options_data);
    free_shared_options_data(shared_options_data);
    arg_freetable(argtable, NUM_ANNOT_OPTIONS);

    return 0;

}

annot_options_t *new_annot_cli_options() {
    annot_options_t *options = (annot_options_t*) malloc (sizeof(annot_options_t));
    options->bam_directory = arg_str1(NULL, "bamdir", NULL, "The folder where the input BAM files are stored");
    options->missing = arg_lit0(NULL, "missing", "Annotate the missings values");
    options->dbsnp = arg_lit0(NULL, "dbsnp", "Annotate the dbSNP id");
    options->effect = arg_lit0(NULL, "effect", "Annotate the Effect");
    //options->phase = arg_lit0(NULL, "phase", "Annotate the Phase");   // TODO To implement
    options->all = arg_lit0(NULL, "all", "Activate all annotations");
    options->num_options = NUM_ANNOT_OPTIONS;
    return options;
}

annot_options_data_t *new_annot_options_data(annot_options_t *options) {
    annot_options_data_t *options_data = (annot_options_data_t*) malloc (sizeof(annot_options_data_t));
    options_data->bam_directory = strdup(*(options->bam_directory->sval));
    options_data->missing = options->missing->count;
    options_data->dbsnp = options->dbsnp->count;
    options_data->effect = options->effect->count;
    //options_data->phase = options->phase->count;   // TODO To implement
    
    if (options->all->count > 0) {
        options_data->missing = 1;
        //options_data->phase = 1;
        options_data->effect = 1;
        options_data->dbsnp = 1;
    }
    
    return options_data;
}

void free_annot_options_data(annot_options_data_t *options_data) {
    free(options_data->bam_directory);
    free(options_data);
}

