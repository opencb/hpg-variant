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

#include "epistasis.h"
#include "epistasis_runner.h"

int epistasis(int argc, char *argv[], const char *configuration_file) {

    /* ******************************
     *      Modifiable options      *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(0);
    epistasis_options_t *epistasis_options = new_epistasis_cli_options();

    // If no arguments or only -h / --help are provided, show usage
    void **argtable;
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_epistasis_options(epistasis_options, shared_options, arg_end(epistasis_options->num_options + shared_options->num_options));
        show_usage("hpg-var-gwas epi", argtable, epistasis_options->num_options + shared_options->num_options);
        arg_freetable(argtable, 12);
        return 0;
    }

    
    /* ******************************
     *       Execution steps        *
     * ******************************/

#ifdef _USE_MPI
    
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    
    // Only one rank reads the configuration file (libconfig is not compatible with MPI I/O)
    if (mpi_rank == 0) {
#endif
        
    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_epistasis_configuration(configuration_file, epistasis_options, shared_options);

    if (config_errors) {
        LOG_ERROR("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    argtable = parse_epistasis_options(argc, argv, epistasis_options, shared_options);

    // Step 3: check that all options are set with valid values
    // Mandatory options that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_epistasis_opts = verify_epistasis_options(epistasis_options, shared_options);
    if (check_epistasis_opts > 0) {
        return check_epistasis_opts;
    }
        
#ifdef _USE_MPI
    }
#endif
    
    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    // * rank 0 creates the structures and transmits info about them
    // * rank != 0 receives data to be associated with empty structures
    shared_options_data_t *shared_options_data;
    epistasis_options_data_t *epistasis_options_data;
    
#ifdef _USE_MPI
    if (mpi_rank == 0) {
        shared_options_data = new_shared_options_data(shared_options);
        epistasis_options_data = new_epistasis_options_data(epistasis_options);
    } else {
        shared_options_data = calloc (1, sizeof(shared_options_data_t));
        epistasis_options_data = calloc (1, sizeof(epistasis_options_data_t));
    }
    
    bcast_shared_options_data_mpi(shared_options_data, 0, MPI_COMM_WORLD);
    bcast_epistasis_options_data_mpi(epistasis_options_data, 0, MPI_COMM_WORLD);
#else
    shared_options_data = new_shared_options_data(shared_options);
    epistasis_options_data = new_epistasis_options_data(epistasis_options);
#endif
    
    // Step 5: Execute request and manage its response (as CURL request callback function)
    int result = run_epistasis(shared_options_data, epistasis_options_data);
    
    free_epistasis_options_data(epistasis_options_data);
    free_shared_options_data(shared_options_data);
#ifdef _USE_MPI
    if (mpi_rank == 0) {
#endif

    arg_freetable(argtable, 12);
    
#ifdef _USE_MPI
    }
#endif

    return 0;
}

epistasis_options_t *new_epistasis_cli_options(void) {
    epistasis_options_t *options = (epistasis_options_t*) malloc (sizeof(epistasis_options_t));
    options->num_options = NUM_EPISTASIS_OPTIONS;
    options->dataset_filename = arg_file1("d", "dataset", NULL, "Binary dataset used as input");
    options->max_ranking_size = arg_int0(NULL, "rank-size", NULL, "Number of best models saved");
    options->num_cv_repetitions = arg_int0(NULL, "num-cv-runs", NULL, "Number of times the k-fold cross-validation process is run");
    options->num_folds = arg_int0(NULL, "num-folds", NULL, "Number of folds in a k-fold cross-validation");
    options->stride = arg_int0(NULL, "stride", NULL, "Number of SNPs per block partition of the dataset");
    options->order = arg_int1(NULL, "order", NULL, "Number of SNPs to be combined at the same time");
    options->evaluation_subset = arg_str0(NULL, "eval-subset", NULL, "Whether to used training (default) or testing partitions when evaluating the best models");
    options->evaluation_mode = arg_str0(NULL, "eval-mode", NULL, "Whether to rank risky combinations by their CV-C or CV-A (values can be 'count' or 'accu')");
    return options;
}

epistasis_options_data_t *new_epistasis_options_data(epistasis_options_t *options) {
    epistasis_options_data_t *options_data = (epistasis_options_data_t*) malloc (sizeof(epistasis_options_data_t));
    options_data->dataset_filename = strdup(*(options->dataset_filename->filename));
    options_data->eval_subset = strcmp(*(options->evaluation_subset->sval), "testing") ? TRAINING : TESTING ;
    options_data->eval_mode = strcmp(*(options->evaluation_mode->sval), "count") ? CV_A : CV_C ;
    options_data->max_ranking_size = *(options->max_ranking_size->ival);
    options_data->num_cv_repetitions = *(options->num_cv_repetitions->ival);
    options_data->num_folds = *(options->num_folds->ival);
    options_data->order = *(options->order->ival);
    options_data->stride = *(options->stride->ival);
    return options_data;
}

void free_epistasis_options_data(epistasis_options_data_t *options_data) {
    free(options_data->dataset_filename);
    free(options_data);
}
