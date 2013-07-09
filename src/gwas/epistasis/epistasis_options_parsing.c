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


int read_epistasis_configuration(const char *filename, epistasis_options_t *epistasis_options, shared_options_t *shared_options) {
    if (filename == NULL || epistasis_options == NULL || shared_options == NULL) {
        return -1;
    }
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }

    const char *tmp_string;
    
    // Read number of threads that will run epistasis checks in parallel
    ret_code = config_lookup_int(config, "gwas.epistasis.num-threads", shared_options->num_threads->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of threads not found in config file, must be set via command-line");
    } else {
        LOG_DEBUG_F("num-threads = %ld\n", *(shared_options->num_threads->ival));
    }

    // Read maximum number of SNPs per block
    ret_code = config_lookup_int(config, "gwas.epistasis.stride", epistasis_options->stride->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of SNPs per block partition not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("stride = %ld\n", *(epistasis_options->stride->ival));
    }
    
    // Read number of folds per k-fold cross-validation
    ret_code = config_lookup_int(config, "gwas.epistasis.num-folds", epistasis_options->num_folds->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of folds per k-fold cross-validation not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("num-folds = %ld\n", *(epistasis_options->num_folds->ival));
    }
    
    // Read number of times cross-validation will be run
    ret_code = config_lookup_int(config, "gwas.epistasis.num-cv-repetitions", epistasis_options->num_cv_repetitions->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of cross-validation repetitions not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("num-cv-repetitions = %ld\n", *(epistasis_options->num_cv_repetitions->ival));
    }
    
    // Read maximum number of best models recorded
    ret_code = config_lookup_int(config, "gwas.epistasis.max-ranking-size", epistasis_options->max_ranking_size->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Maximum number of best models recorded not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("max-ranking-size = %ld\n", *(epistasis_options->max_ranking_size->ival));
    }
    
    // Read whether the training of testing partitions will be used for evaluating best models
    ret_code = config_lookup_string(config, "gwas.epistasis.evaluation-subset", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Evaluation subset not found in configuration file, must be set via command-line");
    } else {
        *(epistasis_options->evaluation_subset->sval) = strdup(tmp_string);
        LOG_DEBUG_F("evaluation subset = %s (%zu chars)\n",
                   *(epistasis_options->evaluation_subset->sval), strlen(*(epistasis_options->evaluation_subset->sval)));
    }

    // Read whether the CV-c or CV-a will be used to rank risky combinations
    ret_code = config_lookup_string(config, "gwas.epistasis.evaluation-mode", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Evaluation mode not found in configuration file, must be set via command-line");
    } else {
        *(epistasis_options->evaluation_mode->sval) = strdup(tmp_string);
        LOG_DEBUG_F("evaluation mode = %s (%zu chars)\n",
                   *(epistasis_options->evaluation_mode->sval), strlen(*(epistasis_options->evaluation_mode->sval)));
    }

    config_destroy(config);
    free(config);

    return 0;
}

void **parse_epistasis_options(int argc, char *argv[], epistasis_options_t *epistasis_options, shared_options_t *shared_options) {
    struct arg_end *end = arg_end(epistasis_options->num_options + shared_options->num_options);
    void **argtable = merge_epistasis_options(epistasis_options, shared_options, end);
    
    int num_errors = arg_parse(argc, argv, argtable);
    if (num_errors > 0) {
        arg_print_errors(stdout, end, "hpg-var-gwas");
    }
    
    return argtable;
}

void **merge_epistasis_options(epistasis_options_t *epistasis_options, shared_options_t *shared_options, struct arg_end *arg_end) {
    void **tool_options = malloc (12 * sizeof(void*));
    // Input/output files
    tool_options[0] = epistasis_options->dataset_filename;
    tool_options[1] = shared_options->output_directory;
    
    // Epistasis arguments
    tool_options[2] = epistasis_options->order;
    tool_options[3] = epistasis_options->num_folds;
    tool_options[4] = epistasis_options->num_cv_repetitions;
    tool_options[5] = epistasis_options->max_ranking_size;
    tool_options[6] = epistasis_options->evaluation_subset;
    tool_options[7] = epistasis_options->evaluation_mode;
    tool_options[8] = epistasis_options->stride;
    
    // Configuration file
    tool_options[9] = shared_options->config_file;
    
    // Advanced configuration
    tool_options[10] = shared_options->num_threads;
    
    tool_options[11] = arg_end;
    
    return tool_options;
}


int verify_epistasis_options(epistasis_options_t *epistasis_options, shared_options_t *shared_options) {
    // Check whether the input dataset file is defined
    if (epistasis_options->dataset_filename->count == 0) {
        LOG_ERROR("Please specify the dataset file.\n");
        return EPISTASIS_DATASET_NOT_SPECIFIED;
    }
    
    // Check whether the order of combinations is defined
    if (epistasis_options->order->count == 0 || *(epistasis_options->order->ival) == 0) {
        LOG_ERROR("Please specify the number of SNPs to be combined at the same time.\n");
        return EPISTASIS_ORDER_NOT_SPECIFIED;
    }
    
    // Checker whether the number of folds in k-fold CV is defined
    if (*(epistasis_options->num_folds->ival) == 0) {
        LOG_ERROR("Please specify the number of folds in a k-fold cross-validation.\n");
        return EPISTASIS_FOLDS_NOT_SPECIFIED;
    }
    
    // Checker whether the number of CV runs is defined
    if (*(epistasis_options->num_cv_repetitions->ival) == 0) {
        LOG_ERROR("Please specify the times the cross-validation will be run.\n");
        return EPISTASIS_CV_RUNS_NOT_SPECIFIED;
    }
    
    // Checker whether to use the training or testing partition for evaluating the best models
    if (*(epistasis_options->evaluation_subset->sval) == NULL || 
        (strcmp(*(epistasis_options->evaluation_subset->sval), "training") && strcmp(*(epistasis_options->evaluation_subset->sval), "testing")) ) {
        LOG_ERROR("Please specify the dataset partition for evaluating the best models (training/testing).\n");
        return EPISTASIS_EVAL_SUBSET_NOT_SPECIFIED;
    }
    
    // Checker whether the number of SNPs per block partition of the dataset
    if (*(epistasis_options->stride->ival) == 0) {
        LOG_ERROR("Please specify the number of SNPs per block partition of the dataset.\n");
        return EPISTASIS_STRIDE_NOT_SPECIFIED;
    }
    
    return 0;
}
