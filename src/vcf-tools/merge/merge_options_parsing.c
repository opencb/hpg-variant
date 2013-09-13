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


int read_merge_configuration(const char *filename, merge_options_t *options, shared_options_t *shared_options) {
    if (filename == NULL || options == NULL || shared_options == NULL) {
        return -1;
    }
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("config file error: %s\n", config_error_text(config));
        return ret_code;
    }
    
    char *tmp_string;
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "vcf-tools.merge.num-threads", shared_options->num_threads->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of threads not found in config file, must be set via command-line\n");
    } else {
        LOG_DEBUG_F("num-threads = %ld\n", *(shared_options->num_threads->ival));
    }

    // Read maximum number of batches that can be stored at certain moment
    ret_code = config_lookup_int(config, "vcf-tools.merge.max-batches", shared_options->max_batches->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Maximum number of batches not found in configuration file, must be set via command-line\n");
    } else {
        LOG_DEBUG_F("max-batches = %ld\n", *(shared_options->max_batches->ival));
    }
    
    // Read size of every batch read
    ret_code = config_lookup_int(config, "vcf-tools.merge.batch-lines", shared_options->batch_lines->ival);
    ret_code |= config_lookup_int(config, "vcf-tools.merge.batch-bytes", shared_options->batch_bytes->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Neither batch lines nor bytes found in configuration file, must be set via command-line\n");
    } 
    /*else {
        LOG_DEBUG_F("batch-lines = %ld\n", *(shared_options->batch_size->ival));
    }*/
    
    // Read missing mode
    ret_code = config_lookup_string(config, "vcf-tools.merge.missing-mode", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Missing mode not found in configuration file, must be set via command-line\n");
    } else {
        *(options->missing_mode->sval) = strdup(tmp_string);
        LOG_DEBUG_F("missing mode = %s (%zu chars)\n",
                   *(options->missing_mode->sval), strlen(*(options->missing_mode->sval)));
    }

    config_destroy(config);
    free(config);

    return 0;
}

void **parse_merge_options(int argc, char *argv[], merge_options_t *merge_options, shared_options_t *shared_options) {
    struct arg_end *end = arg_end(NUM_MERGE_OPTIONS);
    void **argtable = merge_merge_options(merge_options, shared_options, end);
    
    assert(argv);
    assert(argtable);
    int num_errors = arg_parse(argc, argv, argtable);
    if (num_errors > 0) {
        arg_print_errors(stdout, end, "hpg-var-vcf");
    }
    
    return argtable;
}

void **merge_merge_options(merge_options_t *merge_options, shared_options_t *shared_options, struct arg_end *arg_end) {
    void **tool_options = malloc (NUM_MERGE_OPTIONS * sizeof(void*));
    // Input/output files
    tool_options[0] = merge_options->input_files;
    tool_options[1] = shared_options->output_filename;
    tool_options[2] = shared_options->output_directory;
    
    // Species
    tool_options[3] = shared_options->species;
    
    // Merge options
    tool_options[4] = merge_options->missing_mode;
    tool_options[5] = merge_options->strict_reference;
    tool_options[6] = merge_options->copy_filter;
    tool_options[7] = merge_options->copy_info;
    tool_options[8] = merge_options->info_fields;
    
    // Configuration file
    tool_options[9] = shared_options->log_level;
    tool_options[10] = shared_options->config_file;
    
    // Advanced configuration
    tool_options[11] = shared_options->host_url;
    tool_options[12] = shared_options->version;
    tool_options[13] = shared_options->max_batches;
    tool_options[14] = shared_options->batch_lines;
    tool_options[15] = shared_options->batch_bytes;
    tool_options[16] = shared_options->num_threads;
    tool_options[17] = shared_options->mmap_vcf_files;
    
    tool_options[18] = arg_end;
    
    return tool_options;
}


int verify_merge_options(merge_options_t *merge_options, shared_options_t *shared_options) {
    // Check whether the input VCF files are defined
    if (merge_options->input_files->count == 0) {
        LOG_ERROR("Please specify the input VCF files.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    // Check whether the missing mode is defined
    if (merge_options->missing_mode->count == 0 && strlen(*(merge_options->missing_mode->sval)) == 0) {
        LOG_ERROR("Please specify how to fill missing samples information (mark as missing/reference).\n");
        return MISSING_MODE_NOT_SPECIFIED;
    }
    
    return 0;
}
