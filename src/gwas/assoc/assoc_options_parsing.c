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

#include "assoc.h"


int read_assoc_configuration(const char *filename, assoc_options_t *assoc_options, shared_options_t *shared_options) {
    if (filename == NULL || assoc_options == NULL || shared_options == NULL) {
        return -1;
    }

    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }

    // Read number of threads that will make request to the web service
    ret_code = config_lookup_int(config, "gwas.assoc.num-threads", shared_options->num_threads->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of threads not found in config file, must be set via command-line\n");
    } else {
        LOG_DEBUG_F("num-threads = %ld\n", *(shared_options->num_threads->ival));
    }

    // Read maximum number of batches that can be stored at certain moment
    ret_code = config_lookup_int(config, "gwas.assoc.max-batches", shared_options->max_batches->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Maximum number of batches not found in configuration file, must be set via command-line\n");
    } else {
        LOG_DEBUG_F("max-batches = %ld\n", *(shared_options->max_batches->ival));
    }
    
    // Read size of a batch (in lines or bytes)
    ret_code = config_lookup_int(config, "gwas.assoc.batch-lines", shared_options->batch_lines->ival);
    ret_code |= config_lookup_int(config, "gwas.assoc.batch-bytes", shared_options->batch_bytes->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Neither batch lines nor bytes found in configuration file, must be set via command-line\n");
    }
    
    config_destroy(config);
    free(config);

    return 0;
}

void **parse_assoc_options(int argc, char *argv[], assoc_options_t *assoc_options, shared_options_t *shared_options) {
    struct arg_end *end = arg_end(assoc_options->num_options + shared_options->num_options);
    void **argtable = merge_assoc_options(assoc_options, shared_options, end);
    
    int num_errors = arg_parse(argc, argv, argtable);
    if (num_errors > 0) {
        arg_print_errors(stdout, end, "hpg-var-gwas");
    }
    
    return argtable;
}

void **merge_assoc_options(assoc_options_t *assoc_options, shared_options_t *shared_options, struct arg_end *arg_end) {
    void **tool_options = malloc (30 * sizeof(void*));
    
    // Input/output files
    tool_options[0] = shared_options->vcf_filename;
    tool_options[1] = shared_options->ped_filename;
    tool_options[2] = shared_options->output_filename;
    tool_options[3] = shared_options->output_directory;
    
    // Species
    tool_options[4] = shared_options->species;
    
    // Association test arguments
    tool_options[5] = assoc_options->chisq;
    tool_options[6] = assoc_options->fisher;

    // Filter arguments
    tool_options[7] = shared_options->num_alleles;
    tool_options[8] = shared_options->coverage;
    tool_options[9] = shared_options->quality;
    tool_options[10] = shared_options->maf;
    tool_options[11] = shared_options->missing;
    tool_options[12] = shared_options->gene;
    tool_options[13] = shared_options->region;
    tool_options[14] = shared_options->region_file;
    tool_options[15] = shared_options->region_type;
    tool_options[16] = shared_options->snp;
    tool_options[17] = shared_options->indel;
    tool_options[18] = shared_options->dominant;
    tool_options[19] = shared_options->recessive;
    
    // Configuration file
    tool_options[20] = shared_options->log_level;
    tool_options[21] = shared_options->config_file;
    
    // Advanced configuration
    tool_options[22] = shared_options->host_url;
    tool_options[23] = shared_options->version;
    tool_options[24] = shared_options->max_batches;
    tool_options[25] = shared_options->batch_lines;
    tool_options[26] = shared_options->batch_bytes;
    tool_options[27] = shared_options->num_threads;
    tool_options[28] = shared_options->mmap_vcf_files;
    
    tool_options[29] = arg_end;
    
    return tool_options;
}


int verify_assoc_options(assoc_options_t *assoc_options, shared_options_t *shared_options) {
    // Check whether the input VCF file is defined
    if (shared_options->vcf_filename->count == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    // Check whether the task to perform is defined
    if (assoc_options->chisq->count + assoc_options->fisher->count == 0) {
        LOG_ERROR("Please specify the task to perform.\n");
        return GWAS_TASK_NOT_SPECIFIED;
    }

    // Check whether more than one task is specified
    if (assoc_options->chisq->count + assoc_options->fisher->count > 1) {
        LOG_ERROR("Please specify only one task to perform.\n");
        return GWAS_MANY_TASKS_SPECIFIED;
    }
    
    // Check whether the input PED file is defined
    if (shared_options->ped_filename->filename == NULL || strlen(*(shared_options->ped_filename->filename)) == 0) {
        LOG_ERROR("Please specify the input PED file.\n");
        return PED_FILE_NOT_SPECIFIED;
    }
    
    // Checker whether batch lines or bytes are defined
    if (*(shared_options->batch_lines->ival) == 0 && *(shared_options->batch_bytes->ival) == 0) {
        LOG_ERROR("Please specify the size of the reading batches (in lines or bytes).\n");
        return BATCH_SIZE_NOT_SPECIFIED;
    }
    
    // Checker if both batch lines or bytes are defined
    if (*(shared_options->batch_lines->ival) > 0 && *(shared_options->batch_bytes->ival) > 0) {
        LOG_WARN("The size of reading batches has been specified both in lines and bytes. The size in bytes will be used.\n");
        return 0;
    }
    
    return 0;
}
