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

#include "filter.h"


int read_filter_configuration(const char *filename, filter_options_t *options, shared_options_t *shared_options) {
    if (filename == NULL || options == NULL || shared_options == NULL) {
        return -1;
    }
    
    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("config file error: %s\n", config_error_text(config));
        return ret_code;
    }
    
    const char *tmp_string;
    
    // Read number of threads that will make request to the web service
    ret_code = config_lookup_int(config, "vcf-tools.filter.num-threads", shared_options->num_threads->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of threads not found in config file, must be set via command-line\n");
    } else {
        LOG_DEBUG_F("num-threads = %ld\n", *(shared_options->num_threads->ival));
    }

    // Read maximum number of batches that can be stored at certain moment
    ret_code = config_lookup_int(config, "vcf-tools.filter.max-batches", shared_options->max_batches->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Maximum number of batches not found in configuration file, must be set via command-line\n");
    } else {
        LOG_DEBUG_F("max-batches = %ld\n", *(shared_options->max_batches->ival));
    }
    
    // Read size of every batch read
    ret_code = config_lookup_int(config, "vcf-tools.filter.batch-lines", shared_options->batch_lines->ival);
    ret_code |= config_lookup_int(config, "vcf-tools.filter.batch-bytes", shared_options->batch_bytes->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Neither batch lines nor bytes found in configuration file, must be set via command-line\n");
    } 
    /*else {
        LOG_DEBUG_F("batch-lines = %ld\n", *(shared_options->batch_size->ival));
    }*/
    
    config_destroy(config);
    free(config);

    return 0;
}


void **parse_filter_options(int argc, char *argv[], filter_options_t *filter_options, shared_options_t *shared_options) {
    struct arg_end *end = arg_end(NUM_FILTER_OPTIONS);
    void **argtable = merge_filter_options(filter_options, shared_options, end);
    
    int num_errors = arg_parse(argc, argv, argtable);
    if (num_errors > 0) {
        arg_print_errors(stdout, end, "hpg-var-vcf");
    }
    
    return argtable;
}

void **merge_filter_options(filter_options_t *filter_options, shared_options_t *shared_options, struct arg_end *arg_end) {
    void **tool_options = malloc (NUM_FILTER_OPTIONS * sizeof(void*));
    // Input/output files
    tool_options[0] = shared_options->vcf_filename;
    tool_options[1] = shared_options->ped_filename;
    tool_options[2] = shared_options->output_directory;
    
    // Species
    tool_options[3] = shared_options->species;
    
    // Filter options
    tool_options[4] = shared_options->num_alleles;
    tool_options[5] = shared_options->coverage;
    tool_options[6] = shared_options->quality;
    tool_options[7] = shared_options->maf;
    tool_options[8] = shared_options->missing;
    tool_options[9] = shared_options->gene;
    tool_options[10] = shared_options->region;
    tool_options[11] = shared_options->region_file;
    tool_options[12] = shared_options->region_type;
    tool_options[13] = shared_options->snp;
    tool_options[14] = shared_options->variant_type;
    tool_options[15] = shared_options->indel;
    tool_options[16] = shared_options->dominant;
    tool_options[17] = shared_options->recessive;
    tool_options[18] = filter_options->save_rejected;
    
    // Configuration file
    tool_options[19] = shared_options->log_level;
    tool_options[20] = shared_options->config_file;
    
    // Advanced configuration
    tool_options[21] = shared_options->host_url;
    tool_options[22] = shared_options->version;
    tool_options[23] = shared_options->max_batches;
    tool_options[24] = shared_options->batch_lines;
    tool_options[25] = shared_options->batch_bytes;
    tool_options[26] = shared_options->num_threads;
    tool_options[27] = shared_options->mmap_vcf_files;
    
    tool_options[28] = arg_end;
    
    return tool_options;
}


int verify_filter_options(filter_options_t *filter_options, shared_options_t *shared_options) {
    // Check whether the input VCF file is defined
    if (shared_options->vcf_filename->count == 0) {
        LOG_ERROR("Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
    // Check whether a filter or more has been specified
    if (!shared_options->coverage->count && !shared_options->num_alleles->count && !shared_options->quality->count && 
        !shared_options->snp->count && !shared_options->region->count && !shared_options->region_file->count && 
        !shared_options->maf->count && !shared_options->gene->count && !shared_options->indel->count &&
        !shared_options->variant_type->count && !shared_options->dominant->count && !shared_options->recessive->count) {
        LOG_ERROR("Please specify at least one filter\n");
        return EMPTY_LIST_OF_FILTERS;
    }

    // Check whether the input PED file is defined when needed (for certain filters)
    if ((shared_options->dominant->count || shared_options->recessive->count) &&
        shared_options->ped_filename->count == 0) {
        LOG_ERROR("Please specify the input PED file.\n");
        return PED_FILE_NOT_SPECIFIED;
    }
    
    // Check whether the host URL is defined
    if (shared_options->host_url->sval == NULL || strlen(*(shared_options->host_url->sval)) == 0) {
        LOG_ERROR("Please specify the host URL to the web service.\n");
        return HOST_URL_NOT_SPECIFIED;
    }

    // Check whether the version is defined
    if (shared_options->version->sval == NULL || strlen(*(shared_options->version->sval)) == 0) {
        LOG_ERROR("Please specify the version.\n");
        return VERSION_NOT_SPECIFIED;
    }

    // Check whether the species is defined
    if (shared_options->species->sval == NULL || strlen(*(shared_options->species->sval)) == 0) {
        LOG_ERROR("Please specify the species to take as reference.\n");
        return SPECIES_NOT_SPECIFIED;
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
