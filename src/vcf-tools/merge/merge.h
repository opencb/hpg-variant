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

#ifndef VCF_TOOLS_MERGE_H
#define VCF_TOOLS_MERGE_H

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bioformats/vcf/vcf_annotation.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_stats.h>
#include <bioformats/vcf/vcf_util.h>
#include <commons/log.h>
#include <commons/string_utils.h>
#include <commons/config/libconfig.h>
#include <containers/khash.h>
#include <containers/list.h>
#include <containers/cprops/linked_list.h>

#include "error.h"
#include "hpg_variant_utils.h"
#include "shared_options.h"

#define NUM_MERGE_OPTIONS   19


#define MERGED_RECORD       1
#define MERGED_HEADER       2
#define MERGED_DELIMITER    3

KHASH_SET_INIT_STR(names);

enum missing_mode { MISSING, REFERENCE };

typedef struct merge_options {
    struct arg_str *input_files;        /**< List of files used as input */
    struct arg_str *missing_mode;       /**< How to fill a missing sample field whenever its data is missing */
    struct arg_str *info_fields;        /**< Attributes of the new INFO fields generated */
    struct arg_lit *strict_reference;   /**< Whether to reject variants whose reference allele is not the same in all files */
    struct arg_lit *copy_filter;        /**< Whether to copy the contents of the original FILTER field into the samples */
    struct arg_lit *copy_info;          /**< Whether to copy the contents of the original INFO field into the samples */
} merge_options_t;

typedef struct merge_options_data {
    char **input_files;     /**< List of files used as input */
    char **info_fields;     /**< List of attributes of the new INFO fields generated */
    
    int num_files;          /**< Number of files used as input */
    int num_info_fields;    /**< Number of attributes of the new INFO fields generated */ 
    
    int strict_reference;   /**< Whether to reject variants whose reference allele is not the same in all files */
    int copy_filter;        /**< Whether to copy the contents of the original FILTER field into the samples */
    int copy_info;          /**< Whether to copy the contents of the original INFO field into the samples */
    
    enum missing_mode missing_mode;   /**< How to fill a missing sample field whenever its data is missing */
    
    array_list_t *config_search_paths; /**< Paths to search for the configuration files specific to the merge tool */

} merge_options_data_t;

typedef struct {
    vcf_record_t *record;
    vcf_file_t *file;
} vcf_record_file_link;



static merge_options_t *new_merge_cli_options(void);

/**
 * Initialize a merge_options_data_t structure mandatory fields.
 */
static merge_options_data_t *new_merge_options_data(merge_options_t *options, array_list_t *config_search_paths);

/**
 * Free memory associated to a merge_options_data_t structure.
 */
static void free_merge_options_data(merge_options_data_t *options_data);


vcf_record_file_link *vcf_record_file_link_new(vcf_record_t *record, vcf_file_t *file);

void vcf_record_file_link_free(vcf_record_file_link *link);


/* ******************************
 *       Tool execution         *
 * ******************************/

int merge_vcf_headers(vcf_file_t **files, int num_files, merge_options_data_t *options, list_t *output_list);

array_list_t *merge_vcf_sample_names(vcf_file_t **files, int num_files);

int merge_vcf_records(array_list_t **records_by_position, int num_positions, vcf_file_t **files, int num_files, merge_options_data_t *options, list_t *output_list);

vcf_record_t *merge_position(vcf_record_file_link **position_in_files, int position_occurrences,
                             vcf_file_t **files, int num_files, merge_options_data_t *options, int *err_code);


char *merge_id_field(vcf_record_file_link **position_in_files, int position_occurrences);

float merge_quality_field(vcf_record_file_link **position_in_files, int position_occurrences);

char *merge_alternate_field(vcf_record_file_link **position_in_files, int position_occurrences, cp_hashtable *alleles_table);

char *merge_filter_field(vcf_record_file_link **position_in_files, int position_occurrences);

char *merge_info_field(vcf_record_file_link **position_in_files, int position_occurrences, char **info_fields, int num_fields,
                       vcf_record_t *output_record, cp_hashtable *alleles, char *empty_sample);

char *merge_format_field(vcf_record_file_link **position_in_files, int position_occurrences, merge_options_data_t *options, array_list_t *format_fields);

array_list_t *merge_samples(vcf_record_file_link **position_in_files, int position_occurrences, vcf_file_t **files, int num_files, 
                            cp_hashtable *alleles_table, array_list_t *format_fields, int *format_indices, char *empty_sample, 
                            int gt_pos, int filter_pos, int info_pos, merge_options_data_t *options);


/* ******************************
 *      Options parsing         *
 * ******************************/

/**
 * Read the basic configuration parameters of the tool. If the configuration
 * file can't be read, these parameters should be provided via the command-line
 * interface.
 * 
 * @param filename File the options data are read from
 * @param options_data Local options values (filtering, sorting...) 
 * 
 * @return If the configuration has been successfully read
 */
int read_merge_configuration(const char *filename, merge_options_t *options_data, shared_options_t *shared_options);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param global_options_data
 */
void **parse_merge_options(int argc, char *argv[], merge_options_t *options_data, shared_options_t *shared_options_data);

void **merge_merge_options(merge_options_t *merge_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * 
 * @param options_data
 */
int verify_merge_options(merge_options_t *options_data, shared_options_t *shared_options_data);



#endif
