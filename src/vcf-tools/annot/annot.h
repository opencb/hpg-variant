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

#ifndef VCF_TOOLS_ANNOT_H
#define VCF_TOOLS_ANNOT_H

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>

#include <bioformats/db/db_utils.h>
#include <bioformats/db/cellbase_connector.h>
#include <bioformats/features/variant/variant_effect.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/bam/bam_file.h>
#include <bioformats/vcf/vcf_db.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_stats.h>
#include <bioformats/vcf/vcf_stats_report.h>
#include <bioformats/vcf/vcf_util.h>
#include <commons/file_utils.h>
#include <commons/log.h>
#include <config/libconfig.h>
#include <sqlite/sqlite3.h>
#include <containers/khash.h>

#include "error.h"
#include "shared_options.h"
#include "hpg_variant_utils.h"

#define NUM_ANNOT_OPTIONS       16
#define MAX_VARIANTS_PER_QUERY  1000
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


typedef struct annot_options {
    struct arg_str *bam_directory; /**< BAM DIRECTORY */
    struct arg_lit *missing;
    struct arg_lit *dbsnp;
    struct arg_lit *effect;
//    struct arg_lit *phase;    // TODO To implement
    struct arg_lit *all;
    
    int num_options;
} annot_options_t;

/**
 * @struct annot_options_data
 * 
 */
typedef struct annot_options_data {
    char *bam_directory; /**< BAM DIRECTORY */
    int missing;
    int dbsnp;
    int effect;
//    int phase;    // TODO To implement
} annot_options_data_t;


typedef struct vcf_annot_sample {
    char *name;
    array_list_t *chromosomes;
} vcf_annot_sample_t;

typedef struct vcf_annot_chr {
    char *name;
    array_list_t *positions;
} vcf_annot_chr_t;

typedef struct vcf_annot_pos {
    unsigned int pos;
    int dp;
} vcf_annot_pos_t;

typedef struct vcf_annot_bam {
    char* bam_filename;
    char* bai_filename;
} vcf_annot_bam_t;


// When counting records instead of printing them,
// data passed to the bam_fetch callback is encapsulated in this struct.
typedef struct {
	bam_header_t *header;
	int *count;
} count_func_data_t;

static annot_options_t *new_annot_cli_options(void);

/**
 * Initialize a annot_options_data_t structure mandatory fields.
 */
static annot_options_data_t *new_annot_options_data(annot_options_t *options);

/**
 * Free memory associated to a annot_options_data_t structure.
 */
static void free_annot_options_data(annot_options_data_t *options_data);


/* ******************************
 *       Tool execution         *
 * ******************************/

int run_annot(char** urls, shared_options_data_t *shared_options_data, annot_options_data_t *options_data);


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
int read_annot_configuration(const char *filename, annot_options_t *options, shared_options_t *shared_options);

/**
 * 
 * @param argc
 * @param argv
 * @param options_data
 * @param shared_options_data
 */
void **parse_annot_options(int argc, char *argv[], annot_options_t *annot_options, shared_options_t *shared_options);

void **merge_annot_options(annot_options_t *annot_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * 
 * @param options
 */
int verify_annot_options(annot_options_t *annot_options, shared_options_t *shared_options);


KHASH_MAP_INIT_STR(bams, char*);



int vcf_annot_chr_cmp(const vcf_annot_chr_t *v1, const vcf_annot_chr_t *v2);

void vcf_annot_sample_free(vcf_annot_sample_t *sample);

void vcf_annot_chr_free(vcf_annot_chr_t *chr);

void vcf_annot_pos_free(vcf_annot_pos_t *pos);

int vcf_annot_process_chunk(vcf_record_t **variants, int num_variants, array_list_t *sample_list, vcf_file_t *vcf_file);

static void vcf_annot_sort_sample(vcf_annot_sample_t* annot_sample);

int vcf_annot_check_bams(vcf_annot_sample_t* annot_sample, khash_t(bams)* sample_bams);

int vcf_annot_edit_chunk(vcf_record_t **variants, int num_variants, array_list_t *sample_list, khash_t(bams)* sample_bams, vcf_file_t *vcf_file);

static vcf_annot_pos_t * vcf_annot_get_pos(char *sample_name, char *chr, size_t pos, array_list_t *sample_list);

static int vcf_annot_pos_cmp (const void * a, const void * b);

static int binary_search_pos(array_list_t* array_list , unsigned int f,unsigned int l, unsigned int target);
void set_field_value_in_sample(char **sample, int position, char* value);

vcf_annot_chr_t * vcf_annot_get_chr(char* chr, array_list_t* array_chr);

#endif
