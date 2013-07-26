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

#ifndef GLOBAL_OPTIONS_H
#define GLOBAL_OPTIONS_H

/** 
 * @file shared_options.h 
 * @brief Structures and functions associated to application-wide options
 * 
 * This file defines the structures which store the options for the whole application, and also 
 * functions to read their value from a configuration file or the command line.
 */ 

#include <stdlib.h>
#include <string.h>

#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_util.h>
#include <commons/log.h>
#include <commons/argtable/argtable2.h>
#include <commons/config/libconfig.h>

#include "error.h"

/**
 * Number of options applicable to the whole application.
 */
#define NUM_GLOBAL_OPTIONS  27

typedef struct shared_options {
    struct arg_file *vcf_filename;      /**< VCF file used as input. */
    struct arg_file *ped_filename;      /**< PED file used as input. */
    struct arg_file *output_filename;   /**< Filename template for the main output file. */
    struct arg_str *output_directory;   /**< Directory where the output files will be stored. */
    
    struct arg_str *host_url;           /**< URL of the host where the web service runs. */
    struct arg_str *version;            /**< Version of the WS to query. */
    struct arg_str *species;            /**< Species whose genome is taken as reference. */
    
    struct arg_int *max_batches;        /**< Maximum number of batches stored at the same time. */
    struct arg_int *batch_lines;        /**< Maximum size of a batch (in lines). */
    struct arg_int *batch_bytes;        /**< Maximum size of a batch (in bytes). */
    struct arg_int *num_threads;        /**< Number of threads when a task runs in parallel. */
    
    struct arg_int *coverage;           /**< Filter by coverage. */
    struct arg_dbl *maf;                /**< Filter by minimum allele frequency (MAF). */
    struct arg_dbl *missing;            /**< Filter by missing values. */
    struct arg_int *num_alleles;        /**< Filter by number of alleles. */
    struct arg_int *quality;            /**< Filter by quality. */
    struct arg_str *gene;               /**< Filter by gene */
    struct arg_str *region;             /**< Filter by region */
    struct arg_file *region_file;       /**< Filter by region (using a GFF file) */
    struct arg_str *region_type;        /**< Filter by region and type (used along with the 'region_file' arg) */
    struct arg_str *snp;                /**< Filter by SNP */
    struct arg_str *indel;              /**< Filter by indel */
    struct arg_dbl *dominant;           /**< Filter by samples following a dominant inheritance model */
    struct arg_dbl *recessive;          /**< Filter by samples following a recessive inheritance model */
    
    struct arg_str *log_level;          /**< Level to register in the log file */
    struct arg_file *config_file;       /**< Path to the configuration file */
    struct arg_lit *mmap_vcf_files;     /**< Whether to map VCF files to virtual memory or use the I/O API. */
    
    int num_options;
} shared_options_t;


/**
 * @brief Values for the application-wide options.
 * 
 * This struct contains the values for all the options of the application that can be applied to 
 * all its tools, such as filenames for different formats (VCF, GFF, BAM...) and the output files 
 * and folders.
 */
typedef struct shared_options_data {
    char *vcf_filename;                 /**< VCF file used as input. */
    char *ped_filename;                 /**< PED file used as input. */
    char *output_filename;              /**< Filename template for the main output file. */
    char *output_directory;             /**< Directory where the output files will be stored. */
    
    char *host_url;                     /**< URL of the host where the web service runs. */
    char *version;                      /**< Version of the WS to query. */
    char *species;                      /**< Species whose genome is taken as reference. */
    
    int max_batches;                    /**< Maximum number of batches stored at the same time. */
    int batch_lines;                    /**< Maximum size of a batch (in lines). */
    int batch_bytes;                    /**< Maximum size of a batch (in bytes). */
    int num_threads;                    /**< Number of threads when a task runs in parallel. */
    int entries_per_thread;             /**< Number of entries in a batch each thread processes. */
    
    filter_chain *chain;                /**< Chain of filters to apply to the VCF records, if that is the case. */
    
    int log_level;                      /**< Level to register in the log file */
} shared_options_data_t;


/**
 * @brief Initializes an global_options_t structure mandatory members.
 * @return A new global_options_t structure.
 * 
 * Initializes the only mandatory member of a global_options_t, which is the output directory.
 */
shared_options_t *new_shared_cli_options(int ped_required);

/**
 * @brief Initializes an global_options_data_t structure mandatory members.
 * @return A new global_options_data_t structure.
 * 
 * Initializes the only mandatory member of a global_options_data_t, which is the output directory.
 */
shared_options_data_t *new_shared_options_data(shared_options_t *options);

/**
 * @brief Free memory associated to a global_options_data_t structure.
 * @param options_data the structure to be freed
 * 
 * Free memory associated to a global_options_data_t structure, including its text buffers.
 */
void free_shared_options_data(shared_options_data_t *options_data);


/* **********************************************
 *                Options parsing               *
 * **********************************************/

/**
 * @brief Reads the configuration parameters of the application.
 * @param filename file the options data are read from
 * @param options_data options values (vcf filename, output directory...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 * 
 * Reads the basic configuration parameters of the application. If the configuration file can't be 
 * read, these parameters should be provided via the command-line interface.
 */
int read_shared_configuration(const char *filename, shared_options_t *options_data);


#endif
