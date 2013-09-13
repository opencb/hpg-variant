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

#include "hpg_variant_utils.h"


/* ***********************
 *     Initialization    *
 * ***********************/

array_list_t *get_configuration_search_paths(int argc, char *argv[]) {
    char *c = get_config_path_from_args(argc, argv);
    char *config_dirpaths[3] = { c, getcwd(NULL, 0), strdup("/etc/hpg-variant") };
    char *home = get_config_home_folder(config_dirpaths, 3);

    array_list_t *paths = sort_config_paths_by_priority(c, home);
    
    free(config_dirpaths[0]);
    free(config_dirpaths[1]);
    free(config_dirpaths[2]);
    free(home);
    
    return paths;
}

char *get_config_path_from_args(int argc, char *argv[]) {
    char *config_dirpath = NULL;
    struct stat sb;

    for (int i = 0; i < argc-1; i++) {
        if (!strcmp("--config", argv[i])) {
            config_dirpath = strdup(argv[i+1]);
            break;
        }
    }

    if (config_dirpath && stat(config_dirpath, &sb) == -1 && errno == ENOENT) {
        LOG_WARN("The folder specified to store configuration files does not exist.\n");
        return NULL;
    }

    if (config_dirpath && !S_ISDIR(sb.st_mode)) {
        LOG_WARN("The path specified to store configuration files is not a folder.\n");
        return NULL;

    }
    
    return config_dirpath;
}

char *get_config_home_folder(char *config_dirpaths[], int num_dirpaths) {
    if (!getenv("HOME")) {	// Job mode, no HOME directory
        return NULL;
    }

    char *home_config_dirpath = malloc (1024 * sizeof(char));
    char hpg_variant_conf_path_src[1024];
    char hpg_variant_conf_path_dest[1024];
    char vcf_info_fields_conf_path_src[1024];
    char vcf_info_fields_conf_path_dest[1024];

    struct stat sb;

    // Get home folder path
    sprintf(home_config_dirpath, "%s/.hpg-variant", getenv("HOME"));

    if (stat(home_config_dirpath, &sb) == -1 && errno == ENOENT) {
        // Create non-existing folder
        create_directory(home_config_dirpath);

        // Populate with the files from config_dirpaths
        FILE *from, *to;
        char ch;
        int hpg_variant_conf_copied = 0, vcf_info_fields_conf_copied = 0;

        for (int i = 0; i < num_dirpaths && !hpg_variant_conf_copied && !vcf_info_fields_conf_copied; i++) {
            if (!config_dirpaths[i]) {
                continue;
            }

            sprintf(hpg_variant_conf_path_src, "%s/hpg-variant.conf", config_dirpaths[i]);
            sprintf(vcf_info_fields_conf_path_src, "%s/vcf-info-fields.conf", config_dirpaths[i]);

            if (!hpg_variant_conf_copied && !stat(hpg_variant_conf_path_src, &sb)) {
                // Copy file hpg-variant.conf
                sprintf(hpg_variant_conf_path_dest, "%s/hpg-variant.conf", home_config_dirpath);
                from = fopen(hpg_variant_conf_path_src, "r");
                to = fopen(hpg_variant_conf_path_dest, "w");

                while((ch = fgetc(from)) != EOF) {
                    fputc(ch, to);
                }

                fclose(from);
                fclose(to);
                hpg_variant_conf_copied = 1;
            }

            if (!vcf_info_fields_conf_copied && !stat(vcf_info_fields_conf_path_src, &sb)) {
                // Copy file hpg-variant.conf
                sprintf(vcf_info_fields_conf_path_dest, "%s/vcf-info-fields.conf", home_config_dirpath);
                from = fopen(vcf_info_fields_conf_path_src, "r");
                to = fopen(vcf_info_fields_conf_path_dest, "w");

                while((ch = fgetc(from)) != EOF) {
                    fputc(ch, to);
                }

                fclose(from);
                fclose(to);
                vcf_info_fields_conf_copied = 1;
            }
        }
    }

    return home_config_dirpath;
}

array_list_t *sort_config_paths_by_priority(char *config_arg_path, char *home_path) {
    int num_paths = 2;
    num_paths += config_arg_path ? 1 : 0;
    num_paths += home_path ? 1 : 0;

    // Priority is: --config, current folder, home folder, /etc folder

    array_list_t *paths = array_list_new(num_paths, 1.1, COLLECTION_MODE_ASYNCHRONIZED);
    if (config_arg_path) {
        array_list_insert(strdup(config_arg_path), paths);
    }

    array_list_insert(getcwd(NULL, 0), paths);

    if (home_path) {
        array_list_insert(strdup(home_path), paths);
    }

    array_list_insert(strndup("/etc/hpg-variant", 16), paths);

    return paths;
}


/* **********************
 *     Config files     *
 * **********************/

char *retrieve_config_file(char *filename, array_list_t *paths_to_search) {
    assert(filename);
    assert(paths_to_search);
    
    char *filepath = NULL;
    struct stat sb;

    char aux_filepath[1024];
    for (int i = 0; i < paths_to_search->size; i++) {
        sprintf(aux_filepath, "%s/%s", (char*) array_list_get(i, paths_to_search), filename);
        if (!stat(aux_filepath, &sb)) {
            filepath = strdup(aux_filepath);
            break;
        }
    }

    return filepath;
}


/* **********************
 *    Job management    *
 * **********************/

FILE* new_job_status_file(char* path) {
    return fopen(path, "w");
}


void update_job_status_file(int percentage, FILE* file) {
    if (percentage < 100) {
        fprintf(file, "%d\n", percentage);
    } else {
        fprintf(file, "100\tDone\n");
    }
}

void close_job_status_file(FILE* file) {
    if (file) {
        fclose(file);
    }
}


/* ***********************
 *        Filtering      *
 * ***********************/

int get_filtering_output_files(shared_options_data_t *shared_options, FILE** passed_file, FILE** failed_file) {
    assert(shared_options);
    
    int filename_len = 0;
    int dirname_len = strlen(shared_options->output_directory);
    char prefix_filename[strlen(shared_options->vcf_filename)];
    
    if (shared_options->chain != NULL) {
        get_filename_from_path(shared_options->vcf_filename, prefix_filename);
        filename_len = strlen(prefix_filename);
    } else {
        // Output files are not created
        return 0;
    }
    
    LOG_DEBUG_F("prefix filename = %s\n", prefix_filename);
    
    char passed_filename[dirname_len + filename_len + 11];
    char failed_filename[dirname_len + filename_len + 11];
    
    sprintf(passed_filename, "%s/%s.filtered", shared_options->output_directory, prefix_filename);
    sprintf(failed_filename, "%s/%s.rejected", shared_options->output_directory, prefix_filename);
    
    *passed_file = fopen(passed_filename, "w");
    *failed_file = fopen(failed_filename, "w");
    
    LOG_DEBUG_F("passed filename = %s\nfailed filename = %s\n", passed_filename, failed_filename);
    
    return 0;
}

int write_filtering_output_files(array_list_t *passed_records, array_list_t *failed_records, FILE* passed_file, FILE* failed_file) {
    int ret_code = 0;
    if (passed_file) {
        if (passed_records != NULL && passed_records->size > 0) {
        #pragma omp critical 
            {
                for (int r = 0; r < passed_records->size; r++) {
                    ret_code |= write_vcf_record(passed_records->items[r], passed_file);
                }
            }
        }
    }
    
    if (failed_file) {
        if (failed_records != NULL && failed_records->size > 0) {
        #pragma omp critical 
            {
                for (int r = 0; r < failed_records->size; r++) {
                    ret_code |= write_vcf_record(failed_records->items[r], failed_file);
                }
            }
        }
    }
    
    return ret_code;
}

array_list_t *filter_records(filter_t** filters, int num_filters, individual_t **individuals, khash_t(ids) *sample_ids, int num_variables,
                             array_list_t *input_records, array_list_t **failed_records) {
    array_list_t *passed_records = NULL;
    if (filters == NULL || num_filters == 0) {
        passed_records = input_records;
    } else {
        *failed_records = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
        passed_records = run_filter_chain(input_records, *failed_records, individuals, sample_ids, num_variables, filters, num_filters);
    }
    return passed_records;
}

void free_filtered_records(array_list_t *passed_records, array_list_t *failed_records, array_list_t *input_records) {
    static int i = 0;
    // Free items in both lists (not their internal data)
    if (passed_records != input_records) {
        LOG_DEBUG_F("[Batch %d] %zu passed records\n", i, passed_records->size);
        array_list_free(passed_records, NULL);
    }
    if (failed_records) {
        LOG_DEBUG_F("[Batch %d] %zu failed records\n", i, failed_records->size);
        array_list_free(failed_records, NULL);
    }
    i++;
}


/* ***********************
 *         Output        *
 * ***********************/

FILE *get_output_file(shared_options_data_t *shared_options_data, char *default_name, char **path) {
    char *output_directory = (shared_options_data->output_directory && strlen(shared_options_data->output_directory) > 0) ? 
                              shared_options_data->output_directory : "." ;
    char *output_filename = (shared_options_data->output_filename && strlen(shared_options_data->output_filename) > 0) ? 
                             shared_options_data->output_filename : default_name;
    
    *path = (char*) malloc ((strlen(output_directory) + strlen(output_filename) + 2) * sizeof(char));
    sprintf(*path, "%s/%s", output_directory, output_filename);
    return fopen(*path, "w");
}


/* ***********************
 *      Miscellaneous    *
 * ***********************/

void show_usage(char *tool, void **argtable) {
    printf("Usage: %s", tool);
    arg_print_syntaxv(stdout, argtable, "\n");
    arg_print_glossary(stdout, argtable, " %-40s %s\n");
}

void show_version(char *tool) {
    printf("HPG Variant %s %s\nCopyright (C) 2013 Institute for Computational Medicine (CIPF).\n\
This is free software; see the source for copying conditions.\n", tool, HPG_VARIANT_VERSION);
}


int *create_chunks(int length, int max_chunk_size, int *num_chunks, int **chunk_sizes) {
    *num_chunks = (int) ceil((float) length / max_chunk_size);
    int *chunk_starts = (int*) calloc (*(num_chunks), sizeof(int));
    *chunk_sizes = (int*) calloc (*(num_chunks), sizeof(int));
    LOG_DEBUG_F("%d chunks of %d elements top\n", *num_chunks, max_chunk_size);
    
    for (int j = 0; j < *num_chunks; j++) {
        chunk_starts[j] = j * max_chunk_size;
        (*chunk_sizes)[j] = max_chunk_size;
    }
    int mod = length % max_chunk_size;
    if (mod > 0) {
        (*chunk_sizes)[*num_chunks-1] = mod;
    }

    return chunk_starts;
}
