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

#include "hpg_variant_utils.h"


/* ***********************
 *     Initialization    *
 * ***********************/

//char *find_configuration_file(int argc, char *argv[]) {
//    FILE *config_file = NULL;
//    char *config_filepath = NULL;
//    for (int i = 0; i < argc-1; i++) {
//        if (!strcmp("--config", argv[i])) {
//            config_filepath = argv[i+1];
//        }
//    }
//    if (!config_filepath) {
//        config_filepath = "hpg-variant.conf";
//    }
//
//    config_file = fopen(config_filepath, "r");
//    if (!config_file) {
//        LOG_FATAL("Configuration file can't be loaded!");
//    } else {
//        fclose(config_file);
//    }
//
//    LOG_DEBUG_F("Configuration file = %s\n", config_filepath);
//    return config_filepath;
//}

array_list_t *get_configuration_search_paths(int argc, char *argv[]) {
	char *c = get_config_path_from_args(argc, argv);
	char *config_dirpaths[3] = { c, getcwd(NULL, 0), "/etc/hpg-variant" };

	char *h = get_config_home_folder(config_dirpaths, 3);

	return sort_config_paths_by_priority(c, h);
}

char *get_config_path_from_args(int argc, char *argv[]) {
    FILE *config_dir = NULL;
    char *config_dirpath = NULL;
    struct stat sb;

    for (int i = 0; i < argc-1; i++) {
        if (!strcmp("--config", argv[i])) {
            config_dirpath = argv[i+1];
            break;
        }
    }

    if (config_dirpath && stat(config_dirpath, &sb) == -1 && errno == ENOENT) {
    	LOG_WARN("The folder specified to store configuration files does not exist.");
    	return NULL;
	}

    if (config_dirpath && !S_ISDIR(sb.st_mode)) {
    	LOG_WARN("The path specified to store configuration files is not a folder.");
    	return NULL;

    }

    return config_dirpath;
}

char *get_config_home_folder(char *config_dirpaths[], int num_dirpaths) {
	if (!getenv("HOME")) {	// Job mode, no HOME directory
		return NULL;
	}

	FILE *home_config_dir = NULL;
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
		array_list_insert(config_arg_path, paths);
	}

	array_list_insert(getcwd(NULL, 0), paths);

	if (home_path) {
		array_list_insert(home_path, paths);
	}

	array_list_insert(strndup("/etc/hpg-variant", 16), paths);

	return paths;
}


/* **********************
 *     Config files     *
 * **********************/

char *retrieve_config_file(char *filename, array_list_t *paths_to_search) {
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
 *      Miscellaneous    *
 * ***********************/

void show_usage(char *tool, void **argtable, int num_arguments) {
    printf("Usage: %s", tool);
    arg_print_syntaxv(stdout, argtable, "\n");
    arg_print_glossary(stdout, argtable, " %-40s %s\n");
}

int get_output_files(shared_options_data_t *shared_options, FILE** passed_file, FILE** failed_file) {
    if (shared_options == NULL) {
        return 1;
    }
    
    char *prefix_filename, *passed_filename, *failed_filename;
    int prefix_filename_len = 0;
    int dirname_len = strlen(shared_options->output_directory);
    
    if (shared_options->output_filename != NULL && strlen(shared_options->output_filename) > 0) {
        prefix_filename = shared_options->output_filename;
    } else if (shared_options->chain != NULL) {
        prefix_filename = calloc(strlen(shared_options->vcf_filename), sizeof(char));
        get_filename_from_path(shared_options->vcf_filename, prefix_filename);
    } else {
        // Output files are not created
        return 0;
    }
    
    prefix_filename_len = strlen(prefix_filename);
    
    passed_filename = (char*) calloc (dirname_len + prefix_filename_len + 11, sizeof(char));
    sprintf(passed_filename, "%s/%s.filtered", shared_options->output_directory, prefix_filename);
    *passed_file = fopen(passed_filename, "w");

    failed_filename = (char*) calloc (dirname_len + prefix_filename_len + 11, sizeof(char));
    sprintf(failed_filename, "%s/%s.rejected", shared_options->output_directory, prefix_filename);
    *failed_file = fopen(failed_filename, "w");
    
    LOG_DEBUG_F("passed filename = %s\nfailed filename = %s\n", passed_filename, failed_filename);
    
    free(passed_filename);
    free(failed_filename);
    
    return 0;
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
