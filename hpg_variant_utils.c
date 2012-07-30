#include "hpg_variant_utils.h"


char *find_configuration_file(int argc, char *argv[]) {
    FILE *config_file = NULL;
    char *config_filepath = NULL;
    for (int i = 0; i < argc-1; i++) {
        if (!strcmp("--config", argv[i])) {
            config_filepath = argv[i+1];
        }
    }
    if (!config_filepath) {
        config_filepath = "hpg-variant.cfg";
    }
    
    config_file = fopen(config_filepath, "r");
    if (!config_file) {
        LOG_FATAL("Configuration file can't be loaded!");
    } else {
        fclose(config_file);
    }
    
    LOG_INFO_F("Configuration file = %s\n", config_filepath);
    return config_filepath;
}

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



int get_output_files(shared_options_data_t *global_options_data, FILE** passed_file, FILE** failed_file) {
    if (global_options_data == NULL) {
        return 1;
    }
    
    char *prefix_filename, *passed_filename, *failed_filename;
    int prefix_filename_len = 0;
    int dirname_len = strlen(global_options_data->output_directory);
    
    if (global_options_data->output_filename != NULL && strlen(global_options_data->output_filename) > 0) {
        prefix_filename = global_options_data->output_filename;
    } else if (global_options_data->chain != NULL) {
        char input_filename[strlen(global_options_data->vcf_filename)];
        memset(input_filename, 0, strlen(global_options_data->vcf_filename) * sizeof(char));
        get_filename_from_path(global_options_data->vcf_filename, input_filename);
        prefix_filename = strdup(input_filename);
    } else {
        // Output files are not created
        return 0;
    }
    
    prefix_filename_len = strlen(prefix_filename);
    
    passed_filename = (char*) calloc (dirname_len + prefix_filename_len + 11, sizeof(char));
    sprintf(passed_filename, "%s/%s.filtered", global_options_data->output_directory, prefix_filename);
    *passed_file = fopen(passed_filename, "w");

    failed_filename = (char*) calloc (dirname_len + prefix_filename_len + 11, sizeof(char));
    sprintf(failed_filename, "%s/%s.rejected", global_options_data->output_directory, prefix_filename);
    *failed_file = fopen(failed_filename, "w");
    
//     LOG_INFO_F("passed filename = %s\nfailed filename = %s\n", passed_filename, failed_filename);
    
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
