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


list_item_t** create_chunks(list_t* records, int max_chunk_size, int *num_chunks) {
    *num_chunks = (int) ceil((float) records->length / max_chunk_size);
    LOG_DEBUG_F("%d chunks of %d elements top\n", *num_chunks, max_chunk_size);
    
    list_item_t **chunk_starts = (list_item_t**) malloc ((*num_chunks) * sizeof(list_item_t*));
    list_item_t *current = records->first_p;
    for (int j = 0; j < *num_chunks; j++) {
        chunk_starts[j] = current;
        for (int k = 0; k < max_chunk_size && current->next_p != NULL; k++) {
            current = current->next_p;
        }
    }

    return chunk_starts;
}
