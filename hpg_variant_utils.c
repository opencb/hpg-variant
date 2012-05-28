#include "hpg_variant_utils.h"


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
