#ifndef HPG_VARIANT_UTIL_H
#define HPG_VARIANT_UTIL_H

#include <math.h>
#include <stdlib.h>

#include <list.h>
#include <log.h>


// inline int int_cmp(int *a, int *b) {
//     return *a - *b;
// }

/**
 * @brief Given a list of records, distributes them in chunks of similar size
 * @param records list of records to separate in chunks
 * @param max_chunk_size maximum size of a chunk
 * @param[out] num_chunks number of chunks created
 * @return The list of chunks
 * 
 * Given a list of records, defines another list of chunks whose elements point to records separated 
 * by a maximum distance of max_chunk_size. These records will mark the beginning of each chunk.
 */
list_item_t** create_chunks(list_t* records, int max_chunk_size, int *num_chunks);

#endif
