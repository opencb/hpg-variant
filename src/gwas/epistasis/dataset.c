#include "dataset.h"

#define NUM_GENOTYPES   3


/* ***************************
 *  Whole dataset management *
 * ***************************/

#ifdef _USE_MPI

uint8_t *epistasis_dataset_load_mpi(char *filename, int *num_affected, int *num_unaffected, size_t *num_variants, 
                                    size_t *file_len, size_t *genotypes_offset, MPI_File *fd) {
    MPI_Offset len;
    MPI_Status status;
    
    // Check if file exists
    FILE *fp = fopen(filename, "rb");
    if (!fp) { return NULL; }
    
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fd);
    
    MPI_File_get_size(*fd, &len);
    
    LOG_DEBUG_F("File %s length = %llu bytes / %zu bytes\n", filename, len);
    
    uint8_t *map = malloc(len * sizeof(uint8_t));
    
    MPI_File_seek(*fd, 0, MPI_SEEK_SET);
    MPI_Barrier(MPI_COMM_WORLD);

    *genotypes_offset = sizeof(size_t) + sizeof(uint32_t) + sizeof(uint32_t);
    
    MPI_File_read(*fd, map, len, MPI_BYTE, &status);
    
    size_t *map_size_t = (size_t*) map;
    uint32_t *map_uint32 = (uint32_t*) (map + sizeof(size_t));
    *num_variants = map_size_t[0];
    *num_affected = map_uint32[0];
    *num_unaffected = map_uint32[1];
    LOG_DEBUG_F("num variants = %zu, aff = %d, unaff = %d\n", *num_variants, *num_affected, *num_unaffected);
    
    *file_len = len;
    
    return map;
}

void epistasis_dataset_close_mpi(uint8_t *contents, MPI_File fd) {
    free(contents);
    MPI_File_close(&fd);
}

#else

uint8_t *epistasis_dataset_load(int *num_affected, int *num_unaffected, size_t *num_variants, size_t *file_len, size_t *genotypes_offset, char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) { return NULL; }
    
    fread(num_variants, 1, sizeof(size_t), fp);
    fread(num_affected, 1, sizeof(uint32_t), fp);
    fread(num_unaffected, 1, sizeof(uint32_t), fp);
    fclose(fp);
    
    *genotypes_offset = sizeof(size_t) + sizeof(uint32_t) + sizeof(uint32_t);
    
    return mmap_file(file_len, filename);
}

int epistasis_dataset_close(uint8_t *contents, size_t file_len) {
    assert(contents);
    assert(file_len > 0);
    return munmap((void*) contents, file_len);
}

#endif

/* *********************************************
 *  Combinations of blocks, SNPs and genotypes *
 * *********************************************/

int get_block_stride(size_t block_operations, int order) {
    return ceil(pow(block_operations, ((double) 1/order)));
}

int get_next_block(int num_blocks, int order, int block_coordinates[order]) {
    // TODO not using next_combination because get_next_block tolerates repetition!
    for (int i = order - 1; i >= 0; i--) {
        if (block_coordinates[i] + 1 < num_blocks) {
            // Increment coordinate, for example: (0,1,2) -> (0,1,3)
            (block_coordinates[i])++;
            
            // The following coordinates must be, at least, the same as the current one
            // Let num_blocks=4, (0,1,3) --increment--> (0,2,3) --correct--> (0,2,2)
            if (i < order - 1) {
                for (int j = i + 1; j < order; j++) {
                    block_coordinates[j] = block_coordinates[i];
                }
            }
            
            return 1; // New valid block found
        }
    }
    
    return 0; // No more blocks available
}

void get_first_combination_in_block(int order, int init_coordinates[order], int block_coordinates[order], int stride) {
    init_coordinates[0] = block_coordinates[0] * stride;
    
    for (int i = 1; i < order; i++) {
        // Set values in base to the block coordinates and the stride of each block
        init_coordinates[i] = block_coordinates[i] * stride;
        
        // If two coordinates belong to the same block they must have consecutive positions
        // Let stride=100 and block_coordinates=(0,1,1), then combination (0,100,100) -> (0,100,101)
        if (init_coordinates[i] <= init_coordinates[i-1]) {
            init_coordinates[i] = init_coordinates[i-1] + 1;
        }
    }
}


/**
 * @brief  Generates the next combination of n elements as k after comb
 *
 * @param order the size of the subsets to generate
 * @param comb the previous combination ( use (0, 1, 2, ..., k) for first)
 * @param block_coordinates 
 * @param stride the size of the original set
 * @return 1 if a valid combination was found, 0 otherwise
 **/
int get_next_combination_in_block(int order, int comb[order], int block_coordinates[order], int stride, int num_variants) {
    int i = order - 1;
    ++comb[i];
    
    // comb[i] compared against last position allowed in its block,
    // or in the whole dataset, if the last index in the block exceeds the number of variants
    while (i > 0 && comb[i] >= MIN((block_coordinates[i] + 1) * stride - order + 1 + i, num_variants)) {
//         printf("** comb[i] = %d\tlimit = %d\n", comb[i], ((block_coordinates[i] + 1) * stride) - order + 1 + i);
        --i;
        ++comb[i];
    }
    
    /* Combination (n-k, n-k+1, ..., n) reached */
    if (comb[0] > (block_coordinates[0] + 1) * stride - 1 || comb[0] >= num_variants) {
//         printf("-- comb[i] = %d\tlimit = %d\n", comb[i], ((block_coordinates[i] + 1) * stride) - order + 1 + i);
        return 0; /* No more combinations can be generated */
    }

    /* comb now looks like (..., x, n, n, n, ..., n).
    Turn it into (..., x, x + 1, x + 2, ...) */
//     for (i = i + 1; i < k; ++i) {
//         comb[i] = comb[i - 1] + 1;
//         print_combination(comb, i, k);
//     }
    for (i = i + 1; i < order; ++i) {
        if (block_coordinates[i-1] == block_coordinates[i]) {
            comb[i] = comb[i - 1] + 1;
        } else {
            comb[i] = block_coordinates[i] * stride;
        }
    }
    
    if (comb[order - 1] > (block_coordinates[order - 1] + 1) * stride - 1 || comb[order - 1] >= num_variants) {
//         printf("++ comb[i] = %d\tlimit = %d\n", comb[order - 1], ((block_coordinates[order - 1] + 1) * stride) - order + 1 + i);
        return 0; /* No more combinations can be generated */
    }

//     print_combination(comb, 1, order);

    return 1;
}

uint8_t **get_genotype_combinations(int order, int *num_combinations) {
    *num_combinations = pow(NUM_GENOTYPES, order);
    uint8_t **combinations = malloc (*num_combinations * sizeof(uint8_t*));
    combinations[0] = calloc(order, sizeof(uint8_t));
    
    int has_next = 1;
    for (int i = 1; i < *num_combinations && has_next; i++) {
        combinations[i] = malloc(order * sizeof(uint8_t));
        memcpy(combinations[i], combinations[i-1], order * sizeof(uint8_t));
        has_next = get_next_genotype_combination(order, combinations[i]);
    }
    
    return combinations;
}

uint8_t get_next_genotype_combination(int order, uint8_t comb[order]) {
    for (int i = order - 1; i >= 0; i--) {
        ++comb[i];
        if (comb[i] < NUM_GENOTYPES) {
            return 1;
        } else if (comb[0] >= NUM_GENOTYPES) {
            return 0;
        } else {
            comb[i] = 0;
        }
    }
    
    return comb[0] < NUM_GENOTYPES;
}



/* ***************************
 *        Input/Output       *
 * ***************************/

void print_combination(int comb[], unsigned long idx, int k) {
    printf("%lu -> {", idx);
    int i;
    for (i = 0; i < k; ++i)
//         printf("%d, ", comb[i] + 1);
        printf("%d, ", comb[i]);
    printf("\b\b}\n");
}

void print_gt_combination(uint8_t comb[], unsigned long idx, int k) {
    printf("%lu -> {", idx);
    int i;
    for (i = 0; i < k; ++i)
//         printf("%d, ", comb[i] + 1);
        printf("%d, ", comb[i]);
    printf("\b\b}\n");
}
