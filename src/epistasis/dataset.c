#include "dataset.h"

#define NUM_GENOTYPES   3


/* ***************************
 *  Whole dataset management *
 * ***************************/

uint8_t *epistasis_dataset_process_records(vcf_record_t **variants, size_t num_variants, int *destination,
                                           int num_samples, int threads) {
    uint8_t *genotypes = malloc (num_variants * num_samples * sizeof(uint8_t));
    
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < num_variants; i++) {
        vcf_record_t *record = variants[i];
        int gt_position = get_field_position_in_format("GT", strndup(record->format, record->format_len));
        
        // For each sample, get genotype and increment counters in the dataset
        for (int k = 0; k < num_samples; k++) {
//             printf("%d\tbase = %d\tindex = %d\n", k, destination[k], i * num_samples + destination[k]);
            char *sample = strdup(array_list_get(k, record->samples));
            int allele1, allele2, gt_dataset_index;
            if (get_alleles(sample, gt_position, &allele1, &allele2)) {
                genotypes[i * num_samples + destination[k]] = 255;
                
            } else {
                if (!allele1 && !allele2) { // Homozygous in first allele
                    genotypes[i * num_samples + destination[k]] = 0;
                } else if (allele1 != allele2) { // Heterozygous
                    genotypes[i * num_samples + destination[k]] = 1;
                } else if (allele1 && allele1 == allele2) { // Homozygous in second allele
                    genotypes[i * num_samples + destination[k]] = 2;
                }
            }
            free(sample);
        }
    }
    
    return genotypes;
}

uint8_t *epistasis_dataset_load(int *num_affected, int *num_unaffected, size_t *num_variants, size_t *file_len, size_t *genotypes_offset, char *filename) {
    uint8_t *contents = mmap_file(file_len, filename);
    
    *num_variants = contents[0];
    *num_affected = contents[sizeof(size_t)];
    *num_unaffected = contents[sizeof(size_t) + sizeof(uint32_t)];
    *genotypes_offset = sizeof(size_t) + sizeof(uint32_t) + sizeof(uint32_t);
    
    return contents;
}

int epistasis_dataset_close(uint8_t *contents, size_t file_len) {
    assert(contents);
    assert(file_len > 0);
    return munmap((void*) contents, file_len);
}



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

int* get_first_combination_in_block(int order, int block_coordinates[order], int stride) {
    int *init_coordinates = malloc(order * sizeof(int));
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
    
    return init_coordinates;
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
int get_next_combination_in_block(int order, int comb[order], int block_coordinates[order], int stride) {
    int i = order - 1;
    ++comb[i];
    
    // comb[i] compared against last position allowed in its block
    while ((i > 0) && (comb[i] >= ((block_coordinates[i] + 1) * stride) - order + 1 + i)) {
        --i;
        ++comb[i];
    }
    
    /* Combination (n-k, n-k+1, ..., n) reached */
    if (comb[0] > (block_coordinates[0] + 1) * stride - 1) {
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
    
    if (comb[order - 1] > (block_coordinates[order - 1] + 1) * stride - 1) {
        return 0; /* No more combinations can be generated */
    }

//     print_combination(comb, (unsigned long) 3000, order);

    return 1;
}

uint8_t **get_genotype_combinations(int order, int *num_combinations) {
    *num_combinations = pow(3, order);
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
