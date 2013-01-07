#include "dataset.h"

// epistasis_dataset* epistasis_dataset_new() {
//     epistasis_dataset *dataset = malloc (sizeof(epistasis_dataset));
//     dataset->num_affected = 0;
//     dataset->num_unaffected = 0;
//     dataset->genotype_counts = array_list_new(1000, 1.5, COLLECTION_MODE_SYNCHRONIZED);
//     return dataset;
// }
// 
// void epistasis_dataset_free(epistasis_dataset *dataset) {
//     array_list_free(dataset->genotype_counts, free);
//     free(dataset);
// }

uint8_t *epistasis_dataset_process_records(vcf_record_t **variants, size_t num_variants, int *destination, int num_samples) {
    uint8_t *genotypes = malloc (num_variants * num_samples * sizeof(uint8_t));
    for (int i = 0; i < num_variants; i++) {
        vcf_record_t *record = variants[i];
        int gt_position = get_field_position_in_format("GT", strndup(record->format, record->format_len));
        
        // For each sample, get genotype and increment counters in the dataset
//         bool missing_found = false;
        for (int k = 0; k < num_samples; k++) {
//             printf("%d\tbase = %d\tindex = %d\n", k, destination[k], i * num_samples + destination[k]);
            char *sample = strdup(array_list_get(k, record->samples));
            int allele1, allele2, gt_dataset_index;
            if (get_alleles(sample, gt_position, &allele1, &allele2)) {
//                 if (!missing_found) {
//                     LOG_INFO_F("Missing samples found in %.*s:%ld (will be ignored)\n",
//                                 record->chromosome_len, record->chromosome, record->position);
//                     missing_found = true;
//                 }
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

// void epistasis_dataset_add_entry(uint16_t *genotype_count, epistasis_dataset *dataset) {
//     array_list_insert(genotype_count, dataset->genotype_counts);    
// }
// 
// 
// size_t epistasis_dataset_get_num_variants(epistasis_dataset *dataset) {
//     return dataset->genotype_counts->size;
// }
// 
// uint8_t *epistasis_dataset_get_variant_counts(size_t index, epistasis_dataset *dataset) {
//     return array_list_get(index, dataset->genotype_counts);
// }




int get_block_stride(size_t block_operations, int order) {
    return ceil(pow(block_operations, ((double) 1/order)));
}

int get_next_block(int num_blocks, int order, int block_coordinates[order]) {
    for (int i = order - 1; i >= 0; i--) {
        if (block_coordinates[i] + 1 < num_blocks) {
            // Increment coordinate, for example: (0,1,2) -> (0,1,3)
            (block_coordinates[i])++;
            
            // The following coordinates must be, at least, the same as the current one
            // Let num_blocks=4, (0,1,3) -> (0,2,3) -> (0,2,2)
            if (i < order - 1) {
                for (int j = i + 1; j < order; j++) {
                    block_coordinates[j] = block_coordinates[i];
                }
            }
            
            return 0; // New valid block found
        }
    }
    
    return 1; // No more blocks available
}

int* get_first_combination_in_block(int order, int block_coordinates[order], int stride) {
    int init_coordinates[order];
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
