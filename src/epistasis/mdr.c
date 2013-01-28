#include "mdr.h"

bool mdr_high_risk_combinations(unsigned int count_affected, unsigned int count_unaffected, 
                                unsigned int samples_affected, unsigned int samples_unaffected, void **aux_return_values) {
    // Simplest MDR
    // return count_affected > count_unaffected;
    
    if (count_affected == 0 && count_unaffected == 0) {
        return false;
    }
    
    // Normalization applied (from MDR source code, file Model.java)
    int total_in_cell = count_affected + count_unaffected;
    double affected_unaffected_ratio = (double) samples_affected / samples_unaffected;
    double proportional_unaffected = count_unaffected * affected_unaffected_ratio;
    double reduction_ratio = total_in_cell / (proportional_unaffected + count_affected);
    
    double normalized_unaffected = proportional_unaffected * reduction_ratio;
    double normalized_affected = total_in_cell - normalized_unaffected;
    
    LOG_DEBUG_F("affected = %d/%d, unaffected = %d/%d -> %f\n", count_affected, samples_affected, 
                count_unaffected, samples_unaffected, normalized_affected - normalized_unaffected);
    return normalized_affected - normalized_unaffected > 1;
}



int** get_k_folds(unsigned int num_samples_affected, unsigned int num_samples_unaffected, unsigned int k, unsigned int **sizes) {
    if (num_samples_affected < k) {
        LOG_WARN("There are less affected samples than folds and they won't be properly distributed");
    }
    if (num_samples_unaffected < k) {
        LOG_WARN("There are less unaffected samples than folds and they won't be properly distributed");
    }
    
    // Fill an array of samples identifiers
    unsigned int num_samples = num_samples_affected + num_samples_unaffected;
    int samples[num_samples];
    for (int i = 0; i < num_samples; i++) {
        samples[i] = i;
    }
    
    // Shuffle affected and unaffected samples separately in order to guarantee they are not mixed
    array_shuffle_int(samples, num_samples_affected);
    array_shuffle_int(samples + num_samples_affected, num_samples_unaffected);
    
    // Get size of each fold (total, affected, unaffected) and initialize them
    int **folds = malloc (k * sizeof(unsigned int*));
    unsigned int *fold_sizes = calloc(3 * k, sizeof(unsigned int));
    unsigned int fold_size = num_samples / k;
    unsigned int mod_folds_size = num_samples % k;
    unsigned int all_folds_size = fold_size * k;
    for (int i = 0; i < k; i++) {
        // If the fold size can't be exactly the same for all folds, make distribution a bit uneven
        fold_sizes[3 * i] = (i < mod_folds_size) ? fold_size + 1 : fold_size ;
    }
    
    // Temporary buffers for storing affected and unaffected samples,
    // that will be merged so they are properly grouped for counting
    int total_affected_assigned = 0, total_unaffected_assigned = 0;
    int *samples_affected[k];
    int *samples_unaffected[k];
    for (int i = 0; i < k; i++) {
        samples_affected[i] = malloc(fold_sizes[3 * i] * sizeof(int));
        samples_unaffected[i] = malloc(fold_sizes[3 * i] * sizeof(int));
    }
    
    // Fill k-folds array
    // While the number of assigned samples is less than the total, assign one to each fold
    int i = 0;
    while (total_affected_assigned + total_unaffected_assigned < num_samples) {
        for (i = 0; i < k && total_affected_assigned + total_unaffected_assigned < num_samples; i++) {
            int my_offset_aff = fold_sizes[3 * i + 1];
            int my_offset_unaff = fold_sizes[3 * i + 2];
            
            if (total_affected_assigned < num_samples_affected) {
                samples_affected[i][my_offset_aff] = samples[total_affected_assigned];
                (fold_sizes[3 * i + 1])++ ;
                total_affected_assigned++ ;
            }
            
            if (total_unaffected_assigned < num_samples_unaffected) {
                samples_unaffected[i][my_offset_unaff] = samples[num_samples_affected + total_unaffected_assigned];
                (fold_sizes[3 * i + 2])++ ;
                total_unaffected_assigned++ ;
            }
        }
        
        LOG_DEBUG_F("total [%u] affected assigned = %u\tunaffected assigned = %u\n", i, total_affected_assigned, total_unaffected_assigned);
        
        i = (i < k) ? i + 1 : 0;
    }
    
    // Adjust to true values, compacting the array and removing -1 values
    for (int i = 0; i < k; i++) {
        fold_sizes[3 * i] = fold_sizes[3 * i + 1] + fold_sizes[3 * i + 2];
        folds[i] = malloc(fold_sizes[3 * i] * sizeof(int));
        
        // Insert affected samples
        memcpy(folds[i], samples_affected[i], fold_sizes[3 * i + 1] * sizeof(int));
        // Insert unaffected samples
        memcpy(folds[i] + fold_sizes[3 * i + 1], samples_unaffected[i], fold_sizes[3 * i + 2] * sizeof(int));
        // Sort them
        qsort(folds[i], fold_sizes[3 * i], sizeof(int), compare_int);
    }
    
    // Free temporary data
    for (int i = 0; i < k; i++) {
        free(samples_affected[i]);
        free(samples_unaffected[i]);
    }
    
//     printf("FOLDS:\n");
//     for (int i = 0; i < k; i++) {
//         printf("[%d] size (%d,%d,%d) -> ", i, fold_sizes[3 * i], fold_sizes[3 * i + 1], fold_sizes[3 * i + 2]);
//         
//         for (int j = 0; j < fold_sizes[3 * i]; j++) {
//             printf("%d ", folds[i][j]);
//         }
//         printf("\n");
//     }
    
    *sizes = fold_sizes;
    return folds;
}
