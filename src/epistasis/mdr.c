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

int** get_k_folds(unsigned int samples_affected, unsigned int samples_unaffected, unsigned int k, unsigned int **sizes) {
    if (samples_affected < k) {
        LOG_WARN("There are less affected samples than folds and they won't be properly distributed");
    }
    if (samples_unaffected < k) {
        LOG_WARN("There are less unaffected samples than folds and they won't be properly distributed");
    }
    
    // Fill an array of samples identifiers
    unsigned int num_samples = samples_affected + samples_unaffected;
//     int *samples = malloc (num_samples * sizeof(int));
    int samples[num_samples];
    for (int i = 0; i < num_samples; i++) {
        samples[i] = i;
    }
    
//     printf("before = { ");
//     for (int i = 0; i < samples_affected; i++) {
//         printf("%d ", samples[i]);
//     }
//     printf(" }\n");
    
    // Shuffle affected and unaffected samples separately in order to guarantee they are not mixed
    array_shuffle_int(samples, samples_affected);
    array_shuffle_int(samples + samples_affected, samples_unaffected);
    
//     printf("shuffled = { ");
//     for (int i = 0; i < num_samples; i++) {
//         printf("%d ", samples[i]);
//     }
//     printf(" }\n");
    
    // Get size of each fold (total, affected, unaffected)
    unsigned int *fold_sizes = malloc(3 * k * sizeof(unsigned int));
    unsigned int fold_size = num_samples / k;
    unsigned int mod_folds_size = num_samples % k;
    unsigned int all_folds_size = fold_size * k;
    for (int i = 0; i < k; i++) {
        // If the fold size can't be exactly the same for all folds, make distribution a bit uneven
        fold_sizes[3 * i] = (i < mod_folds_size) ? fold_size + 1 : fold_size ;
        fold_sizes[3 * i + 1] = 0 ;
        fold_sizes[3 * i + 2] = 0 ;
    }
    
//     printf("fold sizes = { ");
//     for (int i = 0; i < k; i++) {
//         printf("%u ", fold_sizes[3*i]);
//     }
//     printf(" }\n");
    
    // Fill k-folds array
    int **folds = malloc (k * sizeof(unsigned int*));
    
    unsigned int total_affected_assigned = 0, total_unaffected_assigned = 0;
    double dbl_affected_per_fold, dbl_unaffected_per_fold;
    unsigned int affected_per_fold, unaffected_per_fold;
    for (int i = 0; i < k; i++) {
        dbl_affected_per_fold = (double) samples_affected * fold_sizes[3 * i] / num_samples;
        dbl_unaffected_per_fold = (double) samples_unaffected * fold_sizes[3 * i] / num_samples;
        affected_per_fold = round(dbl_affected_per_fold);
        unaffected_per_fold = round(dbl_unaffected_per_fold);
        
        if (affected_per_fold + unaffected_per_fold < fold_sizes[3 * i]) {
            // If the number of samples assigned is too small, check which stratum is more 
            // reasonable to have rounded down, and increment the size of the other one
            if (fmod(dbl_affected_per_fold, 1.0) < fmod(dbl_unaffected_per_fold, 1.0)) {
                unaffected_per_fold++;
            } else {
                affected_per_fold++;
            }
        } else if (affected_per_fold + unaffected_per_fold > fold_sizes[3 * i]) {
            // If the number of samples assigned is too large, check which stratum is more 
            // reasonable to have rounded up, and decrement the size of the other one
            if (fmod(dbl_affected_per_fold, 1.0) > fmod(dbl_unaffected_per_fold, 1.0)) {
                unaffected_per_fold--;
            } else {
                affected_per_fold--;
            }
        }
        
//         printf("[%u] affected = %u\tunaffected = %u\n", i, affected_per_fold, unaffected_per_fold);
        
        // Assign affected and unaffected indices to the fold
        folds[i] = malloc((affected_per_fold + unaffected_per_fold) * sizeof(unsigned int));
        for (int j = 0; j < affected_per_fold + unaffected_per_fold; j++) {
            folds[i][j] = -1;
        }
        
        unsigned int fold_affected_assigned = 0, fold_unaffected_assigned = 0;
        
        for (int j = 0; j < affected_per_fold && total_affected_assigned + j < samples_affected; j++) {
            // Affected samples are formerly positioned in the samples array
            folds[i][j] = samples[total_affected_assigned + j];
            fold_affected_assigned++;
        }
        for (int j = 0; j < unaffected_per_fold && total_unaffected_assigned + j < samples_unaffected; j++) {
            // Unaffected samples are positioned after affected ones
            folds[i][affected_per_fold + j] = samples[samples_affected + total_unaffected_assigned + j];
            fold_unaffected_assigned++;
        }
        
        total_affected_assigned += fold_affected_assigned;
        total_unaffected_assigned += fold_unaffected_assigned;
        
        fold_sizes[3 * i + 1] = fold_affected_assigned;
        fold_sizes[3 * i + 2] = fold_unaffected_assigned;
        
//         printf("[%u] affected assigned = %u\tunaffected assigned = %u\n", i, fold_affected_assigned, fold_unaffected_assigned);
//         printf("total [%u] affected assigned = %u\tunaffected assigned = %u\n", i, total_affected_assigned, total_unaffected_assigned);
        
    }
    
    // Just in case the number of assigned un/affected doesn't match the total,
    // distribute the missing ones among the k-folds
    int i = 0;
    while (total_affected_assigned < samples_affected) {
        folds[i] = realloc(folds[i], (fold_sizes[3 * i] + 1) * sizeof(unsigned int));
        folds[i][fold_sizes[3 * i]] = samples[total_affected_assigned];
        
        (fold_sizes[3 * i])++;
        (fold_sizes[3 * i + 1])++;
        total_affected_assigned++;
        
//         printf("tmp%d -> affected assigned = %u\tunaffected assigned = %u\n", i, total_affected_assigned, total_unaffected_assigned);
        i = (i < k) ? i + 1 : 0;
    }
    
    i = 0;
    while (total_unaffected_assigned < samples_unaffected) {
        folds[i] = realloc(folds[i], (fold_sizes[3 * i] + 1) * sizeof(unsigned int));
        folds[i][fold_sizes[3 * i]] = samples[samples_affected + total_unaffected_assigned];
        
        (fold_sizes[3 * i])++;
        (fold_sizes[3 * i + 2])++;
        total_unaffected_assigned++;
        
        i = (i < k) ? i + 1 : 0;
    }
    
//     // Adjust to true values (TODO consider -1 values!)
//     for (int i = 0; i < k; i++) {
//         fold_sizes[3 * i] = fold_sizes[3 * i + 1] + fold_sizes[3 * i + 2];
//     }
    
    
//     printf("final -> affected assigned = %u\tunaffected assigned = %u\n", total_affected_assigned, total_unaffected_assigned);
//         
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
