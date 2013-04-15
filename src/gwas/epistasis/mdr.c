#include "mdr.h"

bool mdr_high_risk_combinations(unsigned int count_affected, unsigned int count_unaffected, 
                                unsigned int samples_affected, unsigned int samples_unaffected, void **aux_return_values) {
    // Simplest MDR
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
    
    LOG_DEBUG_F("affected = %d/%d, unaffected = %d/%d -> %f - %f = %f\n", count_affected, samples_affected, 
                count_unaffected, samples_unaffected, normalized_affected, normalized_unaffected, normalized_affected - normalized_unaffected);
    return normalized_affected >= normalized_unaffected;
}


int *mdr_high_risk_combinations2(int *counts_affected, int *counts_unaffected, int num_counts,
                                 unsigned int num_affected, unsigned int num_unaffected,
                                 void **aux_return_values) {
    int max_num_high_risk = 16 * (int) ceil(((double) num_counts) / 16);
    int *high_risk = _mm_malloc(max_num_high_risk * sizeof(int), 16);
    
    __m128 affected_unaffected_ratio = _mm_set1_ps((float) num_affected / num_unaffected);
    __m128 counts_affected_sse, counts_unaffected_sse, total_in_cell;
    __m128 proportional_unaffected, reduction_ratio;
    __m128 normalized_affected, normalized_unaffected;
    __m128i result;
    
    for (int i = 0; i < num_counts; i += 4) {
        counts_affected_sse = _mm_set_ps(counts_affected[i+3], counts_affected[i+2], counts_affected[i+1], counts_affected[i]);
        counts_unaffected_sse = _mm_set_ps(counts_unaffected[i+3], counts_unaffected[i+2], counts_unaffected[i+1], counts_unaffected[i]);
        
        total_in_cell = _mm_add_ps(counts_affected_sse, counts_unaffected_sse);
        
        proportional_unaffected = _mm_mul_ps(counts_unaffected_sse, affected_unaffected_ratio);
        reduction_ratio = _mm_div_ps(total_in_cell, _mm_add_ps(proportional_unaffected, counts_affected_sse));
        
        normalized_unaffected = _mm_mul_ps(proportional_unaffected, reduction_ratio);
        normalized_affected = _mm_sub_ps(total_in_cell, normalized_unaffected);
        
        result = _mm_cvtps_epi32(_mm_cmpge_ps(normalized_affected, normalized_unaffected));
        _mm_store_si128(high_risk + i, result);
    }
    
    return high_risk;
}

