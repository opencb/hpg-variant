#include "mdr.h"

bool mdr_high_risk_combinations(int count_affected, int count_unaffected, 
                                int samples_affected, int samples_unaffected, void **aux_return_values) {
    // Most simple MDR
    // return count_affected > count_unaffected;
    
    // Normalization applied (from MDR source code, file Model.java)
    int total_in_cell = count_affected + count_unaffected;
    double affected_unaffected_ratio = (double) samples_affected / samples_unaffected;
    double proportional_unaffected = count_unaffected * affected_unaffected_ratio;
    double reduction_ratio = total_in_cell / (proportional_unaffected + count_affected);
    
    double normalized_unaffected = proportional_unaffected * reduction_ratio;
    double normalized_affected = total_in_cell - normalized_unaffected;
    
    LOG_DEBUG_F("affected = %d/%d, unaffected =%d/%d -> %f\n", count_affected, samples_affected, 
                count_unaffected, samples_unaffected, normalized_affected - normalized_unaffected);
    return normalized_affected - normalized_unaffected > 1;
}
