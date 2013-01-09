#ifndef EPISTASIS_MDR
#define EPISTASIS_MDR

#include <stdbool.h>

#include <commons/log.h>

bool mdr_high_risk_combinations(int count_affected, int count_unaffected, 
                                int samples_affected, int samples_unaffected, void **aux_return_values);



#endif
