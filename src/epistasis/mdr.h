#ifndef EPISTASIS_MDR
#define EPISTASIS_MDR

#include <stdbool.h>

#include <commons/log.h>

bool mdr_high_risk_combinations(unsigned int count_affected, unsigned int count_unaffected, 
                                unsigned int samples_affected, unsigned int samples_unaffected, void **aux_return_values);



#endif
