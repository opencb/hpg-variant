#ifndef EPISTASIS_MDR
#define EPISTASIS_MDR

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include <xmmintrin.h>
#include <smmintrin.h>

#include <commons/log.h>

bool mdr_high_risk_combinations(unsigned int count_affected, unsigned int count_unaffected, 
                                unsigned int samples_affected, unsigned int samples_unaffected, void **aux_return_values);

int *mdr_high_risk_combinations2(int *counts_affected, int *counts_unaffected, int num_counts,
                                 unsigned int num_affected, unsigned int num_unaffected,
                                 void **aux_return_values);
#endif
