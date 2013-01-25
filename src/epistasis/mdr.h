#ifndef EPISTASIS_MDR
#define EPISTASIS_MDR

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include <commons/log.h>
#include <data/array_utils.h>

bool mdr_high_risk_combinations(unsigned int count_affected, unsigned int count_unaffected, 
                                unsigned int samples_affected, unsigned int samples_unaffected, void **aux_return_values);

int** get_k_folds(unsigned int samples_affected, unsigned int samples_unaffected, unsigned int k, unsigned int **sizes);

int* compact_array(int *array, size_t n);

#endif
