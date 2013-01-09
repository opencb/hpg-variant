#ifndef EPISTASIS_MODEL
#define EPISTASIS_MODEL

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dataset.h"

int* get_counts(int order, uint8_t *genotypes, int num_affected, int num_unaffected, int *num_counts);

uint8_t* get_masks(int order, uint8_t *genotypes, int num_samples, int *num_masks);

#endif
