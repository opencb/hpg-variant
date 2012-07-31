#ifndef ASSOCIATION_BASIC_TEST_H
#define ASSOCIATION_BASIC_TEST_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_cdf.h>

typedef struct {
    char *chromosome;
    char *reference;
    char *alternate;
    
    unsigned long int position;
    
    int affected1;
    int affected2;
    int unaffected1;
    int unaffected2;
    
    double odds_ratio;
    double chi_square;
    double p_value;
} assoc_basic_result_t;

double assoc_basic_test(int a, int b, int c, int d);

assoc_basic_result_t *assoc_basic_result_new(char *chromosome, int chromosome_len, unsigned long int position, 
                                             char *reference, int reference_len, char *alternate, int alternate_len, 
                                             int affected1, int affected2, int unaffected1, int unaffected2, double chi_square);

void assoc_basic_result_free(assoc_basic_result_t *result);

#endif
