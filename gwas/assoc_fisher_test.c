#include "assoc_fisher_test.h"


double assoc_fisher_test(int a, int b, int c, int d, double *factorial_logarithms) {
    return fisher_test(a, b, c, d, TWO_SIDED, factorial_logarithms);
}



assoc_fisher_result_t* assoc_fisher_result_new(char *chromosome, unsigned long int position, char *reference, char *alternate, 
                                             int affected1, int affected2, int unaffected1, int unaffected2, double p_value) {
    assoc_fisher_result_t *result = (assoc_fisher_result_t*) malloc (sizeof(assoc_fisher_result_t));
    
    result->chromosome = strdup(chromosome);
    result->position = position;
    result->reference = strdup(reference);
    result->alternate = strdup(alternate);
    result->affected1 = affected1;
    result->affected2 = affected2;
    result->unaffected1 = unaffected1;
    result->unaffected2 = unaffected2;
    result->odds_ratio = (affected2 == 0 || unaffected1 == 0) ? NAN : 
                         ((double) affected1 / affected2) * ((double) unaffected2 / unaffected1);
    result->p_value = p_value;
    
    return result;
}

void assoc_fisher_result_free(assoc_fisher_result_t* result) {
    free(result->chromosome);
    free(result->reference);
    free(result->alternate);
    free(result);
}
