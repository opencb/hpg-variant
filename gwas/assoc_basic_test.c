#include "assoc_basic_test.h"

double assoc_basic_test(int a, int b, int c, int d) {
    double total_alleles = a + c + b + d;
    
    double total_affected = a + c;
    double total_unaffected = b + d;
    
    double total_allele1 = a + b;
    double total_allele2 = c + d;
    
    double expected_affected_allele1    = (total_affected   * total_allele1) / total_alleles;
    double expected_affected_allele2    = (total_affected   * total_allele2) / total_alleles;
    double expected_unaffected_allele1  = (total_unaffected * total_allele1) / total_alleles;
    double expected_unaffected_allele2  = (total_unaffected * total_allele2) / total_alleles;

    return ((a - expected_affected_allele1)   * (a - expected_affected_allele1))   / expected_affected_allele1 + 
           ((c - expected_affected_allele2)   * (c - expected_affected_allele2))   / expected_affected_allele2 +
           ((b - expected_unaffected_allele1) * (b - expected_unaffected_allele1)) / expected_unaffected_allele1 + 
           ((d - expected_unaffected_allele2) * (d - expected_unaffected_allele2)) / expected_unaffected_allele2 ;
}


assoc_basic_result_t* assoc_basic_result_new(char *chromosome, unsigned long int position, char *reference, char *alternate, 
                                             int affected1, int affected2, int unaffected1, int unaffected2, double chi_square) {
    assoc_basic_result_t *result = (assoc_basic_result_t*) malloc (sizeof(assoc_basic_result_t));
    
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
    result->chi_square = chi_square;
    
    return result;
}

void assoc_basic_result_free(assoc_basic_result_t* result) {
    free(result->chromosome);
    free(result->reference);
    free(result->alternate);
    free(result);
}
