/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * hpg-variant is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-variant. If not, see <http://www.gnu.org/licenses/>.
 */

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


assoc_basic_result_t* assoc_basic_result_new(char *chromosome, int chromosome_len, unsigned long int position, char *id, int id_len,
                                             char *reference, int reference_len, char *alternate, int alternate_len, 
                                             int affected1, int affected2, int unaffected1, int unaffected2, double chi_square) {
    assoc_basic_result_t *result = (assoc_basic_result_t*) malloc (sizeof(assoc_basic_result_t));
    
    result->chromosome = strndup(chromosome, chromosome_len);
    result->position = position;
    result->id = strndup(id, id_len);
    result->reference = strndup(reference, reference_len);
    result->alternate = strndup(alternate, alternate_len);
    result->affected1 = affected1;
    result->affected2 = affected2;
    result->unaffected1 = unaffected1;
    result->unaffected2 = unaffected2;
    result->odds_ratio = (affected2 == 0 || unaffected1 == 0) ? NAN : 
                         ((double) affected1 / affected2) * ((double) unaffected2 / unaffected1);
    result->chi_square = chi_square;
    result->p_value = 1 - gsl_cdf_chisq_P(chi_square, 1);
    
    return result;
}

void assoc_basic_result_free(assoc_basic_result_t* result) {
    free(result->chromosome);
    free(result->reference);
    free(result->alternate);
    free(result);
}
