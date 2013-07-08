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

#include "assoc_fisher_test.h"


double assoc_fisher_test(int a, int b, int c, int d, double *factorial_logarithms) {
    return fisher_test(a, b, c, d, TWO_SIDED, factorial_logarithms);
}



assoc_fisher_result_t* assoc_fisher_result_new(char *chromosome, int chromosome_len, unsigned long int position, char *id, int id_len, 
                                               char *reference, int reference_len, char *alternate, int alternate_len, 
                                               int affected1, int affected2, int unaffected1, int unaffected2, double p_value) {
    assoc_fisher_result_t *result = (assoc_fisher_result_t*) malloc (sizeof(assoc_fisher_result_t));
    
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
    result->p_value = p_value;
    
    return result;
}

void assoc_fisher_result_free(assoc_fisher_result_t* result) {
    free(result->chromosome);
    free(result->reference);
    free(result->alternate);
    free(result);
}
