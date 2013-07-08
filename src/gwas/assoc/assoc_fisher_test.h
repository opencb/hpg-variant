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

#ifndef ASSOCIATION_FISHER_TEST_H
#define ASSOCIATION_FISHER_TEST_H

#include <stdlib.h>
#include <string.h>

#include <stats/fisher.h>

typedef struct {
    char *chromosome;
    char *id;
    char *reference;
    char *alternate;
    
    unsigned long int position;
    
    int affected1;
    int affected2;
    int unaffected1;
    int unaffected2;
    
    double odds_ratio;
    double p_value;
} assoc_fisher_result_t;

double assoc_fisher_test(int a, int b, int c, int d, double *factorial_logarithms);

assoc_fisher_result_t *assoc_fisher_result_new(char *chromosome, int chromosome_len, unsigned long int position,  char *id, int id_len,
                                               char *reference, int reference_len, char *alternate, int alternate_len, 
                                               int affected1, int affected2, int unaffected1, int unaffected2, double p_value);

void assoc_fisher_result_free(assoc_fisher_result_t *result);

#endif
