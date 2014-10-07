/*
 * Copyright (c) 2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2013 Ignacio Medina (ICM-CIPF)
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
