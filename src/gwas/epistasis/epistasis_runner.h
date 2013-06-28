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

#ifndef EPISTASIS_RUNNER_H
#define EPISTASIS_RUNNER_H

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>
#include <omp.h>

#include <commons/log.h>
#include <containers/heap.h>
#include <containers/khash.h>

#include "shared_options.h"
#include "hpg_variant_utils.h"

#include "cross_validation.h"
#include "dataset.h"
#include "epistasis.h"
#include "model.h"
#include "mpi/mpi_epistasis_helper.h"

int run_epistasis(shared_options_data_t *global_options_data, epistasis_options_data_t* options_data);

#endif
