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

#ifndef MPI_EPISTASIS_HELPER_H
#define	MPI_EPISTASIS_HELPER_H

#include <stdlib.h>

#include <mpi.h>

#include "../epistasis.h"
#include "../model.h"
#include "shared_options.h"

#define TAG_RANKING_RISKY_SIZE  0
#define TAG_RANKING_RISKY_ELEM  1

typedef struct {
    int lengths[2];
    MPI_Datatype types[2];
    MPI_Aint offsets[2];
    MPI_Datatype datatype;
} risky_combination_mpi_t;

void bcast_shared_options_data_mpi(shared_options_data_t *options_data, int root, MPI_Comm comm);

void bcast_epistasis_options_data_mpi(epistasis_options_data_t *options_data, int root, MPI_Comm comm);

void risky_combination_mpi_init(risky_combination_mpi_t *type);

void risky_combination_mpi_free(risky_combination_mpi_t *type);

void send_ranking_risky_size_mpi(struct heap *ranking, int dest, MPI_Comm comm);

size_t receive_ranking_risky_size_mpi(int src, MPI_Comm comm, MPI_Status stat);

void send_risky_combination_mpi(risky_combination *risky, risky_combination_mpi_t type, int dest, MPI_Comm comm);

risky_combination* receive_risky_combination_mpi(int order, risky_combination_mpi_t type, masks_info info, int src, MPI_Comm comm, MPI_Status stat);


#endif
