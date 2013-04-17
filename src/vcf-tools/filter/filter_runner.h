/*
 * Copyright (c) 2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
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

#ifndef FILTER_RUNNER_H
#define	FILTER_RUNNER_H

#include "filter.h"
#include "shared_options.h"


typedef struct {
    vcf_batch_t *batch;
    array_list_t *passed_records;
    array_list_t *failed_records;
} filter_temp_output_t;


/* ******************************
 *       Tool execution         *
 * ******************************/

int run_filter(shared_options_data_t *shared_options_data, filter_options_data_t *options_data);


/* ******************************
 *           Auxiliary          *
 * ******************************/

filter_temp_output_t *filter_temp_output_new(vcf_batch_t *batch, array_list_t *passed_records, array_list_t *failed_records);

void filter_temp_output_free(filter_temp_output_t *temp);

#endif
