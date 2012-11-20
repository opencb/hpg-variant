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

#ifndef MAIN_H
#define MAIN_H

/** 
 * @file main.h
 * @brief Entry point of the application
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <commons/log.h>
#include <containers/array_list.h>

#include "error.h"
#include "hpg_variant_utils.h"
#include "filter/filter.h"
#include "merge/merge.h"
#include "split/split.h"
#include "stats/stats.h"

int vcf_tool_filter(int argc, char *argv[], const char *configuration_file);

int vcf_tool_merge(int argc, char *argv[], const char *configuration_file, array_list_t *config_search_paths);

int vcf_tool_split(int argc, char *argv[], const char *configuration_file);

int vcf_tool_stats(int argc, char *argv[], const char *configuration_file);


#endif
