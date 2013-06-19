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

#ifndef HPG_VARIANT_ERROR_H
#define HPG_VARIANT_ERROR_H

// General errors (range 1-100 available)
#define NOT_IMPLEMENTED_TOOL                    1
#define CANT_READ_CONFIG_FILE                   2

#define VCF_FILE_NOT_SPECIFIED                  10
#define PED_FILE_NOT_SPECIFIED                  11

#define BATCH_SIZE_NOT_SPECIFIED                20

#define HOST_URL_NOT_SPECIFIED                  30
#define VERSION_NOT_SPECIFIED                   31
#define SPECIES_NOT_SPECIFIED                   32

// Effect tool errors
#define EFFECT_REGIONS_NOT_SPECIFIED            103

// GWAS tool errors
#define GWAS_TASK_NOT_SPECIFIED                 200
#define GWAS_MANY_TASKS_SPECIFIED               201

#define EPISTASIS_DATASET_NOT_SPECIFIED         210
#define EPISTASIS_ORDER_NOT_SPECIFIED           211
#define EPISTASIS_FOLDS_NOT_SPECIFIED           212
#define EPISTASIS_CV_RUNS_NOT_SPECIFIED         213
#define EPISTASIS_EVAL_MODE_NOT_SPECIFIED       214
#define EPISTASIS_OPS_PER_BLOCK_NOT_SPECIFIED   215

// VCF tools errors
// -- Filter tool errors
#define EMPTY_LIST_OF_FILTERS                   300

// -- Merge tool errors
#define MISSING_MODE_NOT_SPECIFIED              400
#define INFO_FIELDS_NOT_SPECIFIED               401
#define DISCORDANT_CHROMOSOME                   410
#define DISCORDANT_POSITION                     411
#define DISCORDANT_REFERENCE                    412

// -- Split tool errors
#define CRITERION_NOT_SPECIFIED                 500
#define INTERVALS_NOT_SPECIFIED                 501
// -- Stats tool errors

#endif
