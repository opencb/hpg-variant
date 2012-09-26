/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (CGI-CIPF)
 * Copyright (c) 2012 Ignacio Medina (CGI-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
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

// Effect tool errors
#define EFFECT_HOST_URL_NOT_SPECIFIED           100
#define EFFECT_VERSION_NOT_SPECIFIED            101
#define EFFECT_SPECIES_NOT_SPECIFIED            102
#define EFFECT_REGIONS_NOT_SPECIFIED            103


// GWAS tool errors
#define GWAS_TASK_NOT_SPECIFIED                 200
#define GWAS_MANY_TASKS_SPECIFIED               201

#endif
