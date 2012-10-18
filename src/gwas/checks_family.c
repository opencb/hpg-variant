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

#include "checks_family.h"

int check_mendel(char *chromosome, int father_allele1, int father_allele2, int mother_allele1, int mother_allele2, 
                 int child_allele1, int child_allele2, enum Sex child_sex) {
    // Ignore if any allele is missing
    if (father_allele1 < 0 || father_allele2 < 0 || mother_allele1 < 0 || mother_allele2 < 0 || child_allele1 < 0 || child_allele2 < 0) {
        return -1;
    }
    
    // Ignore haploid chromosomes
    if (!strcasecmp(chromosome, "Y") || !strcasecmp(chromosome, "MT")) {
        return -2;
    }
    
    int mendel_type = 0;
    
    if (strcasecmp(chromosome, "X") || child_sex == FEMALE) {
        
        if ( (!child_allele1 && child_allele2) || 
             (child_allele1 && !child_allele2) ) {
            // KID = 01/10
            // 00x00 -> 01  (m1)
            // 11x11 -> 01  (m2)

            if ( (!father_allele1 && !father_allele2) &&
                 (!mother_allele1 && !mother_allele2) ) {
                mendel_type = 1;
            } else if ( (father_allele1 && father_allele2) &&
                        (mother_allele1 && mother_allele2) ) {
                mendel_type = 2;
            }
            
        } else if ( !child_allele1 && !child_allele2 ) {
            // KID = 00
            // 00x11 -> 00 (m3) P11->00
            // 01x11 -> 00 (m3)
            // ??x11 -> 00 (m3)

            // 11x00 -> 00 (m4) M11->00
            // 11x01 -> 00 (m4)
            // 11x?? -> 00 (m4)

            // 11x11 -> 00 (m5) P11+M11->00

            // Hom parent can't breed opposite hom child

            // rule = at least one '11' parent

            if ( (father_allele1 && father_allele2) ||
                 (mother_allele1 && mother_allele2) ) {
                
                if ( father_allele1 && father_allele2 &&
                     mother_allele1 && mother_allele2 ) {
                    mendel_type = 5;
                } else if ( father_allele1 && father_allele2 ) {
                    mendel_type = 4;
                } else {
                    mendel_type = 3;
                }
            }
            
        } else {
            // KID = 11

            // 00x01 -> 11 (m6)
            // 00x11 -> 11
            // 00x?? -> 11

            // 01x00 -> 11 (m7)
            // 11x00 -> 11
            // ??x00 -> 11

            // 00x00 -> 11 (m8) P00+M00->11

            // rule = at least one '00' parent

            if ( (!father_allele1 && !father_allele2) ||
                 (!mother_allele1 && !mother_allele2) ) {
                
                if ( !father_allele1 && !father_allele2 &&
                     !mother_allele1 && !mother_allele2 ) {
                    mendel_type = 8;
                } else if ( !father_allele1 && !father_allele2 ) {
                    mendel_type = 6;
                } else {
                    mendel_type = 7;
                }
            }

        }

    } else {
        if ( child_allele1 && child_allele2 &&
             !mother_allele1 && !mother_allele2 ) {
            mendel_type = 9;
        }

        if ( !child_allele1 && !child_allele2 &&
             mother_allele1 && mother_allele2 ) {
            mendel_type = 10;
        }

    }
        
    return mendel_type;
}
