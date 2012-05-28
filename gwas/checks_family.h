#ifndef CHECKS_FAMILY_H
#define CHECKS_FAMILY_H

#include <string.h>

#include <family.h>
#include <region.h>

int check_mendel(char *chromosome, int father_allele1, int father_allele2, int mother_allele1, int mother_allele2, 
                 int child_allele1, int child_allele2, enum Sex child_sex);

#endif
