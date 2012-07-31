#ifndef ASSOCIATION_FISHER_TEST_H
#define ASSOCIATION_FISHER_TEST_H

#include <stdlib.h>
#include <string.h>

// #include <cprops/hashtable.h>
// #include <cprops/linked_list.h>
// #include <omp.h>
// 
// #include <bioformats/ped/ped_file.h>
// #include <bioformats/ped/ped_file_structure.h>
// #include <bioformats/vcf/vcf_file.h>
// #include <bioformats/vcf/vcf_file_structure.h>
// #include <bioformats/vcf/vcf_util.h>
// #include <commons/log.h>
// #include <containers/list.h>
#include <math/fisher.h>
// 
// #include "checks_family.h"

typedef struct {
    char *chromosome;
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

assoc_fisher_result_t *assoc_fisher_result_new(char *chromosome, int chromosome_len, unsigned long int position, 
                                               char *reference, int reference_len, char *alternate, int alternate_len, 
                                               int affected1, int affected2, int unaffected1, int unaffected2, double p_value);

void assoc_fisher_result_free(assoc_fisher_result_t *result);

#endif
