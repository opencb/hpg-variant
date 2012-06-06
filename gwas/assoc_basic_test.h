#ifndef ASSOCIATION_BASIC_TEST_H
#define ASSOCIATION_BASIC_TEST_H

#include <stdlib.h>
#include <string.h>

#include <cprops/hashtable.h>
#include <cprops/linked_list.h>
#include <omp.h>

#include <bioformats/ped/ped_file.h>
#include <bioformats/ped/ped_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_util.h>
#include <commons/log.h>
#include <containers/list.h>

#include "checks_family.h"

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
    double chi_square;
} assoc_basic_result_t;

int assoc_basic_test(ped_file_t *ped_file, list_item_t *variants, int num_variants, cp_hashtable *sample_ids, list_t *output_list);
// int basic_assoc_test(list_item_t *variants, int num_variants, cp_hashtable *sample_ids, list_t *output_list);

// void basic_assoc_individual(individual_t *individual, vcf_record_t *record, int allele1, int allele2, 
//                            int *affected1, int *affected2, int *unaffected1, int *unaffected2, double *chi_square);
void assoc_basic_individual(individual_t *individual, vcf_record_t *record, int allele1, int allele2, 
                           int *affected1, int *affected2, int *unaffected1, int *unaffected2);

double chi_square(int a, int b, int c, int d);

assoc_basic_result_t *assoc_basic_result_new(char *chromosome, unsigned long int position, char *reference, char *alternate, 
                                             int affected1, int affected2, int unaffected1, int unaffected2, double chi_square);

void assoc_basic_result_free(assoc_basic_result_t *result);

#endif
