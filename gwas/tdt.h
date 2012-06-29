#ifndef TRANSMISSION_DISEQUILIBRIUM_TEST_H
#define TRANSMISSION_DISEQUILIBRIUM_TEST_H

#include <stdlib.h>
#include <string.h>

#include <cprops/hashtable.h>
#include <cprops/linked_list.h>
#include <gsl/gsl_cdf.h>
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
    
    int t1;
    int t2;
    double odds_ratio;
    double chi_square;
    double p_value;
} tdt_result_t;

// int tdt_test(ped_file_t *ped_file, list_item_t *variants, int num_variants, cp_hashtable *sample_ids, list_t *output_list);
int tdt_test(ped_file_t *ped_file, vcf_record_t **variants, int num_variants, cp_hashtable *sample_ids, list_t *output_list);


tdt_result_t *tdt_result_new(char *chromosome, unsigned long int position, char *reference, char *alternate, double t1, double t2, double chi_square);

void tdt_result_free(tdt_result_t *result);

#endif
