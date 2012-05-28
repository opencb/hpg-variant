#ifndef TRANSMISSION_DISEQUILIBRIUM_TEST_H
#define TRANSMISSION_DISEQUILIBRIUM_TEST_H

#include <stdlib.h>
#include <string.h>

#include <cprops/hashtable.h>
#include <cprops/linked_list.h>
#include <omp.h>

#include <list.h>
#include <log.h>
#include <ped_file.h>
#include <ped_file_structure.h>
#include <vcf_file.h>
#include <vcf_file_structure.h>
#include <vcf_util.h>

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
} tdt_result_t;

int tdt_test(ped_file_t *ped_file, list_item_t *variants, int num_variants, cp_hashtable *sample_ids, list_t *output_list);


tdt_result_t *tdt_result_new(char *chromosome, unsigned long int position, char *reference, char *alternate, double t1, double t2, double chi_square);

void tdt_result_free(tdt_result_t *result);

#endif
