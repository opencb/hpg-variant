#ifndef TRANSMISSION_DISEQUILIBRIUM_TEST_H
#define TRANSMISSION_DISEQUILIBRIUM_TEST_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cprops/hashtable.h>
#include <omp.h>

#include <family.h>
#include <ped_file.h>
#include <ped_file_structure.h>
#include <string_utils.h>
#include <vcf_batch.h>
#include <vcf_file.h>
#include <vcf_filters.h>
#include <vcf_util.h>

#include "global_options.h"
#include "gwas.h"
#include "hpg_variant_utils.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

typedef struct {
    char *chromosome;
    char *reference;
    char *alternate;
    
    unsigned long int position;
    
    int t1;
    int t2;
    double o_range;
    double chi_square;
} tdt_result_t;


int run_tdt_test(global_options_data_t *global_options_data, gwas_options_data_t *options_data);

int tdt_test(ped_file_t *ped_file, list_item_t *variants, int num_variants, cp_hashtable *sample_ids, list_t *output_list);


tdt_result_t *tdt_result_new(char *chromosome, unsigned long int position, char *reference, char *alternate, double t1, double t2, double chi_square);

void tdt_result_free(tdt_result_t *result);


static cp_hashtable *associate_samples_and_positions(vcf_file_t *file);

#endif
