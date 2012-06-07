#ifndef ASSOC_H
#define ASSOC_H

#include <cprops/hashtable.h>
#include <omp.h>

#include <bioformats/family.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/ped/ped_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_util.h>
#include <commons/log.h>
#include <containers/list.h>

#include "gwas.h"
#include "assoc_basic_test.h"
// #include "assoc_fisher_test.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

void assoc_test(enum GWAS_task test_type, ped_file_t *ped_file, list_item_t *variants, int num_variants, cp_hashtable *sample_ids, list_t *output_list);

void assoc_count_individual(individual_t *individual, vcf_record_t *record, int allele1, int allele2, 
                           int *affected1, int *affected2, int *unaffected1, int *unaffected2);

#endif
