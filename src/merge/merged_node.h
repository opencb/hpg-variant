#ifndef MERGED_NODE
#define MERGED_NODE

#include <stdlib.h>

#include <bioformats/features/region/region.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <containers/array_list.h>


typedef struct {
    region_t *region;
    array_list_t *links;
    
    char **chromosome_order;
    int num_chromosomes;
} merged_node;

typedef struct {
    vcf_record_t *record;
    vcf_file_t *file;
} vcf_record_file_link;



merged_node* merged_node_new(char* chromosome, long int position, char **chromosome_order, int num_chromosomes);

void merged_node_free(merged_node *node);


vcf_record_file_link *vcf_record_file_link_new(vcf_record_t *record, vcf_file_t *file);

void vcf_record_file_link_free(vcf_record_file_link *link);
// void vcf_record_file_link_free(vcf_record_file_link *link, int free_record, int free_file);


int merged_node_add_link(vcf_record_t *record, vcf_file_t *file, merged_node *node);

vcf_record_file_link *merged_node_remove_link(vcf_record_file_link *link, merged_node* node);


void set_cmp_criteria(char **chromosome_order, int num_chromosomes, merged_node *node);

int merged_node_cmp(merged_node *node1, merged_node *node2);


#endif
