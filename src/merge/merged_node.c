#include "merged_node.h"

merged_node* merged_node_new(char* chromosome, long int position, char **chromosome_order, int num_chromosomes) {
    assert(chromosome);
    assert(chromosome_order);
    merged_node *node = malloc (sizeof(merged_node));
    node->region = region_new(chromosome, position, position);
    node->links = array_list_new(8, 1.5, COLLECTION_MODE_SYNCHRONIZED);
    node->chromosome_order = chromosome_order;
    node->num_chromosomes = num_chromosomes;
    return node;
}

void merged_node_free(merged_node* node) {
    assert(node);
    region_free(node->region);
    array_list_free(node->links, NULL);
}

vcf_record_file_link *vcf_record_file_link_new(vcf_record_t *record, vcf_file_t *file) {
    vcf_record_file_link *link = malloc(sizeof(vcf_record_file_link));
    link->record = record;
    link->file = file;
    return link;
}

void vcf_record_file_link_free(vcf_record_file_link *link) {
    assert(link);
    vcf_record_free_deep(link->record);
    free(link);
}

// void vcf_record_file_link_free(vcf_record_file_link *link, int free_record, int free_file) {
//     assert(link);
//     if (free_record) {
//         vcf_record_free_deep(link->record);
//     }
//     if (free_file) {
//         vcf_close(link->file);
//     }
//     free(link);
// }

int merged_node_add_link(vcf_record_t* record, vcf_file_t* file, merged_node* node) {
    assert(record);
    assert(file);
    assert(node);
    vcf_record_file_link *link = malloc(sizeof(vcf_record_file_link));
    link->record = record;
    link->file = file;
    return array_list_insert(link, node->links);
}

vcf_record_file_link *merged_node_remove_link(vcf_record_file_link *link, merged_node* node) {
    assert(link);
    assert(node);
    return array_list_remove(link, node->links);
}


void set_cmp_criteria(char** chromosome_order, int num_chromosomes, merged_node *node) {
    assert(chromosome_order);
    assert(node);
    node->chromosome_order = chromosome_order;
    node->num_chromosomes = num_chromosomes;
}


int merged_node_cmp(merged_node *node1, merged_node *node2) {
    return compare_regions(node1->region, node2->region, node1->chromosome_order, node1->num_chromosomes);
}
