#include "dataset.h"

epistasis_dataset* epistasis_dataset_new() {
    epistasis_dataset *dataset = malloc (sizeof(epistasis_dataset));
    dataset->num_affected = 0;
    dataset->num_unaffected = 0;
    dataset->genotype_counts = array_list_new(1000, 1.5, COLLECTION_MODE_SYNCHRONIZED);
    return dataset;
}

void epistasis_dataset_free(epistasis_dataset *dataset) {
    array_list_free(dataset->genotype_counts, free);
    free(dataset);
}


void epistasis_dataset_process_records(vcf_record_t **variants, int num_variants, int num_samples, int *phenotypes, epistasis_dataset *dataset) {
    for (int i = 0; i < num_variants; i++) {
        vcf_record_t *record = variants[i];
        int gt_position = get_field_position_in_format("GT", strndup(record->format, record->format_len));
        uint16_t *genotypes = calloc (6, sizeof(uint16_t));
        
        // For each sample, get genotype and increment counters in the dataset
        for (int k = 0; k < num_samples; k++) {
            char *sample = strdup(array_list_get(k, record->samples));
            int allele1, allele2, gt_dataset_index;
            if (get_alleles(sample, gt_position, &allele1, &allele2)) {
                LOG_WARN_F("Sample #k in variant %.*s:%ld is missing\n", 
                            k, record->chromosome_len, record->chromosome, record->position);
            } else {
                LOG_DEBUG_F("alleles = %d/%d, phenotype = %d\n", allele1, allele2, phenotypes[k]);
                // Increment genotype count depending on phenotype
                if (!allele1 && allele1 == allele2) { // Homozygous in first allele
                    gt_dataset_index = 0;
                } else if ((!allele1 && allele2) || (allele1 && !allele2)) { // Heterozygous
                    gt_dataset_index = 2;
                } else if (allele1 && allele1 == allele2) { // Homozygous in second allele
                    gt_dataset_index = 4;
                }
                
                if (phenotypes[k]) { // Case
                    (genotypes[gt_dataset_index])++;
                } else { // Control
                    (genotypes[gt_dataset_index + 1])++;
                }
                
//                 printf("genotypes %d { ", i);
//                 for (int i = 0; i < 6; i++) {
//                     printf("%d ", genotypes[i]);
//                 }
//                 printf(" }\n");
    
            }
            free(sample);
        }
        
        array_list_insert(genotypes, dataset->genotype_counts);
    }
}

void epistasis_dataset_add_entry(uint16_t *genotype_count, epistasis_dataset *dataset) {
    array_list_insert(genotype_count, dataset->genotype_counts);    
}


uint16_t epistasis_dataset_get_num_variants(epistasis_dataset *dataset) {
    return dataset->genotype_counts->size;
}

uint16_t *epistasis_dataset_get_variant_counts(size_t index, epistasis_dataset *dataset) {
    return array_list_get(index, dataset->genotype_counts);
}
