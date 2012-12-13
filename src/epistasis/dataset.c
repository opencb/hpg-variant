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


void epistasis_dataset_process_records(vcf_record_t **variants, size_t num_variants, int num_samples, int *phenotypes, epistasis_dataset *dataset) {
    for (int i = 0; i < num_variants; i++) {
        vcf_record_t *record = variants[i];
        int gt_position = get_field_position_in_format("GT", strndup(record->format, record->format_len));
        uint8_t *genotypes = malloc (num_samples * sizeof(uint8_t));
        
        // For each sample, get genotype and increment counters in the dataset
        for (int k = 0; k < num_samples; k++) {
            char *sample = strdup(array_list_get(k, record->samples));
            int allele1, allele2, gt_dataset_index;
            if (get_alleles(sample, gt_position, &allele1, &allele2)) {
                LOG_INFO_F("Missing sample #%d in variant %.*s:%ld will be ignored\n",
                            k, record->chromosome_len, record->chromosome, record->position);
                genotypes[k] = 255;
            } else {
                LOG_DEBUG_F("alleles = %d/%d, phenotype = %d\n", allele1, allele2, phenotypes[k]);
                if (!allele1 && !allele2) { // Homozygous in first allele
                    genotypes[k] = 0;
                } else if (allele1 != allele2) { // Heterozygous
                    genotypes[k] = 1;
                } else if (allele1 && allele1 == allele2) { // Homozygous in second allele
                    genotypes[k] = 2;
                }
            }
            free(sample);
        }
        
        array_list_insert(genotypes, dataset->genotype_counts);
    }
}

void epistasis_dataset_add_entry(uint16_t *genotype_count, epistasis_dataset *dataset) {
    array_list_insert(genotype_count, dataset->genotype_counts);    
}


size_t epistasis_dataset_get_num_variants(epistasis_dataset *dataset) {
    return dataset->genotype_counts->size;
}

uint8_t *epistasis_dataset_get_variant_counts(size_t index, epistasis_dataset *dataset) {
    return array_list_get(index, dataset->genotype_counts);
}