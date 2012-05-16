#include "stats.h"

variant_stats_t* new_variant_stats(char *chromosome, unsigned long position, char *ref_allele) {
    variant_stats_t *stats = (variant_stats_t*) malloc (sizeof(variant_stats_t));
    
    stats->chromosome = chromosome;
    stats->position = position;
    stats->ref_allele = ref_allele;
    stats->alternates = NULL;
    stats->num_alleles = 1;
    stats->alleles_count = NULL;
    stats->genotypes_count = NULL;
    
    stats->missing_alleles = 0;
    stats->missing_genotypes = 0;
    
    return stats;
}

void free_variant_stats(variant_stats_t* stats) {
    if (stats->ref_allele) { free(stats->ref_allele); }
    if (stats->alternates) {
        for (int i = 0; i < stats->num_alleles; i++) {
            free(stats->alternates[i]);
        }
        free(stats->alternates);
    }
    if (stats->alleles_count) { free(stats->alleles_count); }
    if (stats->genotypes_count) { free(stats->genotypes_count); }
    free(stats);
}

int get_variants_stats(list_item_t* variants, int num_variants, list_t* output_list) {
    char *copy_buf, *copy_buf2, *token, *sample;
    char *save_strtok;
    
    int num_alternates, gt_pos, cur_pos;
    int allele1, allele2, alleles_code;
    
    vcf_record_t *record;
    variant_stats_t *stats;
    list_item_t *cur_variant = variants;
    
    for (int i = 0; i < num_variants && cur_variant != NULL; i++) {
        record = (vcf_record_t*) cur_variant->data_p;
        copy_buf = (char*) calloc (strlen(record->chromosome)+1, sizeof(char));
        strncat(copy_buf, record->chromosome, strlen(record->chromosome));
        copy_buf2 = (char*) calloc (strlen(record->reference)+1, sizeof(char));
        strncat(copy_buf2, record->reference, strlen(record->reference));
        stats = new_variant_stats(copy_buf, record->position, copy_buf2);
        
        // Create list of alternates
        copy_buf = (char*) calloc (strlen(record->alternate)+1, sizeof(char));
        strcat(copy_buf, record->alternate);
        stats->alternates = split(copy_buf, ",", &num_alternates);
        if (!strcmp(stats->alternates[0], ".")) {
            stats->num_alleles = 1;
        } else {
            stats->num_alleles = num_alternates + 1;
        }
        LOG_INFO_F("num alternates = %d\tnum_alleles = %d\n", num_alternates, stats->num_alleles);
        
        // Create lists of allele and genotypes counters
        stats->alleles_count = (int*) calloc (stats->num_alleles, sizeof(int));
        stats->genotypes_count = (int*) calloc (stats->num_alleles * stats->num_alleles, sizeof(int));
        
        // Get position where GT is in sample
        copy_buf = (char*) calloc (strlen(record->format)+1, sizeof(char));
        strcat(copy_buf, record->format);
        gt_pos = get_genotype_position_in_format(copy_buf);
        LOG_INFO_F("Genotype position = %d\n", gt_pos);
        if (gt_pos < 0) { continue; }   // This variant has no GT field
        
        // Traverse samples and find the present and missing alleles
        for(list_item_t *cur_sample = record->samples->first_p; cur_sample != NULL; cur_sample = cur_sample->next_p) {
            sample = (char*) cur_sample->data_p;
            
            // Get to GT position
            copy_buf = (char*) calloc (strlen(sample)+1, sizeof(char));
            strcat(copy_buf, sample);
            alleles_code = get_alleles(copy_buf, gt_pos, &allele1, &allele2);
            LOG_DEBUG_F("sample = %s, alleles = %d/%d\n", sample, allele1, allele2);
            
            if (allele1 < 0 || allele2 < 0) {
                // Missing genotype (one or both alleles missing)
                stats->missing_genotypes++;
                if (allele1 < 0) { 
                    stats->missing_alleles++; 
                } else {
                    stats->alleles_count[allele1]++;
                }
                    
                if (allele2 < 0) { 
                    stats->missing_alleles++;
                } else {
                    stats->alleles_count[allele2]++;
                }
            } else {
                // Both alleles set
                cur_pos = allele1 * (stats->num_alleles) + allele2;
                stats->alleles_count[allele1]++;
                stats->alleles_count[allele2]++;
                stats->genotypes_count[cur_pos]++;
            }
        }
        
        // Insert results in output list
        list_item_t *variant_result = list_item_new(i, 0, stats);
        list_insert_item(variant_result, output_list);
    }
    
    return 0;
}

