/*
 * Copyright (c) 2012-2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * hpg-variant is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-variant. If not, see <http://www.gnu.org/licenses/>.
 */

#include "stats.h"


static int write_output_variant_alleles_stats(variant_stats_t *var_stats, FILE *stats_fd);    
static int write_output_variant_genotypes_stats(variant_stats_t *var_stats, FILE *stats_fd);
static inline int write_output_variant_missing_data(variant_stats_t *var_stats, FILE *stats_fd);
static inline int write_output_variant_inheritance_data(variant_stats_t *var_stats, FILE *stats_fd);


int run_stats(shared_options_data_t *shared_options_data, stats_options_data_t *options_data) {
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, MIN(10, shared_options_data->max_batches) * shared_options_data->batch_lines, output_list);
    file_stats_t *file_stats = file_stats_new();
    sample_stats_t **sample_stats;

    int ret_code;
    double start, stop, total;
    
    vcf_file_t *vcf_file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
    if (!vcf_file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ped_file_t *ped_file = NULL;
    if (shared_options_data->ped_filename) {
        ped_file = ped_open(shared_options_data->ped_filename);
        if (!ped_file) {
            LOG_FATAL("PED file does not exist!\n");
        }
        LOG_INFO("About to read PED file...\n");
        // Read PED file before doing any processing
        ret_code = ped_read(ped_file);
        if (ret_code != 0) {
            LOG_FATAL_F("Can't read PED file: %s\n", ped_file->filename);
        }
    }
    
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }
    
    LOG_INFO("About to retrieve statistics from VCF file...\n");

#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            if (shared_options_data->batch_bytes > 0) {
                ret_code = vcf_parse_batches_in_bytes(shared_options_data->batch_bytes, vcf_file);
            } else if (shared_options_data->batch_lines > 0) {
                ret_code = vcf_parse_batches(shared_options_data->batch_lines, vcf_file);
            }

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            notify_end_parsing(vcf_file);
        }
        
#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            individual_t **individuals;
            khash_t(ids) *sample_ids = NULL;
            
            start = omp_get_wtime();
            
            int i = 0;
            vcf_batch_t *batch = NULL;
            while ((batch = fetch_vcf_batch(vcf_file)) != NULL) {
                if (i == 0) {
                    sample_stats = malloc (get_num_vcf_samples(vcf_file) * sizeof(sample_stats_t*));
                    for (int j = 0; j < get_num_vcf_samples(vcf_file); j++) {
                        sample_stats[j] = sample_stats_new(array_list_get(j, vcf_file->samples_names));
                    }
                    
                    if (ped_file) {
                        // Create map to associate the position of individuals in the list of samples defined in the VCF file
                        sample_ids = associate_samples_and_positions(vcf_file);
                        // Sort individuals in PED as defined in the VCF file
                        individuals = sort_individuals(vcf_file, ped_file);
                    }
                }
                
                if (i % 50 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->records->size, batch->records->capacity);
                }

                // Divide the list of passed records in ranges of size defined in config file
                int num_chunks;
                int *chunk_sizes;
                array_list_t *input_records = batch->records;
                int *chunk_starts = create_chunks(input_records->size, shared_options_data->entries_per_thread, &num_chunks, &chunk_sizes);
                
                assert(sample_ids);
                
                // OpenMP: Launch a thread for each range
                #pragma omp parallel for num_threads(shared_options_data->num_threads)
                for (int j = 0; j < num_chunks; j++) {
                    LOG_DEBUG_F("[%d] Stats invocation\n", omp_get_thread_num());
                    // Invoke variant stats and/or sample stats when applies
                    if (options_data->variant_stats) {
                        ret_code = get_variants_stats((vcf_record_t**) (input_records->items + chunk_starts[j]), 
                                                      chunk_sizes[j], individuals, sample_ids, output_list, file_stats);
                    }
                    if (options_data->sample_stats) {
                        ret_code |= get_sample_stats((vcf_record_t**) (input_records->items + chunk_starts[j]), 
                                                      chunk_sizes[j], individuals, sample_ids, sample_stats, file_stats);
                    }
                }
                
                free(chunk_starts);
                vcf_batch_free(batch);
                
                i++;
            }
            
            // Write sample statistics
            int dirname_len = strlen(shared_options_data->output_directory);
            char *stats_filename;
            FILE *stats_fd;
            
            if (options_data->sample_stats) {
                if (shared_options_data->output_filename == NULL || strlen(shared_options_data->output_filename) == 0) {
                    stats_filename = (char*) calloc ((dirname_len + strlen("stats-samples") + 2), sizeof(char));
                    sprintf(stats_filename, "%s/stats-samples", shared_options_data->output_directory);
                } else {
                    stats_filename = (char*) calloc ((dirname_len + strlen(shared_options_data->output_filename) + 10), sizeof(char));
                    sprintf(stats_filename, "%s/%s-samples", shared_options_data->output_directory, shared_options_data->output_filename);
                }
                
                stats_fd = fopen(stats_filename, "w");
                free(stats_filename);
                fprintf(stats_fd, "#SAMPLE\t\tMISS GT\t\tMENDEL ERR\n");
                
                sample_stats_t *sam_stats;
                for (int i = 0; i < vcf_file->samples_names->size; i++) {
                    sam_stats = sample_stats[i];
                    fprintf(stats_fd, "%s\t\t%zu\t\t%zu\n", sam_stats->name, sam_stats->missing_genotypes, sam_stats->mendelian_errors);
                    sample_stats_free(sam_stats);
                }
                
                // Close sample stats file
                if (stats_fd != NULL) { fclose(stats_fd); }
            }
            
            stop = omp_get_wtime();
            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
            
            // Decrease list writers count
            for (i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
            
            if (sample_ids) { kh_destroy(ids, sample_ids); }
            if (individuals) { free(individuals); }
        }
        
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d writes the output\n", omp_get_thread_num());
            
            char *stats_filename, *summary_filename;
            FILE *stats_fd, *summary_fd;
    
            // Create file streams (summary and results)
            int dirname_len = strlen(shared_options_data->output_directory);
            summary_filename = (char*) calloc ((dirname_len + strlen("summary-stats") + 2), sizeof(char));
            sprintf(summary_filename, "%s/summary-stats", shared_options_data->output_directory);
            
            // Write variant statistics
            if (options_data->variant_stats) {
                if (shared_options_data->output_filename == NULL || strlen(shared_options_data->output_filename) == 0) {
                    stats_filename = (char*) calloc ((dirname_len + strlen("stats-variants") + 2), sizeof(char));
                    sprintf(stats_filename, "%s/stats-variants", shared_options_data->output_directory);
                } else {
                    stats_filename = (char*) calloc ((dirname_len + strlen(shared_options_data->output_filename) + 11), sizeof(char));
                    sprintf(stats_filename, "%s/%s-variants", shared_options_data->output_directory, shared_options_data->output_filename);
                }
                
                LOG_DEBUG_F("stats filename = %s\nsummary filename = %s\n", stats_filename, summary_filename);
                stats_fd = fopen(stats_filename, "w");
                summary_fd = fopen(summary_filename, "w");
                free(stats_filename);
                free(summary_filename);
                LOG_DEBUG("File streams created\n");
                
                fprintf(stats_fd, 
                        "#CHROM\tPOS\tINDEL?\tList of [ALLELE  COUNT  FREQ]\t\t\tList of [GT  COUNT  FREQ]\t\t\t\t\t\tMISS_AL\tMISS_GT\tMEND_ER\t%% AFF | UNAFF dominant\t%% AFF | UNAFF recessive\n");
                
                // For each variant, generate a new line with the format (block of blanks = tab):
                // chromosome   position   [<allele>   <count>   <freq>]+   [<genotype>   <count>   <freq>]+   miss_all   miss_gt
                list_item_t* item = NULL;
                variant_stats_t *var_stats;
                while ((item = list_remove_item(output_list)) != NULL) {
                    var_stats = item->data_p;
                    
                    fprintf(stats_fd, "%s\t%ld\t",
                            var_stats->chromosome, 
                            var_stats->position);
                    
                    write_output_variant_alleles_stats(var_stats, stats_fd);
                    write_output_variant_genotypes_stats(var_stats, stats_fd);
                    write_output_variant_missing_data(var_stats, stats_fd);
                    write_output_variant_inheritance_data(var_stats, stats_fd);
                    
                    // Free resources
                    variant_stats_free(var_stats);
                    list_item_free(item);
                }
                
                // Close variant stats file
                if (stats_fd != NULL) { fclose(stats_fd); }
            
                // Write whole file stats (data only got when launching variant stats)
                fprintf(summary_fd, 
                        "Number of variants = %d\nNumber of samples = %d\nNumber of biallelic variants = %d\nNumber of multiallelic variants = %d\n\n",
                        file_stats->variants_count, file_stats->samples_count, file_stats->biallelics_count, file_stats->multiallelics_count
                    );
                
                fprintf(summary_fd, 
                        "Number of SNP = %d\nNumber of indels = %d\n\n",
                        file_stats->snps_count, file_stats->indels_count
                    );
                
                fprintf(summary_fd, 
                        "Number of transitions = %d\nNumber of transversions = %d\nTi/TV ratio = %.4f\n\nPercentage of PASS = %.2f%%\nAverage quality = %.2f\n",
                        file_stats->transitions_count, file_stats->transversions_count,
                        (float) file_stats->transitions_count / file_stats->transversions_count,
                        ((float) file_stats->pass_count / file_stats->variants_count) * 100.0,
                        file_stats->accum_quality / file_stats->variants_count
                    );
            }
        }
        
    }
    
    vcf_close(vcf_file);
    if (ped_file) { ped_close(ped_file, 1); }
    free(sample_stats);
    free(file_stats);
//     free(read_list);
    free(output_list);
    
    return 0;
}

static int write_output_variant_alleles_stats(variant_stats_t *var_stats, FILE *stats_fd) {
    int written = 0;
    
    // Is indel?
    written += (var_stats->is_indel) ? fprintf(stats_fd, "Y\t") : fprintf(stats_fd, "N\t");
    
    // Reference allele
    written += fprintf(stats_fd, "%s\t%d\t%.4f\t",
                       var_stats->ref_allele,
                       var_stats->alleles_count[0],
                       var_stats->alleles_freq[0]);

    // Alternate alleles
    for (int i = 1; i < var_stats->num_alleles; i++) {
        written += fprintf(stats_fd, "%s\t%d\t%.4f\t",
                           var_stats->alternates[i-1],
                           var_stats->alleles_count[i],
                           var_stats->alleles_freq[i]);
    }
    
    return written;
} 
    
static int write_output_variant_genotypes_stats(variant_stats_t *var_stats, FILE *stats_fd) {
    int written = 0;
    int gt_count = 0;
    float gt_freq = 0;
    
    for (int i = 0; i < var_stats->num_alleles; i++) {
        for (int j = i; j < var_stats->num_alleles; j++) {
            int idx1 = i * var_stats->num_alleles + j;
            if (i == j) {
                gt_count = var_stats->genotypes_count[idx1];
                gt_freq = var_stats->genotypes_freq[idx1];
            } else {
                int idx2 = j * var_stats->num_alleles + i;
                gt_count = var_stats->genotypes_count[idx1] + var_stats->genotypes_count[idx2];
                gt_freq = var_stats->genotypes_freq[idx1] + var_stats->genotypes_freq[idx2];
            }

            written += fprintf(stats_fd, "%s|%s\t%d\t%.4f\t",
                               i == 0 ? var_stats->ref_allele : var_stats->alternates[i-1],
                               j == 0 ? var_stats->ref_allele : var_stats->alternates[j-1],
                               gt_count, gt_freq);
        }
    }
    
    return written;
}

static inline int write_output_variant_missing_data(variant_stats_t *var_stats, FILE *stats_fd) {
    return fprintf(stats_fd, "%d\t%d\t",
                   var_stats->missing_alleles,
                   var_stats->missing_genotypes);
}

static inline int write_output_variant_inheritance_data(variant_stats_t *var_stats, FILE *stats_fd) {
    return fprintf(stats_fd, "%d\t%.2f | %.2f\t%.2f | %.2f\n",
                   var_stats->missing_alleles,
                   var_stats->cases_percent_dominant,
                   var_stats->controls_percent_dominant,
                   var_stats->cases_percent_recessive,
                   var_stats->controls_percent_recessive);
}
