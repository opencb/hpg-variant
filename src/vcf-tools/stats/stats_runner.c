#include "stats.h"

int run_stats(shared_options_data_t *shared_options_data, stats_options_data_t *options_data) {
//     list_t *read_list = (list_t*) malloc(sizeof(list_t));
//     list_init("batches", 1, shared_options_data->max_batches, read_list);
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, MIN(10, shared_options_data->max_batches) * shared_options_data->batch_lines, output_list);
//     file_stats_t *file_stats = new_file_stats();
    file_stats_t *file_stats = file_stats_new();

    int ret_code;
    double start, stop, total;
    vcf_file_t *file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
    
    if (!file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }
    
#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            if (shared_options_data->batch_bytes > 0) {
                ret_code = vcf_parse_batches_in_bytes(shared_options_data->batch_bytes, file, 1);
            } else if (shared_options_data->batch_lines > 0) {
                ret_code = vcf_parse_batches(shared_options_data->batch_lines, file, 1);
            }

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            notify_end_reading(file);
//             list_decr_writers(read_list);
        }
        
#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            omp_set_num_threads(shared_options_data->num_threads);
            
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            start = omp_get_wtime();
            
            int i = 0;
            vcf_batch_t *batch = NULL;
            while ((batch = fetch_vcf_batch(file)) != NULL) {
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
                
                // OpenMP: Launch a thread for each range
                #pragma omp parallel for
                for (int j = 0; j < num_chunks; j++) {
                    LOG_DEBUG_F("[%d] Stats invocation\n", omp_get_thread_num());
                    ret_code = get_variants_stats((vcf_record_t**) (input_records->items + chunk_starts[j]), 
                                                  chunk_sizes[j], 
                                                  output_list,
                                                  file_stats);
                }
                
                free(chunk_starts);
                vcf_batch_free(batch);
                
                i++;
            }
            
            stop = omp_get_wtime();
            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
            
            // Decrease list writers count
            for (i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }
        
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d writes the output\n", omp_get_thread_num());
            
            char *stats_filename, *summary_filename;
            FILE *stats_fd, *summary_fd;
    
            // Create file streams (summary and results)
            int dirname_len = strlen(shared_options_data->output_directory);
//             int filename_len = strlen("summary-stats");
            
            summary_filename = (char*) calloc ((dirname_len + strlen("summary-stats") + 2), sizeof(char));
            sprintf(summary_filename, "%s/summary-stats", shared_options_data->output_directory);
            
            if (shared_options_data->output_filename == NULL || strlen(shared_options_data->output_filename) == 0) {
//                 filename_len = strlen("stats-tool-output");
            
                stats_filename = (char*) calloc ((dirname_len + strlen("stats-tool-output") + 2), sizeof(char));
                sprintf(stats_filename, "%s/stats-tool-output", shared_options_data->output_directory);
            } else {
//                 filename_len = strlen(shared_options_data->output_filename);
            
                stats_filename = (char*) calloc ((dirname_len + strlen(shared_options_data->output_filename) + 2), sizeof(char));
                sprintf(stats_filename, "%s/%s", shared_options_data->output_directory, shared_options_data->output_filename);
            }
            
            LOG_DEBUG_F("stats filename = %s\nsummary filename = %s\n", stats_filename, summary_filename);
            stats_fd = fopen(stats_filename, "w");
            summary_fd = fopen(summary_filename, "w");
            free(stats_filename);
            free(summary_filename);
            LOG_DEBUG("File streams created\n");
            
            fprintf(stats_fd, "#CHROM\tPOS\tList of [ALLELE  COUNT  FREQ]\tList of [GT  COUNT  FREQ]\tMISS_ALLELES\tMISS_GT\n");
            
            // For each variant, generate a new line with the format (block of blanks = tab):
            // chromosome   position   [<allele>   <count>   <freq>]+   [<genotype>   <count>   <freq>]+   miss_all   miss_gt
            list_item_t* item = NULL;
            variant_stats_t *stats;
            int num_alleles, count_alleles_total, genotypes_count_total;
            FILE *fd = NULL;
       
            while ((item = list_remove_item(output_list)) != NULL) {
                stats = item->data_p;
                num_alleles = stats->num_alleles;
                
                fprintf(stats_fd, "%s\t%ld\t",
                        stats->chromosome, 
                        stats->position);
                
                // Reference allele
                fprintf(stats_fd, "%s\t%d\t%.4f\t",
                        stats->ref_allele,
                        stats->alleles_count[0],
                        stats->alleles_freq[0]
                       );
                
                // Alternate alleles
                for (int i = 1; i < num_alleles; i++) {
                    fprintf(stats_fd, "%s\t%d\t%.4f\t",
                            stats->alternates[i-1],
                            stats->alleles_count[i],
                            stats->alleles_freq[i]
                           );
                }
                
                // Genotypes
                int gt_count = 0;
                float gt_freq = 0;
                for (int i = 0; i < num_alleles; i++) {
                    for (int j = i; j < num_alleles; j++) {
                        int idx1 = i * num_alleles + j;
                        if (i == j) {
                            gt_count = stats->genotypes_count[idx1];
                            gt_freq = stats->genotypes_freq[idx1];
                        } else {
                            int idx2 = j * num_alleles + i;
                            gt_count = stats->genotypes_count[idx1] + stats->genotypes_count[idx2];
                            gt_freq = stats->genotypes_freq[idx1] + stats->genotypes_freq[idx2];
                        }
                        
                    fprintf(stats_fd, "%s|%s\t%d\t%.4f\t",
                            i == 0 ? stats->ref_allele : stats->alternates[i-1],
                            j == 0 ? stats->ref_allele : stats->alternates[j-1],
                            gt_count, 
                            gt_freq
                           );
                    }
                }
                
                // Missing data
                fprintf(stats_fd, "%d\t%d\n",
                        stats->missing_alleles,
                        stats->missing_genotypes
                       );
                
                // Free resources
                variant_stats_free(stats);
                list_item_free(item);
            }
            
            // Write whole file stats
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
                   
            
            // Close file
            if (stats_fd != NULL) { fclose(stats_fd); }
        }
        
    }
    
    vcf_close(file);
    free(file_stats);
//     free(read_list);
    free(output_list);
    
    return 0;
}

