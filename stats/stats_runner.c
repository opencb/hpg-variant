#include "stats.h"

int run_stats(global_options_data_t *global_options_data, stats_options_data_t *options_data) {
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("batches", 1, options_data->max_batches, read_list);
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", options_data->num_threads, MIN(10, options_data->max_batches) * options_data->batch_size, output_list);
    file_stats_t *file_stats = new_file_stats();

    int ret_code;
    double start, stop, total;
    vcf_file_t *file = vcf_open(global_options_data->vcf_filename);
    
    if (!file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ret_code = create_directory(global_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", global_options_data->output_directory);
    }
    
#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            ret_code = vcf_read_batches(read_list, options_data->batch_size, file, 1);

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            list_decr_writers(read_list);
        }
        
#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            omp_set_num_threads(options_data->num_threads);
            
            start = omp_get_wtime();
            
            int i = 0;
            list_item_t* item = NULL;
            while ((item = list_remove_item(read_list)) != NULL) {
                vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
                list_t *input_records = batch;

                if (i % 50 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                                i, omp_get_thread_num(),
                                batch->length, batch->max_length);
                }

                // Divide the list of passed records in ranges of size defined in config file
                int max_chunk_size = options_data->variants_per_thread;
                int num_chunks;
                list_item_t **chunk_starts = create_chunks(input_records, max_chunk_size, &num_chunks);
                
                // OpenMP: Launch a thread for each range
                #pragma omp parallel for
                for (int j = 0; j < num_chunks; j++) {
                    LOG_DEBUG_F("[%d] Stats invocation\n", omp_get_thread_num());
//                     printf("[%d] start in (%s, %ld)\n", j,
//                            ((vcf_record_t*) chunk_starts[j]->data_p)->chromosome, 
//                            ((vcf_record_t*) chunk_starts[j]->data_p)->position);
                    ret_code = get_variants_stats(chunk_starts[j], max_chunk_size, output_list, file_stats);
                }
                if (i % 25 == 0) { LOG_INFO_F("*** %dth stats invocation finished\n", i); }
                
                free(chunk_starts);
                vcf_batch_free(item->data_p);
                list_item_free(item);
                
                i++;
            }
            
            stop = omp_get_wtime();
            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
            
            // Decrease list writers count
            for (i = 0; i < options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }
        
#pragma omp section
        {
            char *stats_filename, *summary_filename;
            FILE *stats_fd, *summary_fd;
    
            // Create file streams (summary and results)
            int dirname_len = strlen(global_options_data->output_directory);
            int filename_len = strlen("summary-stats");
            
            summary_filename = (char*) calloc ((dirname_len + strlen("summary-stats") + 1), sizeof(char));
            strncat(summary_filename, global_options_data->output_directory, dirname_len);
            strncat(summary_filename, "summary-stats", strlen("summary-stats"));
            
            if (global_options_data->output_filename == NULL || strlen(global_options_data->output_filename) == 0) {
                filename_len = strlen("stats-tool-output");
            
                stats_filename = (char*) calloc ((dirname_len + filename_len + 1), sizeof(char));
                strncat(stats_filename, global_options_data->output_directory, dirname_len);
                strncat(stats_filename, "stats-tool-output", filename_len);
            } else {
                filename_len = strlen(global_options_data->output_filename);
            
                stats_filename = (char*) calloc ((dirname_len + filename_len + 1), sizeof(char));
                strncat(stats_filename, global_options_data->output_directory, dirname_len);
                strncat(stats_filename, global_options_data->output_filename, filename_len);
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
                
                // Generate global counters for alleles and genotypes (used for calculating frequencies)
                count_alleles_total = 0;
                genotypes_count_total = 0;
                for (int i = 0; i < num_alleles; i++) {
                    count_alleles_total += stats->alleles_count[i];
                }
                for (int i = 0; i < num_alleles * num_alleles; i++) {
                    genotypes_count_total += stats->genotypes_count[i];
                }
                
                // Chromosome and position
                fprintf(stats_fd, "%s\t%ld\t",
                        stats->chromosome, 
                        stats->position);
                
                // Reference allele
                fprintf(stats_fd, "%s\t%d\t%.4f\t",
                        stats->ref_allele,
                        stats->alleles_count[0],
                        (float) stats->alleles_count[0] / count_alleles_total
                       );
                
                // Alternate alleles
                for (int i = 1; i < num_alleles; i++) {
                    fprintf(stats_fd, "%s\t%d\t%.4f\t",
                            stats->alternates[i-1],
                            stats->alleles_count[i],
                            (float) stats->alleles_count[i] / count_alleles_total
                           );
                }
                
                // Genotypes
                int gt_count = 0;
                for (int i = 0; i < num_alleles; i++) {
                    for (int j = i; j < num_alleles; j++) {
                        int idx1 = i * num_alleles + j;
                        if (i == j) {
                            gt_count = stats->genotypes_count[idx1];
                        } else {
                            int idx2 = j * num_alleles + i;
                            gt_count = stats->genotypes_count[idx1] + stats->genotypes_count[idx2];
                        }
                        
                    fprintf(stats_fd, "%d|%d\t%d\t%.4f\t",
                            i, j, 
                            gt_count, 
                            (float) gt_count / genotypes_count_total
                           );
                    }
                }
                
                // Missing data
                fprintf(stats_fd, "%d\t%d\n",
                        stats->missing_alleles,
                        stats->missing_genotypes
                       );
                
                // Free resources
                free_variant_stats(stats);
                list_item_free(item);
            }
            
            // Write whole file stats
            fprintf(summary_fd, 
                    "Number of samples = %d\nNumber of SNP = %d\nNumber of indels = %d\nNumber of multiallelic variants = %d\n",
                    file_stats->samples_count, file_stats->snps_count, file_stats->indels_count, file_stats->multiallelics_count
                   );
                    
            fprintf(summary_fd, 
                    "Number of transitions = %d\nNumber of transversions = %d\nTi/TV ratio = %.4f\n",
                    file_stats->transitions_count, file_stats->transversions_count,
                    (float) file_stats->transitions_count / file_stats->transversions_count
                   );
                   
            
            // Close file
            if (stats_fd != NULL) { fclose(stats_fd); }
        }
        
    }
    
    vcf_close(file);
    free(file_stats);
    free(read_list);
    free(output_list);
    
    return 0;
}

