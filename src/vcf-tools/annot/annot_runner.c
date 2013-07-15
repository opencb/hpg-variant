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

#include "annot.h"

int vcf_annot_chr_cmp(const vcf_annot_chr_t *v1, const vcf_annot_chr_t *v2);
void vcf_annot_sample_free(vcf_annot_sample_t *sample);
void vcf_annot_chr_free(vcf_annot_chr_t *chr);
void vcf_annot_pos_free(vcf_annot_pos_t *pos);
int vcf_annot_process_chunk(vcf_record_t **variants, int num_variants, array_list_t *sample_list, khash_t(bams)* sample_bams, vcf_file_t *vcf_file);

int run_annot(shared_options_data_t *shared_options_data, annot_options_data_t *options_data) {
    //    file_stats_t *file_stats = file_stats_new();
    //    sample_stats_t **sample_stats;
    //    
    // List that stores the batches of records filtered by each thread
    list_t *output_list[shared_options_data->num_threads];
    // List that stores which thread filtered the next batch to save
    list_t *next_token_list = malloc(sizeof(list_t));

    array_list_t *sample_list = array_list_new(100, 1.25f, COLLECTION_MODE_SYNCHRONIZED);

    int ret_code;
    double start, stop, total;

    vcf_file_t *vcf_file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
    if (!vcf_file) {
        LOG_FATAL("VCF file does not exist!\n");
    }

    for (int i = 0; i < shared_options_data->num_threads; i++) {
        output_list[i] = (list_t*) malloc(sizeof(list_t));
        list_init("input", 1, shared_options_data->num_threads * shared_options_data->batch_lines, output_list[i]);
    }
    list_init("next_token", shared_options_data->num_threads, INT_MAX, next_token_list);

    LOG_INFO("Annotating VCF file...");

    khash_t(bams) *sample_bams = kh_init(bams);


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
            LOG_INFO_F("Thread %d processes data\n", omp_get_thread_num());

            individual_t **individuals = NULL;
            khash_t(ids) *sample_ids = NULL;

            start = omp_get_wtime();

            int i = 0;

            vcf_batch_t *batch = NULL;
            while ((batch = fetch_vcf_batch(vcf_file)) != NULL) {
                if (i % 50 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->records->size, batch->records->capacity);
                }

                // Divide the list of passed records in ranges of size defined in config file

                int num_chunks;
                int *chunk_sizes = NULL;
                array_list_t *input_records = batch->records;
                int *chunk_starts = create_chunks(input_records->size, shared_options_data->entries_per_thread, &num_chunks, &chunk_sizes);
#pragma omp parallel for num_threads(shared_options_data->num_threads) 
                for (int j = 0; j < num_chunks; j++) {
                    vcf_annot_process_chunk((vcf_record_t**)(input_records->items + chunk_starts[j]), chunk_sizes[j], sample_list, sample_bams, vcf_file);
                }

                free(chunk_starts);
                free(chunk_sizes);
                vcf_batch_free(batch);
                i++;
            }

            stop = omp_get_wtime();
            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Decrease list writers count
            for (i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(next_token_list);
                list_decr_writers(output_list[i]);
            }

            if (sample_ids) { kh_destroy(ids, sample_ids); }
            if (individuals) { free(individuals); }
        }
    }


    printf ( "Tam: %d\n", kh_size(sample_bams) );
    khiter_t iter;
    char* aux;
	for (iter = kh_begin(sample_bams); iter != kh_end(sample_bams); ++iter){
		if (kh_exist(sample_bams, iter)){
            aux = kh_value(sample_bams, iter);
            printf ( "%s - %s\n",kh_key(sample_bams, iter),  aux );
        }
    }

    

    vcf_annot_sample_t *annot_sample;
    vcf_annot_chr_t *annot_chr;
    vcf_annot_pos_t *annot_pos;
    int total_missings = 0;

    for (int i = 0; i < array_list_size(sample_list); i++) {
        annot_sample = (vcf_annot_sample_t*) sample_list->items[i];
        printf( "%s (%d)\n", annot_sample->name, array_list_size(annot_sample->chromosomes) );
        total_missings = 0;
        for(int j = 0; j < array_list_size(annot_sample->chromosomes); j++){
            annot_chr = (vcf_annot_chr_t*) array_list_get(j, annot_sample->chromosomes);
            printf ( "\t%s (%d)\n", annot_chr->name, array_list_size(annot_chr->positions));
            total_missings += array_list_size(annot_chr->positions);
        }
        printf ( "\tTotal: %d\n", total_missings );
    }


    /// FREE

    free(next_token_list);
    array_list_free(sample_list, vcf_annot_sample_free);
    for (int i = 0; i < shared_options_data->num_threads; i++) {
        free(output_list[i]);
    }

    vcf_close(vcf_file);
    return 0;
}

int vcf_annot_chr_cmp(const vcf_annot_chr_t *v1, const vcf_annot_chr_t *v2){
    return strcmp(v1->name, v2->name);
}

void vcf_annot_sample_free(vcf_annot_sample_t *sample){
    if(sample != NULL)
    {
        array_list_free(sample->chromosomes, vcf_annot_chr_free);
        free(sample->name);
        free(sample);
    }
}

void vcf_annot_chr_free(vcf_annot_chr_t *chr){
    if(chr != NULL)
    {
        array_list_free(chr->positions, vcf_annot_pos_free);
        free(chr->name);
        free(chr);
    }
}

void vcf_annot_pos_free(vcf_annot_pos_t *pos){
    if(pos != NULL)
    {
        free(pos);
    }
}

int vcf_annot_process_chunk(vcf_record_t **variants, int num_variants, array_list_t *sample_list, khash_t(bams)* sample_bams, vcf_file_t *vcf_file){

    int gt_pos, i;
    int allele1, allele2, alleles_code;
    vcf_record_t *record;
    vcf_annot_sample_t *annot_sample;
    vcf_annot_chr_t * annot_chr;
    vcf_annot_pos_t * annot_pos;
    char * sample_name;
    char * chr;
    char * copy_buf;
    unsigned int pos;
    size_t array_pos;
    khiter_t iter;
    int ret;

    for(int j = 0; j < num_variants; j++){
        record = variants[j]; 
        pos = record->position;
        chr = strndup(record->chromosome, record->chromosome_len);

#pragma omp critical
        {
            if(array_list_size(sample_list) == 0){
                for (int n = 0; n < array_list_size(vcf_file->samples_names); n++) {
                    annot_sample = (vcf_annot_sample_t*) malloc(sizeof(vcf_annot_sample_t));
                    sample_name = (char*) array_list_get(n, vcf_file->samples_names);
                    annot_sample->name = strndup(sample_name, strlen(sample_name));
                    annot_sample->chromosomes = array_list_new(24, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
                    annot_sample->chromosomes->compare_fn = &vcf_annot_chr_cmp;
                    iter = kh_get(ids, sample_bams, annot_sample->name);
                    if (iter != kh_end(sample_bams)) {
                        LOG_FATAL_F("Sample %s appears more than once. File can not be analyzed.\n", annot_sample->name);
                    } else {
                        iter = kh_put(bams, sample_bams, annot_sample->name, &ret);
                        if (ret) {
                            copy_buf = (char*) calloc(strlen(annot_sample->name) + 8 + 1, sizeof(char));
                            strcpy(copy_buf, annot_sample->name);
                            strcat(copy_buf, ".bam.bai");
                            kh_value(sample_bams, iter) = copy_buf; 
                        }
                    }


                    array_list_insert(annot_sample, sample_list);
                }
            }
        }

        for (int n = 0; n < array_list_size(sample_list); n++){
            copy_buf = strndup(record->format, record->format_len);
            gt_pos = get_field_position_in_format("GT", copy_buf); 
            if (copy_buf) {
                free(copy_buf);
                copy_buf = NULL;
            }
            if (gt_pos < 0) { 
                continue; 
            }   // This variant has no GT field

            copy_buf = strdup((char*) array_list_get(n, record->samples)); 
            alleles_code = get_alleles(copy_buf, gt_pos, &allele1, &allele2);
            if (copy_buf) {
                free(copy_buf);
                copy_buf = NULL;
            }

            if (alleles_code == 3) {

                annot_sample = array_list_get(n, sample_list);
#pragma omp critical
                {
                    annot_chr = (vcf_annot_chr_t*) malloc(sizeof(vcf_annot_chr_t));
                    annot_chr->name = chr;
                    array_pos = array_list_index_of(annot_chr, annot_sample->chromosomes);
                    if(array_pos != ULONG_MAX){
                        // liberar annot_chr previo
                        annot_chr = array_list_get(array_pos, annot_sample->chromosomes);
                    }
                    else
                    {
                        annot_chr->positions = array_list_new(100000, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
                        array_list_insert(annot_chr, annot_sample->chromosomes);
                    }
                }

                annot_pos = (vcf_annot_pos_t*) malloc(sizeof(vcf_annot_pos_t));
                annot_pos->pos = pos;
                array_list_insert(annot_pos, annot_chr->positions);
            }
        }
    }

    return 0;
}
