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

int vcf_annot_process_chunk(vcf_record_t **variants, int num_variants, array_list_t *sample_list, vcf_file_t *vcf_file);

static void vcf_annot_sort_sample(vcf_annot_sample_t* annot_sample);

int vcf_annot_check_bams(vcf_annot_sample_t* annot_sample, khash_t(bams)* sample_bams);

int vcf_annot_edit_chunk(vcf_record_t **variants, int num_variants, array_list_t *sample_list, khash_t(bams)* sample_bams, vcf_file_t *vcf_file);

static vcf_annot_pos_t * vcf_annot_get_pos(char *sample_name, char *chr, size_t pos, array_list_t *sample_list);

static int vcf_annot_pos_cmp (const void * a, const void * b);

static int binary_search_pos(array_list_t* array_list , unsigned int f,unsigned int l, unsigned int target);
void set_field_value_in_sample(char **sample, int position, char* value);

// annot.c
static int count_func(const bam1_t *b, void *data)
{
		(*((count_func_data_t*)data)->count)++;
	return 0;
}


int run_annot(shared_options_data_t *shared_options_data, annot_options_data_t *options_data) {

    int ret_code;
    double start, stop, total;
    vcf_annot_sample_t *annot_sample;
    vcf_annot_chr_t *annot_chr;
    vcf_annot_pos_t *annot_pos;
    vcf_annot_bam_t *annot_bam;
    char *sample_name;
    char * copy_buf;
    khiter_t iter;
    int ret;
    list_item_t *output_item;

    char *directory = options_data->bam_directory;

    vcf_file_t *vcf_file = vcf_open(shared_options_data->vcf_filename, shared_options_data->max_batches);
    if (!vcf_file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }

    LOG_INFO("Annotating VCF file...");
    
    array_list_t *sample_list = array_list_new(100, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
    
    list_t *output_list = malloc(sizeof(list_t));
    list_init("output_list", shared_options_data->num_threads, INT_MAX, output_list);

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

            start = omp_get_wtime();

            int i = 0;

            vcf_batch_t *batch = NULL;
            khash_t(bams) *sample_bams = kh_init(bams);
            while ((batch = fetch_vcf_batch(vcf_file)) != NULL) {
                if(i == 0) {
                    for (int n = 0; n < array_list_size(vcf_file->samples_names); n++) {
                        sample_name = (char*) array_list_get(n, vcf_file->samples_names);
                        iter = kh_get(ids, sample_bams, sample_name);
                        if (iter != kh_end(sample_bams)) {
                            LOG_FATAL_F("Sample %s appears more than once. File can not be analyzed.\n", annot_sample->name);
                        } else {
                            iter = kh_put(bams, sample_bams, sample_name, &ret);
                            if (ret) {
                                annot_bam = (vcf_annot_bam_t*) malloc(sizeof(vcf_annot_bam_t));
                                annot_bam->bam_filename = (char*) calloc(strlen(directory) + strlen(sample_name) + 4 + 1, sizeof(char));
                                strcpy(annot_bam->bam_filename, directory);
                                strcat(annot_bam->bam_filename, sample_name);
                                strcat(annot_bam->bam_filename, ".bam");
                                
                                if(!exists(annot_bam->bam_filename)){
                                    LOG_FATAL_F("File %s does not exist\n", annot_bam->bam_filename);
                                }

                                annot_bam->bai_filename = (char*) calloc(strlen(directory) + strlen(sample_name) + 8 + 1, sizeof(char));
                                strcpy(annot_bam->bai_filename, directory);
                                strcat(annot_bam->bai_filename, sample_name);
                                strcat(annot_bam->bai_filename, ".bam.bai");
                                
                                if(!exists(annot_bam->bai_filename)){
                                    LOG_FATAL_F("File %s does not exist\n", annot_bam->bai_filename);
                                }

                                kh_value(sample_bams, iter) = annot_bam; 
                            }
                        }
                    }
                    // Check BAM files
                }

                if(i % 50 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->records->size, batch->records->capacity);
                }
                
                for (int n = 0; n < array_list_size(vcf_file->samples_names); n++) {
                    annot_sample = (vcf_annot_sample_t*) malloc(sizeof(vcf_annot_sample_t));
                    sample_name = (char*) array_list_get(n, vcf_file->samples_names);
                    annot_sample->name = strndup(sample_name, strlen(sample_name));
                    annot_sample->chromosomes = array_list_new(24, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
                    annot_sample->chromosomes->compare_fn = &vcf_annot_chr_cmp;
                    array_list_insert(annot_sample, sample_list);
                }

                // Divide the list of passed records in ranges of size defined in config file

                array_list_t *input_records = batch->records;
                bam_file_t* bam_file;
                bam_header_t *bam_header;
                bam_index_t *idx = 0;
                char* aux;
                char *bam_filename;
                int *chunk_sizes = NULL;
                int result, tid, beg, end;
                int num_chunks;
                khiter_t iter;
                vcf_annot_chr_t *annot_chr;
                vcf_annot_pos_t *annot_pos;
                vcf_annot_sample_t *annot_sample;
                int *chunk_starts = create_chunks(input_records->size, shared_options_data->entries_per_thread, &num_chunks, &chunk_sizes);

#pragma omp parallel for num_threads(shared_options_data->num_threads) 
                for (int j = 0; j < num_chunks; j++) {
                    vcf_annot_process_chunk((vcf_record_t**)(input_records->items + chunk_starts[j]), chunk_sizes[j], sample_list, vcf_file);
                }

#pragma omp parallel for num_threads(shared_options_data->num_threads) 
                for (int j = 0; j < array_list_size(sample_list); j++) {
                    annot_sample = (vcf_annot_sample_t*) sample_list->items[j];
                    vcf_annot_sort_sample(annot_sample);
                    vcf_annot_check_bams(annot_sample, sample_bams);
                }

#pragma omp parallel for num_threads(shared_options_data->num_threads) 
                for (int j = 0; j < num_chunks; j++) {
                    vcf_annot_edit_chunk((vcf_record_t**)(input_records->items + chunk_starts[j]), chunk_sizes[j], sample_list, sample_bams, vcf_file);
                }

                array_list_clear(sample_list, vcf_annot_sample_free);

                output_item = list_item_new(i, 0, batch);
                list_insert_item(output_item, output_list);

                free(chunk_starts);
                free(chunk_sizes);
                i++;
            } // End While


            stop = omp_get_wtime();
            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            if (sample_bams) { kh_destroy(bams, sample_bams); }

            // Decrease list writers count
            for (int i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        array_list_free(sample_list, vcf_annot_sample_free);
        } // end section

#pragma omp section
        {
            LOG_DEBUG_F("Thread %d writes the output\n", omp_get_thread_num());
            
            start = omp_get_wtime();

            char *output_filename = shared_options_data->output_filename;
            
            char *aux_filename;
            FILE *output_fd = get_output_file(shared_options_data, aux_filename, &aux_filename);
            LOG_INFO_F("Output filename = %s\n", aux_filename);
            free(aux_filename);
            
            list_item_t *list_item = NULL;
            vcf_header_entry_t *entry;
            vcf_record_t *record;
            int *num_records;
            vcf_batch_t *batch;

            // Write the header
            write_vcf_header(vcf_file, output_fd);

            while( list_item = list_remove_item(output_list)){
                    batch = (vcf_batch_t*) list_item->data_p;
                    for (int i = 0; i < batch->records->size; i++) {
                        record = (vcf_record_t*) batch->records->items[i];
                        write_vcf_record(record, output_fd);
                    }
                    vcf_batch_free(batch);
                    list_item_free(list_item);
            }
            fclose(output_fd);
            free(output_filename);
        }
    } // END Parallel sections

    list_free_deep(output_list, vcf_batch_free);

    vcf_close(vcf_file);
    return 0;
}

// annot.c
int vcf_annot_chr_cmp(const vcf_annot_chr_t *v1, const vcf_annot_chr_t *v2){
    return strcmp(v1->name, v2->name);
}

// annot.c
void vcf_annot_sample_free(vcf_annot_sample_t *sample){
    if(sample != NULL){
        array_list_free(sample->chromosomes, vcf_annot_chr_free);
        free(sample->name);
        free(sample);
    }
}

// annot.c
void vcf_annot_chr_free(vcf_annot_chr_t *chr){
    if(chr != NULL){
        array_list_free(chr->positions, vcf_annot_pos_free);
        free(chr->name);
        free(chr);
    }
}

// annot.c
void vcf_annot_pos_free(vcf_annot_pos_t *pos){
    if(pos != NULL){
        free(pos);
    }
}

// annot.c
void vcf_annot_bam_free(vcf_annot_bam_t *bam){
    if(bam != NULL)
    {
        free(bam->bam_filename);
        free(bam->bai_filename);
        free(bam);
    }
}

// annot.c
int vcf_annot_process_chunk(vcf_record_t **variants, int num_variants, array_list_t *sample_list, vcf_file_t *vcf_file){

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

        copy_buf = strndup(record->format, record->format_len);
        gt_pos = get_field_position_in_format("GT", copy_buf); 
        if (copy_buf) {
            free(copy_buf);
            copy_buf = NULL;
        }
        if (gt_pos < 0) { 
            continue; 
        }   // This variant has no GT field

        // TODO aaleman:d Comprobar que las samples no son todas missings, en caso de que sí lo
        // sean saltar a la iteración siguiente
        for (int n = 0; n < array_list_size(sample_list); n++){
            copy_buf = strdup((char*) array_list_get(n, record->samples)); 
            alleles_code = get_alleles(copy_buf, gt_pos, &allele1, &allele2);
            if (copy_buf) {
                free(copy_buf);
                copy_buf = NULL;
            }

            if (alleles_code == 3) { //   ./.
                annot_sample = array_list_get(n, sample_list);
                pos = record->position;
                chr = strndup(record->chromosome, record->chromosome_len);


                annot_chr = (vcf_annot_chr_t*) malloc(sizeof(vcf_annot_chr_t));
                annot_chr->name = chr;
#pragma omp critical
                {
                    array_pos = array_list_index_of(annot_chr, annot_sample->chromosomes);
                    if(array_pos != ULONG_MAX){
                        free(annot_chr);
                        free(chr);
                        annot_chr = array_list_get(array_pos, annot_sample->chromosomes);
                    }
                    else{
                        annot_chr->positions = array_list_new(100000, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
                        array_list_insert(annot_chr, annot_sample->chromosomes);
                    }
                    annot_pos = (vcf_annot_pos_t*) malloc(sizeof(vcf_annot_pos_t));
                    annot_pos->pos = pos;
                    array_list_insert(annot_pos, annot_chr->positions);
                }
            }
        }
    }
    return 0;
}

// annot.c
int vcf_annot_check_bams(vcf_annot_sample_t* annot_sample, khash_t(bams)* sample_bams){
    khiter_t iter;
    vcf_annot_bam_t* annot_bam;
    char *bam_filename;
    vcf_annot_chr_t *annot_chr;
    vcf_annot_pos_t *annot_pos;
    int count = 0, result, tid, beg, end;
    bam_file_t* bam_file;
    bam_header_t *bam_header;
    bam_index_t *idx = 0;
    char *query = (char*) calloc(1024, sizeof(char));


    iter = kh_get(ids, sample_bams, annot_sample->name);
    annot_bam = kh_value(sample_bams, iter);
    bam_file = bam_fopen(annot_bam->bam_filename);
    idx = bam_index_load(annot_bam->bam_filename);

    count_func_data_t count_data = { bam_file->bam_header_p, &count };

    for(int j = 0; j < array_list_size(annot_sample->chromosomes); j++){
        annot_chr = (vcf_annot_chr_t*) array_list_get(j, annot_sample->chromosomes);
        //            printf ( "%s (%d)\n", annot_chr->name, array_list_size(annot_chr->positions) );
        for(int k = 0; k < array_list_size(annot_chr->positions); k++){
            annot_pos = (vcf_annot_pos_t*) array_list_get(k, annot_chr->positions);
            sprintf(query, "%s:%d-%d", annot_chr->name, annot_pos->pos, annot_pos->pos); 
            bam_parse_region(bam_file->bam_header_p, query, &tid, &beg, &end);
            result = bam_fetch(bam_file->bam_fd, idx, tid, beg, end, &count_data, count_func);
            annot_pos->dp = *(count_data.count);
            *(count_data.count) = 0;
        }
    }
    bam_index_destroy(idx);
    bam_fclose(bam_file);

    free(bam_filename);
    free(query);
    return 0;
}

// annot.c
int vcf_annot_edit_chunk(vcf_record_t **variants, int num_variants, array_list_t *sample_list, khash_t(bams)* sample_bams, vcf_file_t *vcf_file){

    char * copy_buf;
    char * aux_buf;
    char value[1024];
    int gt_pos, i;
    int dp_pos;
    char *chr;
    size_t pos;
    vcf_record_t *record;
    vcf_annot_pos_t* annot_pos;
    vcf_annot_sample_t* annot_sample;
    
    for(int j = 0; j < num_variants; j++){
        record = variants[j];
        
        copy_buf = strndup(record->format, record->format_len);
        chr = strndup(record->chromosome, record->chromosome_len);
        pos = record->position;
        
        dp_pos = get_field_position_in_format("DP", copy_buf); 
        free(copy_buf);
        
        copy_buf = strndup(record->format, record->format_len);
        gt_pos = get_field_position_in_format("GT", copy_buf); 
        free(copy_buf);
        
        if (dp_pos < 0 || gt_pos < 0) { 
            continue; 
        }   // This variant has no GT/DP field

        for (int n = 0; n < array_list_size(sample_list); n++){
            annot_sample = (vcf_annot_sample_t*) array_list_get(n, sample_list); 
            annot_pos = vcf_annot_get_pos(annot_sample->name, chr, pos, sample_list);
            if(annot_pos && annot_pos->dp > 0){
                copy_buf = (char*) array_list_get(n, record->samples); 
                sprintf(value, "%d", annot_pos->dp);
                set_field_value_in_sample(&copy_buf, dp_pos, value);
                set_field_value_in_sample(&copy_buf, gt_pos, "0/0");
            }
        }
        if(chr)
            free(chr);
    }
    return 0;
}

// vcf_util.c
void set_field_value_in_sample(char **sample, int position, char* value){
    assert(sample);
    assert(value);
    assert(position >= 0);

    int field_pos = 0;
    int num_splits;

    char **splits = split(*sample, ":", &num_splits);
    *sample = realloc(*sample, strlen(*sample) + strlen(value) + 1);
    strcpy(*sample, "");
    for(int i = 0; i < num_splits; i++)
    {
        if(i == position){
            strcat(*sample, value);
        }else{
            strcat(*sample, splits[i]);
        }
        if(i < (num_splits - 1))
            strcat(*sample, ":");
        free(splits[i]);
    }
    free(splits);
}

// annot.c
vcf_annot_pos_t * vcf_annot_get_pos(char *sample_name, char *chr, size_t pos, array_list_t *sample_list){
    vcf_annot_sample_t *annot_sample;
    vcf_annot_pos_t *annot_pos = NULL;
    vcf_annot_chr_t *annot_chr = NULL;
    int res_pos = -1;
    for (int i = 0; i < array_list_size(sample_list); i++) {
        annot_sample = (vcf_annot_sample_t*) array_list_get(i, sample_list);
        if(strcmp(annot_sample->name, sample_name) == 0){
            for (int j = 0; j < array_list_size(annot_sample->chromosomes); j++) {
                annot_chr = (vcf_annot_chr_t*) array_list_get(j, annot_sample->chromosomes);
                if(strcmp(annot_chr->name, chr) == 0){
                    res_pos = binary_search_pos(annot_chr->positions, 0 , array_list_size(annot_chr->positions) - 1, pos);
                    if(res_pos >= 0){
                        annot_pos = array_list_get(res_pos, annot_chr->positions);
                        return annot_pos;
                    }
                }
            }
        }
    }
    return NULL;
}

// annot.c
static int binary_search_pos(array_list_t* array_list , unsigned int f,unsigned int l, unsigned int target) {

    long middle, first, last;
    first = f;
    last = l;
    vcf_annot_pos_t **array = (vcf_annot_pos_t*) array_list->items;
    while (first <= last) {
        if ( array[middle]->pos < target )
            first = middle + 1;    
        else if ( array[middle]->pos == target) 
        {
            return middle;
        }
        else
            last = middle - 1;

        middle = (first + last)/2;
    }
    return -1;
}

static void vcf_annot_sort_sample(vcf_annot_sample_t* annot_sample){
    int j;
    vcf_annot_pos_t ** annot_pos_array;
    vcf_annot_chr_t *annot_chr = NULL;

    int res_pos = -1;
    for (j = 0; j < array_list_size(annot_sample->chromosomes); j++) {
        annot_chr = (vcf_annot_chr_t*) array_list_get(j, annot_sample->chromosomes);
        annot_pos_array = (vcf_annot_pos_t*) annot_chr->positions->items;
        qsort(annot_pos_array, array_list_size(annot_chr->positions), sizeof(vcf_annot_pos_t), vcf_annot_pos_cmp);
    }
}

static int vcf_annot_pos_cmp (const void * a, const void * b){
    long pos_a_val = (*(vcf_annot_pos_t **)a)->pos;  
    long pos_b_val = (*(vcf_annot_pos_t **)b)->pos; 

    return (pos_a_val - pos_b_val );
}
