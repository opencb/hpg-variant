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


static int count_func(const bam1_t *b, void *data)
{
		(*((count_func_data_t*)data)->count)++;
	return 0;
}

int vcf_annot_chr_cmp(const vcf_annot_chr_t *v1, const vcf_annot_chr_t *v2) {
    return strcmp(v1->name, v2->name);
}

void vcf_annot_sample_free(vcf_annot_sample_t *sample) {
    if(sample != NULL) {
        array_list_free(sample->chromosomes, vcf_annot_chr_free);
        free(sample->name);
        free(sample);
    }
}

void vcf_annot_chr_free(vcf_annot_chr_t *chr) {
    if(chr != NULL) {
        array_list_free(chr->positions, vcf_annot_pos_free);
        free(chr->name);
        free(chr);
    }
}

void vcf_annot_pos_free(vcf_annot_pos_t *pos) {
    if(pos != NULL) {
        free(pos);
    }
}

void vcf_annot_bam_free(vcf_annot_bam_t *bam) {
    if(bam != NULL)
    {
        free(bam->bam_filename);
        free(bam->bai_filename);
        free(bam);
    }
}

int vcf_annot_process_chunk(vcf_record_t **variants, int num_variants, array_list_t *sample_list, vcf_file_t *vcf_file) {

    int gt_pos, i;
    int allele1, allele2, alleles_code;
    vcf_record_t *record;
    vcf_annot_sample_t *annot_sample;
    vcf_annot_chr_t * annot_chr;
    vcf_annot_pos_t * annot_pos;
    char * sample_name;
    char * chr;
    char * copy_buf;
    size_t array_pos;
    khiter_t iter;
    int ret;

    for(int j = 0; j < num_variants; j++) {
        record = variants[j]; 

        copy_buf = strndup(record->format, record->format_len);
        gt_pos = get_field_position_in_format("GT", copy_buf); 
        free(copy_buf);
        
        if (gt_pos < 0) { 
            continue; 
        }   // This variant has no GT field

        // TODO aaleman:d Comprobar que las samples no son todas missings, en caso de que sí lo
        // sean saltar a la iteración siguiente
        for (int n = 0; n < array_list_size(sample_list); n++) {
            copy_buf = strdup((char*) array_list_get(n, record->samples)); 
            alleles_code = get_alleles(copy_buf, gt_pos, &allele1, &allele2);
            free(copy_buf);

            if (alleles_code == 3) { //   ./.
                annot_sample = array_list_get(n, sample_list);
                chr = strndup(record->chromosome, record->chromosome_len);
#pragma omp critical
                {
                    annot_chr = vcf_annot_get_chr(chr, annot_sample->chromosomes);
                }
                annot_pos = (vcf_annot_pos_t*) malloc(sizeof(vcf_annot_pos_t));
                annot_pos->pos = record->position;
                array_list_insert(annot_pos, annot_chr->positions);
            }
        }
    }
    return 0;
}

vcf_annot_chr_t * vcf_annot_get_chr(char* chr, array_list_t* array_chr) {
    vcf_annot_chr_t* annot_chr= NULL;
    
    for(int i = 0; i < array_list_size(array_chr); i++) {
        annot_chr = (vcf_annot_chr_t*) array_list_get(i, array_chr);
        if(strcmp(annot_chr->name, chr) == 0) {
            return annot_chr;
        }
    }

    annot_chr = (vcf_annot_chr_t*) malloc(sizeof(vcf_annot_chr_t));
    annot_chr->name = chr;
    annot_chr->positions = array_list_new(100000, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
    array_list_insert(annot_chr, array_chr);

    return annot_chr;
}

int vcf_annot_check_bams(vcf_annot_sample_t* annot_sample, khash_t(bams)* sample_bams) {
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

    for(int j = 0; j < array_list_size(annot_sample->chromosomes); j++) {
        annot_chr = (vcf_annot_chr_t*) array_list_get(j, annot_sample->chromosomes);
        for(int k = 0; k < array_list_size(annot_chr->positions); k++) {
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

int vcf_annot_edit_chunk(vcf_record_t **variants, int num_variants, array_list_t *sample_list, khash_t(bams)* sample_bams, vcf_file_t *vcf_file) {

    char * copy_buf;
    char * aux_buf;
    char value[1024];
    int gt_pos, i;
    int dp_pos;
    char *chr;
    char *info;
    size_t pos;
    vcf_record_t *record;
    vcf_annot_pos_t* annot_pos;
    vcf_annot_sample_t* annot_sample;
    
    for(int j = 0; j < num_variants; j++) {
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

        for (int n = 0; n < array_list_size(sample_list); n++) {
            annot_sample = (vcf_annot_sample_t*) array_list_get(n, sample_list); 
            annot_pos = vcf_annot_get_pos(annot_sample->name, chr, pos, sample_list);
            if(annot_pos && annot_pos->dp > 0) {
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

vcf_annot_pos_t * vcf_annot_get_pos(char *sample_name, char *chr, size_t pos, array_list_t *sample_list) {
    vcf_annot_sample_t *annot_sample;
    vcf_annot_pos_t *annot_pos = NULL;
    vcf_annot_chr_t *annot_chr = NULL;
    int res_pos = -1;
    for (int i = 0; i < array_list_size(sample_list); i++) {
        annot_sample = (vcf_annot_sample_t*) array_list_get(i, sample_list);
        if(strcmp(annot_sample->name, sample_name) == 0) {
            for (int j = 0; j < array_list_size(annot_sample->chromosomes); j++) {
                annot_chr = (vcf_annot_chr_t*) array_list_get(j, annot_sample->chromosomes);
                if(strcmp(annot_chr->name, chr) == 0) {
                    res_pos = binary_search_pos(annot_chr->positions, 0 , array_list_size(annot_chr->positions) - 1, pos);
                    if(res_pos >= 0) {
                        annot_pos = array_list_get(res_pos, annot_chr->positions);
                        return annot_pos;
                    }
                }
            }
        }
    }
    return NULL;
}

static int binary_search_pos(array_list_t* array_list , unsigned int f,unsigned int l, unsigned int target) {

    long middle, first, last;
    first = f;
    last = l;
    vcf_annot_pos_t **array = (vcf_annot_pos_t*) array_list->items;
    while (first <= last) {
        if ( array[middle]->pos < target )
            first = middle + 1;    
        else if ( array[middle]->pos == target)
            return middle;
        else
            last = middle - 1;

        middle = (first + last)/2;
    }
    return -1;
}
