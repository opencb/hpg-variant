/*
 * Copyright (c) 2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2013 Ignacio Medina (ICM-CIPF)
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

#include "epistasis.h"
#include "model.h"

static void show_cross_validation_arguments(int cv_repetition, int order, enum evaluation_mode mode, enum evaluation_subset subset, FILE *fd);
static void show_cross_validation_best_models(int order, struct heap *best_models, int max_ranking_size, compare_risky_heap_func cmp_heap_max, FILE *fd);


void epistasis_report(int order, int cv_repetition, enum evaluation_mode mode, enum evaluation_subset subset, 
                      struct heap *best_models, int max_ranking_size, compare_risky_heap_func cmp_heap_max, FILE *fd) {
    show_cross_validation_arguments(cv_repetition, order, mode, subset, fd);
    show_cross_validation_best_models(order, best_models, max_ranking_size, cmp_heap_max, fd);
}

static void show_cross_validation_arguments(int cv_repetition, int order, enum evaluation_mode mode, enum evaluation_subset subset, FILE *fd) {
    fprintf(fd, "#CROSS VALIDATION %d\n", cv_repetition+1);
    fprintf(fd, "#COMBINATIONS OF: %d SNPs\n", order);
    
    if (mode == CV_C) {
        fprintf(fd, "#EVALUATION MODE: Cross-validation consistency\n");
    } else if (mode == CV_A) {
        fprintf(fd, "#EVALUATION MODE: Cross-validation accuracy\n");
    }
    
    if (subset == TRAINING) {
        fprintf(fd, "#EVALUATION PARTITION: Training\n");
    } else if (subset == TESTING) {
        fprintf(fd, "#EVALUATION PARTITION: Testing\n");
    }
}

static void show_cross_validation_best_models(int order, struct heap *best_models, int max_ranking_size, compare_risky_heap_func cmp_heap_max, FILE *fd) {
    struct heap_node *hn;
    risky_combination *element = NULL;

    int position = 0;
    fprintf(fd, "#POSITION\tSNPs\tGENOTYPES\tCV-C\tCV-A\n");
    while (!heap_empty(best_models) && position < max_ranking_size) {
        hn = heap_take(cmp_heap_max, best_models);
        element = (risky_combination*) hn->value;
        // SNPs
        fprintf(fd, "%d\t(", position+1);
        for (int i = 0; i < order - 1; i++) {
            fprintf(fd, " %d,", element->combination[i]);
        }
        fprintf(fd, " %d )\t", element->combination[order - 1]);
        
        // Genotypes
        for (int i = 0; i < element->num_risky_genotypes; i++) {
            fprintf(fd, "(%d-", element->genotypes[i * order]);
            for (int j = 1; j < order-1; j++) {
                fprintf(fd, "%d, ", element->genotypes[i * order + j]);
            }
            fprintf(fd, "%d), ", element->genotypes[i * order + order - 1]);
        }
        
        // Counts
        fprintf(fd, "%d\t%.3f\n", element->cross_validation_count, element->accuracy);
        risky_combination_free(element);
        free(hn);
        position++;
    }
}
