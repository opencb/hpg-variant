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


static void show_best_models_per_repetition(int order, int num_cv_repetitions, struct heap *best_models[]);
static khash_t(cvc) *prepare_models_for_cvc(int order, int num_cv_repetitions, int max_val_len, struct heap *best_models[]);
static int choose_best_model(int order, int num_cv_repetitions, int max_val_len, struct heap *best_models[], khash_t(cvc) *models_for_cvc, char **bestkey);


void epistasis_report(int order, size_t num_variants, int num_cv_repetitions, struct heap **best_models) {
    // Show the best model of each repetition
    show_best_models_per_repetition(order, num_cv_repetitions, best_models);

    // CVC (get the model that appears more times in the first ranking position)
    int max_val_len = log10f(num_variants);
    khash_t(cvc) *models_for_cvc = prepare_models_for_cvc(order, num_cv_repetitions, max_val_len, best_models);

    char *bestkey;
    int bestvalue = choose_best_model(order, num_cv_repetitions, max_val_len, best_models, models_for_cvc, &bestkey);

    assert(bestkey);
    LOG_INFO_F("Best model is %s with a CVC of %d/%d\n", bestkey, bestvalue, num_cv_repetitions);

    kh_destroy(cvc, models_for_cvc);
}


static void show_best_models_per_repetition(int order, int num_cv_repetitions, struct heap *best_models[]) {
    for (int r = 0; r < num_cv_repetitions; r++) {
        struct heap_node *hn = heap_peek(compare_risky_heap_max, best_models[r]);
        risky_combination *element = (risky_combination*) hn->value;
        assert(element);
        assert(element->combination);
        printf("CV %d\t(", r);
        for (int i = 0; i < order; i++) {
            printf(" %d ", element->combination[i]);
        }
        printf(") - %.3f)\n", element->accuracy);
    }
}

static khash_t(cvc) *prepare_models_for_cvc(int order, int num_cv_repetitions, int max_val_len, struct heap *best_models[]) {
    khash_t(cvc) *models_for_cvc = kh_init(cvc);
    for (int r = 0; r < num_cv_repetitions; r++) {
        struct heap_node *hn = heap_peek(compare_risky_heap_max, best_models[r]);
        risky_combination *risky = (risky_combination*) hn->value;
        
        // key = snp1_snp2_..._snpN
        char *key = calloc(order * (max_val_len + 1), sizeof(char));

        for (int i = 0; i < order-1; i++) {
            sprintf(key + strlen(key), "%d_", risky->combination[i]);
        }
        sprintf(key + strlen(key), "%d", risky->combination[order-1]);

        int ret;
        khiter_t iter = kh_get(cvc, models_for_cvc, key);
        if (iter != kh_end(models_for_cvc)) {
            (kh_value(models_for_cvc, iter))++; // Increment number of occurrences
        } else {
            iter = kh_put(cvc, models_for_cvc, key, &ret);
            if (ret) {
                kh_value(models_for_cvc, iter) = 1;
            }
        }
    }
    
    return models_for_cvc;
}

static int choose_best_model(int order, int num_cv_repetitions, int max_val_len, struct heap *best_models[], khash_t(cvc) *models_for_cvc, char **bestkey) {
    int bestvalue = 0;
    for (int k = kh_begin(models_for_cvc); k < kh_end(models_for_cvc); k++) {
        if (kh_exist(models_for_cvc, k)) {
            char *key = kh_key(models_for_cvc, k);
            int value = kh_value(models_for_cvc, k);
//             printf("%s -> %d\n", key, value);
            if (value > bestvalue) {
                *bestkey = key;
                bestvalue = value;
            } else if (value == bestvalue) {
                // If CVC(best) == CVC(candidate) ---> use CV-a
                double acc_best = 0.0f;
                double acc_candidate = 0.0f;
                
                // Sum all accuracies for the best and the candidate
                for (int r = 0; r < num_cv_repetitions; r++) {
                    struct heap_node *hn = heap_peek(compare_risky_heap_max, best_models[r]);
                    risky_combination *element = (risky_combination*) hn->value;
                    
                    // maybe_key = snp1_snp2_..._snpN
                    char *maybe_key = calloc(order * (max_val_len + 1), sizeof(char));
                    
                    for (int i = 0; i < order-1; i++) {
                        sprintf(maybe_key + strlen(maybe_key), "%d_", element->combination[i]);
                    }
                    sprintf(maybe_key + strlen(maybe_key), "%d", element->combination[order-1]);
                    
                    if (!strcmp(maybe_key, key)) {
                        acc_candidate += element->accuracy;
                    } else if (!strcmp(maybe_key, *bestkey)) {
                        acc_best += element->accuracy;
                    }
                }
                
                // Check which one is greater
                if (acc_candidate > acc_best) {
                    *bestkey = key;
                    bestvalue = value;
                }
            }
        }
    }
    
    return bestvalue;
}
