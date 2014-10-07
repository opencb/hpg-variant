#include "epistasis.h"


void process_set_of_combinations(int num_combinations, int *combs, int order, int stride, 
                                 int num_folds, uint8_t *fold_masks, int *training_sizes, int *testing_sizes,
                                 uint8_t **block_genotypes, uint8_t **genotype_permutations,
                                 uint8_t *masks, enum evaluation_subset subset, masks_info info, 
                                 compare_risky_heap_func cmp_heap_func,
                                 int *counts_aff, int *counts_unaff, unsigned int conf_matrix[4], 
                                 int max_ranking_size, struct heap **ranking_risky_local) {
    // Get genotypes of a row of combinations
    uint8_t *combination_genotypes[info.num_combinations_in_a_row * order];
    for (int c = 0; c < num_combinations; c++) {
        for (int s = 0; s < order; s++) {
            // Derive combination address from block
            combination_genotypes[c * order + s] = block_genotypes[s] +
                                                    (combs[c * order + s] % stride) * info.num_samples_with_padding;
        }
    }

    set_genotypes_masks(order, combination_genotypes, num_combinations, masks, info); // Grouped by SNP

    // Get counts for the provided genotypes
    combination_counts_all_folds(order, fold_masks, num_folds, genotype_permutations, masks, info, counts_aff, counts_unaff);

    // TODO Right now the rest of the pipeline is executed as it previously was, but for the sake of parallelization
    // it could be better to make the choose_high_risk_combinations function work over several combinations
    for (int f = 0; f < num_folds; f++) {
        // Get high risk pairs for those counts
        void *aux_info;
        unsigned int num_risky[info.num_combinations_in_a_row];
        memset(num_risky, 0, info.num_combinations_in_a_row * sizeof(int));

        int *risky_idx = choose_high_risk_combinations2(counts_aff + f * info.num_combinations_in_a_row * info.num_cell_counts_per_combination,
                                                        counts_unaff + f * info.num_combinations_in_a_row * info.num_cell_counts_per_combination,
                                                        info.num_combinations_in_a_row, info.num_cell_counts_per_combination,
                                                        info.num_affected, info.num_unaffected,
                                                        num_risky, &aux_info, mdr_high_risk_combinations2);

/*
        printf("num risky = { ");
        for (int rc = 0; rc < info.num_combinations_in_a_row; rc++) {
            printf("%d ", num_risky[rc]);
        }
        printf("}\n");

        printf("risky gts = { ");
        for (int rc = 0; rc < info.num_combinations_in_a_row * info.num_counts_per_combination; rc++) {
            printf("%d ", risky_idx[rc]);
        }
        printf("}\n");
*/

        int risky_begin_idx = 0;
        for (int rc = 0; rc < num_combinations; rc++) {
            int *comb = combs + rc * order;
            uint8_t **my_genotypes = combination_genotypes + rc * order;

            // ------------------- BEGIN get_model_from_combination_in_fold -----------------------

            risky_combination *risky_comb = NULL;

            // Filter non-risky SNP combinations
            if (num_risky > 0) {
                // Put together the info about the SNP combination and its genotype combinations
                risky_comb = risky_combination_new(order, comb, genotype_permutations,
                                                    num_risky[rc], risky_idx + risky_begin_idx,
                                                    aux_info, info);
            }

            risky_begin_idx += num_risky[rc];

            // ------------------- END get_model_from_combination_in_fold -----------------------

            if (risky_comb) {
                // Check the model against the testing dataset
                double accuracy = test_model(order, risky_comb, my_genotypes, fold_masks + f * info.num_samples_with_padding, subset, 
                                             training_sizes + 3 * f + 1, testing_sizes + 3 * f + 1, info, conf_matrix);
//               printf("*  Balanced accuracy: %.3f\n", accuracy);

                int position = add_to_model_ranking(risky_comb, max_ranking_size, ranking_risky_local[f], cmp_heap_func);

                // If not inserted it means it is not among the most risky combinations, so free it
                if (position < 0) {
                    risky_combination_free(risky_comb);
                }
            }

        }

        free(risky_idx);
    }
}


struct heap* merge_rankings(int num_folds, struct heap **ranking_risky, compare_risky_heap_func heap_min_func, compare_risky_heap_func heap_max_func) {
    size_t repetition_ranking_size = 0;
    for (int i = 0; i < num_folds; i++) {
        repetition_ranking_size += ranking_risky[i]->size;
    }
    risky_combination *repetition_ranking[repetition_ranking_size];
    size_t current_index = 0;

    for (int i = 0; i < num_folds; i++) {
        struct heap_node *hn;
        risky_combination *element = NULL;

//            printf("Ranking fold %d = {\n", i);
        while (!heap_empty(ranking_risky[i])) {
            hn = heap_take(heap_min_func, ranking_risky[i]);
            element = (risky_combination*) hn->value;
            repetition_ranking[current_index] = element;
            current_index++;
            free(hn);

//                printf("(%d ", element->combination[0]);
//                for (int s = 1; s < order; s++) {
//                    printf("%d ", element->combination[s]);
//                }
//                printf("- %.3f) ", element->accuracy);
        }
//            printf("}\n\n");
    }

    assert(current_index == repetition_ranking_size);

    // qsort by coordinates
    qsort(repetition_ranking, repetition_ranking_size, sizeof(risky_combination*), compare_risky);

    // Sum all values of each position and get the mean of accuracies
    struct heap *sorted_repetition_ranking = malloc(sizeof(struct heap)); heap_init(sorted_repetition_ranking);
    risky_combination *current = repetition_ranking[0];

    for (int i = 1; i < repetition_ranking_size; i++) {
        risky_combination *other = repetition_ranking[i];
        if (!compare_risky(&current, &other)) {
            assert(current != other);
            current->accuracy += other->accuracy;
            current->cross_validation_count += other->cross_validation_count;
            risky_combination_free(other);
        } else {
            current->accuracy /= num_folds;
            add_to_model_ranking(current, repetition_ranking_size, sorted_repetition_ranking, heap_max_func);
            current = other;
        }
    }
    // Don't leave last element out!
    current->accuracy /= num_folds;
    add_to_model_ranking(current, repetition_ranking_size, sorted_repetition_ranking, heap_max_func);

    // Save the models ranking
    return sorted_repetition_ranking;
}



/* ******************************
 *      Auxiliary functions     *
 * ******************************/

int compare_risky(const void *risky_1, const void *risky_2) {
    risky_combination *r1 = *((risky_combination**) risky_1);
    risky_combination *r2 = *((risky_combination**) risky_2);
    
    for (int i = 0; i < r1->order; i++) {
        if (r1->combination[i] < r2->combination[i]) {
            return -1;
        }
        if (r1->combination[i] > r2->combination[i]) {
            return 1;
        }
    }
    
    return 0;
}
