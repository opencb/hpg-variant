#include "epistasis.h"


struct heap* merge_rankings(int num_folds, struct heap **ranking_risky) {
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
            hn = heap_take(compare_risky_heap_min, ranking_risky[i]);
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
        risky_combination *element = repetition_ranking[i];
        if (!compare_risky(&current, &element)) {
            assert(current != element);
            current->accuracy += element->accuracy;
            risky_combination_free(element);
        } else {
            current->accuracy /= num_folds;
            add_to_model_ranking(current, repetition_ranking_size, sorted_repetition_ranking, compare_risky_heap_max);
            current = element;
        }
    }
    // Don't leave last element out!
    current->accuracy /= num_folds;
    add_to_model_ranking(current, repetition_ranking_size, sorted_repetition_ranking, compare_risky_heap_max);

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
