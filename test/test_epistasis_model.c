#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <check.h>
#include <xmmintrin.h>
#include <smmintrin.h>

#include <containers/array_list.h>

#include "gwas/epistasis/model.h"


Suite *create_test_suite(void);

int num_samples = 8;
int num_affected = 4;
int num_unaffected = 4;

uint8_t genotypes[4][8] = { { 0, 0, 1, 0, 2, 1, 0, 2 }, 
                            { 0, 1, 1, 0, 0, 0, 1, 1 },
                            { 1, 2, 0, 1, 0, 2, 0, 0 },
                            { 0, 0, 0, 2, 1, 1, 0, 2 } };

/* ******************************
 *      Unchecked fixtures      *
 * ******************************/



/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST(test_get_masks) {
    int order = 4;
    int num_samples_with_padding = 32;
    char msg[64];
    
    // Genotypes of all SNPs
    uint8_t masks_gt0[] = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt1[] = { 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt2[] = { 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt3[] = { 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *masks_genotypes[] = { masks_gt0, masks_gt1, masks_gt2, masks_gt3 };
    
    masks_info info; masks_info_init(order, 1, num_affected, num_unaffected, &info);
    
    uint8_t *masks = _mm_malloc(info.num_combinations_in_a_row * info.num_masks * sizeof(uint8_t), 16);
    set_genotypes_masks(order, masks_genotypes, 1, masks, info);
    
    fail_unless(info.num_masks == 3 * order * num_samples_with_padding, "There should be 3 GTs * 4 SNPs * 4+4 samples (padded) = 384 elements in the masks array");
    
    uint8_t mask_0[] = { 255, 255, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    // Genotype 0
                         0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    // Genotype 1
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };  // Genotype 2
    uint8_t mask_1[] = { 255, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t mask_2[] = { 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         255, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t mask_3[] = { 255, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
/*
     for (int i = 0; i < order; i++) {
         print_gt_combination(masks + i * 3 * num_samples_with_padding, i, 3 * num_samples_with_padding);
     }
     printf("\n");
     print_gt_combination(mask_0, 0, 3 * num_samples_with_padding);
     print_gt_combination(mask_1, 1, 3 * num_samples_with_padding);
     print_gt_combination(mask_2, 2, 3 * num_samples_with_padding);
     print_gt_combination(mask_3, 3, 3 * num_samples_with_padding);
*/

    for (int j = 0; j < num_samples; j++) {
        sprintf(msg, "Mask SNP0: %d should be %d\n", j, mask_0[j]);
        fail_if(masks[j] != mask_0[j], msg);
    }
    
    for (int j = 0; j < num_samples; j++) {
        sprintf(msg, "Mask SNP1: %d should be %d\n", j, mask_1[j]);
        fail_if(masks[num_samples_with_padding * 3 + j] != mask_1[j], msg);
    }
    
    for (int j = 0; j < num_samples; j++) {
        sprintf(msg, "Mask SNP2: %d should be %d\n", j, mask_2[j]);
        fail_if(masks[num_samples_with_padding * 6 + j] != mask_2[j], msg);
    }
    
    for (int j = 0; j < num_samples; j++) {
        sprintf(msg, "Mask SNP3: %d should be %d\n", j, mask_3[j]);
        fail_if(masks[num_samples_with_padding * 9 + j] != mask_3[j], msg);
    }
}
END_TEST

START_TEST(test_get_counts) {
    int order = 2;
    
    // Genotypes of all SNPs
    uint8_t masks_gt0[] = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt1[] = { 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt2[] = { 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *masks_genotypes[] = { masks_gt0, masks_gt1, masks_gt2 };
    
    masks_info info; masks_info_init(order, 1, num_affected, num_unaffected, &info);
    
    uint8_t *masks = _mm_malloc(info.num_combinations_in_a_row * info.num_masks * sizeof(uint8_t), 16);
    set_genotypes_masks(order, masks_genotypes, 1, masks, info);
    
    int num_combinations;
    uint8_t **combinations = get_genotype_combinations(order, &num_combinations);
    
    // Order 2
    int num_counts = pow(NUM_GENOTYPES, order);
    int counts_o2_aff[num_counts], counts_o2_unaff[num_counts];
    combination_counts(order, masks, combinations, num_combinations, counts_o2_aff, counts_o2_unaff, info);
    fail_if(num_counts != 9, "There must be 9 order 2 combinations");
    
    fail_if(counts_o2_aff[0] != 2 || counts_o2_unaff[0] != 0, "A(0,0) -> 2,0");
    fail_if(counts_o2_aff[1] != 1 || counts_o2_unaff[1] != 1, "A(0,1) -> 1,1");
    fail_if(counts_o2_aff[2] != 0 || counts_o2_unaff[2] != 0, "A(0,2) -> 0,0");
    fail_if(counts_o2_aff[3] != 0 || counts_o2_unaff[3] != 1, "A(1,0) -> 0,1");
    fail_if(counts_o2_aff[4] != 1 || counts_o2_unaff[4] != 0, "A(1,1) -> 1,0");
    fail_if(counts_o2_aff[5] != 0 || counts_o2_unaff[5] != 0, "A(1,2) -> 0,0");
    fail_if(counts_o2_aff[6] != 0 || counts_o2_unaff[6] != 1, "A(2,0) -> 0,1");
    fail_if(counts_o2_aff[7] != 0 || counts_o2_unaff[7] != 1, "A(2,1) -> 0,1");
    fail_if(counts_o2_aff[8] != 0 || counts_o2_unaff[8] != 0, "A(2,2) -> 0,0");
    
    free(combinations);
    
    // Order 3
    order = 3;
    combinations = get_genotype_combinations(order, &num_combinations);
    masks_info_init(order, 1, num_affected, num_unaffected, &info);
    masks = _mm_malloc(info.num_combinations_in_a_row * info.num_masks * sizeof(uint8_t), 16);
    set_genotypes_masks(order, masks_genotypes, 1, masks, info);
    num_counts = pow(NUM_GENOTYPES, order);
    int counts_o3_aff[num_counts], counts_o3_unaff[num_counts];
    combination_counts(order, masks, combinations, num_combinations, counts_o3_aff, counts_o3_unaff, info);
    fail_if(num_counts != 27, "There must be 27 order 3 combinations");
    
//     for (int i = 0; i < order; i++) {
//         print_gt_combination(masks + i * 3 * info.num_samples_per_mask, i, 3 * info.num_samples_per_mask);
//     }
    
    int idx = 0;
    fail_if(counts_o3_aff[idx] != 0     || counts_o3_unaff[idx] != 0, "A(0,0,0) -> 0,0");
    fail_if(counts_o3_aff[idx + 1] != 2 || counts_o3_unaff[idx + 1] != 0, "A(0,0,1) -> 2,0");
    fail_if(counts_o3_aff[idx + 2] != 0 || counts_o3_unaff[idx + 2] != 0, "A(0,0,2) -> 0,0");
    
    fail_if(counts_o3_aff[idx + 3] != 0 || counts_o3_unaff[idx + 3] != 1, "A(0,1,0) -> 0,1");
    fail_if(counts_o3_aff[idx + 4] != 0 || counts_o3_unaff[idx + 4] != 0, "A(0,1,1) -> 0,0");
    fail_if(counts_o3_aff[idx + 5] != 1 || counts_o3_unaff[idx + 5] != 0, "A(0,1,2) -> 1,0");
    
    fail_if(counts_o3_aff[idx + 6] != 0 || counts_o3_unaff[idx + 6] != 0, "A(0,2,0) -> 0,0");
    fail_if(counts_o3_aff[idx + 7] != 0 || counts_o3_unaff[idx + 7] != 0, "A(0,2,1) -> 0,0");
    fail_if(counts_o3_aff[idx + 8] != 0 || counts_o3_unaff[idx + 8] != 0, "A(0,2,2) -> 0,0");
    
    idx = 9;
    fail_if(counts_o3_aff[idx] != 0     || counts_o3_unaff[idx] != 0, "A(1,0,0) -> 0,0");
    fail_if(counts_o3_aff[idx + 1] != 0 || counts_o3_unaff[idx + 1] != 0, "A(1,0,1) -> 0,0");
    fail_if(counts_o3_aff[idx + 2] != 0 || counts_o3_unaff[idx + 2] != 1, "A(1,0,2) -> 0,1");
    
    fail_if(counts_o3_aff[idx + 3] != 1 || counts_o3_unaff[idx + 3] != 0, "A(1,1,0) -> 1,0");
    fail_if(counts_o3_aff[idx + 4] != 0 || counts_o3_unaff[idx + 4] != 0, "A(1,1,1) -> 0,0");
    fail_if(counts_o3_aff[idx + 5] != 0 || counts_o3_unaff[idx + 5] != 0, "A(1,1,2) -> 0,0");
    
    fail_if(counts_o3_aff[idx + 6] != 0 || counts_o3_unaff[idx + 6] != 0, "A(1,2,0) -> 0,0");
    fail_if(counts_o3_aff[idx + 7] != 0 || counts_o3_unaff[idx + 7] != 0, "A(1,2,1) -> 0,0");
    fail_if(counts_o3_aff[idx + 8] != 0 || counts_o3_unaff[idx + 8] != 0, "A(1,2,2) -> 0,0");
    
    idx = 18;
    fail_if(counts_o3_aff[idx] != 0     || counts_o3_unaff[idx] != 1, "A(2,0,0) -> 0,1");
    fail_if(counts_o3_aff[idx + 1] != 0 || counts_o3_unaff[idx + 1] != 0, "A(2,0,1) -> 0,0");
    fail_if(counts_o3_aff[idx + 2] != 0 || counts_o3_unaff[idx + 2] != 0, "A(2,0,2) -> 0,0");
    
    fail_if(counts_o3_aff[idx + 3] != 0 || counts_o3_unaff[idx + 3] != 1, "A(2,1,0) -> 0,1");
    fail_if(counts_o3_aff[idx + 4] != 0 || counts_o3_unaff[idx + 4] != 0, "A(2,1,1) -> 0,0");
    fail_if(counts_o3_aff[idx + 5] != 0 || counts_o3_unaff[idx + 5] != 0, "A(2,1,2) -> 0,0");
    
    fail_if(counts_o3_aff[idx + 6] != 0 || counts_o3_unaff[idx + 6] != 0, "A(2,2,0) -> 0,0");
    fail_if(counts_o3_aff[idx + 7] != 0 || counts_o3_unaff[idx + 7] != 0, "A(2,2,1) -> 0,0");
    fail_if(counts_o3_aff[idx + 8] != 0 || counts_o3_unaff[idx + 8] != 0, "A(2,2,2) -> 0,0");
    
    free(combinations);
    _mm_free(masks);
}
END_TEST

START_TEST(test_get_counts_all_folds_order_2) {
    int order = 2, num_folds = 5;
    int num_affected = 5, num_unaffected = 10;
    
    // Genotypes of all SNPs
    uint8_t masks_gt0[]   = { 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt1[]   = { 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt2[]   = { 1, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *masks_genotypes[] = { masks_gt0, masks_gt1, masks_gt2 };
    
    // Masks that associate each sample to its folds
    uint8_t masks_folds[] = { 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                              1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                              1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                              1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                              0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0 };
    
    masks_info info; masks_info_init(order, 1, num_affected, num_unaffected, &info);
    
    uint8_t *masks = _mm_malloc(info.num_combinations_in_a_row * info.num_masks * sizeof(uint8_t), 16);
    set_genotypes_masks(order, masks_genotypes, 1, masks, info);
    
    int num_combinations;
    uint8_t **combinations = get_genotype_combinations(order, &num_combinations);
    
    // Order 2
    int num_counts = pow(NUM_GENOTYPES, order) * num_folds;
    int counts_o2_aff[num_counts], counts_o2_unaff[num_counts];
    combination_counts_all_folds(order, masks_folds, num_folds,
                                 combinations, masks, info, 
                                 counts_o2_aff, counts_o2_unaff);
    
    int base_idx;
    
    // Fold 0
    base_idx = 0;
    fail_if(counts_o2_aff[base_idx + 0] != 2 || counts_o2_unaff[base_idx + 0] != 1, "A(0,0) -> 2,1");
    fail_if(counts_o2_aff[base_idx + 1] != 0 || counts_o2_unaff[base_idx + 1] != 1, "A(0,1) -> 0,1");
    fail_if(counts_o2_aff[base_idx + 2] != 0 || counts_o2_unaff[base_idx + 2] != 1, "A(0,2) -> 0,1");
    fail_if(counts_o2_aff[base_idx + 3] != 0 || counts_o2_unaff[base_idx + 3] != 2, "A(1,0) -> 0,2");
    fail_if(counts_o2_aff[base_idx + 4] != 1 || counts_o2_unaff[base_idx + 4] != 0, "A(1,1) -> 1,0");
    fail_if(counts_o2_aff[base_idx + 5] != 0 || counts_o2_unaff[base_idx + 5] != 0, "A(1,2) -> 0,0");
    fail_if(counts_o2_aff[base_idx + 6] != 1 || counts_o2_unaff[base_idx + 6] != 1, "A(2,0) -> 1,1");
    fail_if(counts_o2_aff[base_idx + 7] != 0 || counts_o2_unaff[base_idx + 7] != 1, "A(2,1) -> 0,1");
    fail_if(counts_o2_aff[base_idx + 8] != 0 || counts_o2_unaff[base_idx + 8] != 1, "A(2,2) -> 0,1");
    
    // Fold 1
    base_idx += 9;
    fail_if(counts_o2_aff[base_idx + 0] != 2 || counts_o2_unaff[base_idx + 0] != 0, "A(0,0) -> 2,0");
    fail_if(counts_o2_aff[base_idx + 1] != 1 || counts_o2_unaff[base_idx + 1] != 1, "A(0,1) -> 1,1");
    fail_if(counts_o2_aff[base_idx + 2] != 0 || counts_o2_unaff[base_idx + 2] != 1, "A(0,2) -> 0,1");
    fail_if(counts_o2_aff[base_idx + 3] != 0 || counts_o2_unaff[base_idx + 3] != 2, "A(1,0) -> 0,2");
    fail_if(counts_o2_aff[base_idx + 4] != 0 || counts_o2_unaff[base_idx + 4] != 0, "A(1,1) -> 0,0");
    fail_if(counts_o2_aff[base_idx + 5] != 0 || counts_o2_unaff[base_idx + 5] != 0, "A(1,2) -> 0,0");
    fail_if(counts_o2_aff[base_idx + 6] != 1 || counts_o2_unaff[base_idx + 6] != 2, "A(2,0) -> 1,2");
    fail_if(counts_o2_aff[base_idx + 7] != 0 || counts_o2_unaff[base_idx + 7] != 1, "A(2,1) -> 0,1");
    fail_if(counts_o2_aff[base_idx + 8] != 0 || counts_o2_unaff[base_idx + 8] != 1, "A(2,2) -> 0,1");
    
    // Fold 2
    base_idx += 9;
    fail_if(counts_o2_aff[base_idx + 0] != 1 || counts_o2_unaff[base_idx + 0] != 1, "A(0,0) -> 1,1");
    fail_if(counts_o2_aff[base_idx + 1] != 1 || counts_o2_unaff[base_idx + 1] != 1, "A(0,1) -> 1,1");
    fail_if(counts_o2_aff[base_idx + 2] != 0 || counts_o2_unaff[base_idx + 2] != 0, "A(0,2) -> 0,0");
    fail_if(counts_o2_aff[base_idx + 3] != 0 || counts_o2_unaff[base_idx + 3] != 3, "A(1,0) -> 0,3");
    fail_if(counts_o2_aff[base_idx + 4] != 1 || counts_o2_unaff[base_idx + 4] != 0, "A(1,1) -> 1,0");
    fail_if(counts_o2_aff[base_idx + 5] != 0 || counts_o2_unaff[base_idx + 5] != 0, "A(1,2) -> 0,0");
    fail_if(counts_o2_aff[base_idx + 6] != 1 || counts_o2_unaff[base_idx + 6] != 2, "A(2,0) -> 1,2");
    fail_if(counts_o2_aff[base_idx + 7] != 0 || counts_o2_unaff[base_idx + 7] != 1, "A(2,1) -> 0,1");
    fail_if(counts_o2_aff[base_idx + 8] != 0 || counts_o2_unaff[base_idx + 8] != 0, "A(2,2) -> 0,0");
    
    // Fold 3
    base_idx += 9;
    fail_if(counts_o2_aff[base_idx + 0] != 2 || counts_o2_unaff[base_idx + 0] != 1, "A(0,0) -> 2,1");
    fail_if(counts_o2_aff[base_idx + 1] != 1 || counts_o2_unaff[base_idx + 1] != 0, "A(0,1) -> 1,0");
    fail_if(counts_o2_aff[base_idx + 2] != 0 || counts_o2_unaff[base_idx + 2] != 1, "A(0,2) -> 0,1");
    fail_if(counts_o2_aff[base_idx + 3] != 0 || counts_o2_unaff[base_idx + 3] != 3, "A(1,0) -> 0,3");
    fail_if(counts_o2_aff[base_idx + 4] != 1 || counts_o2_unaff[base_idx + 4] != 0, "A(1,1) -> 1,0");
    fail_if(counts_o2_aff[base_idx + 5] != 0 || counts_o2_unaff[base_idx + 5] != 0, "A(1,2) -> 0,0");
    fail_if(counts_o2_aff[base_idx + 6] != 0 || counts_o2_unaff[base_idx + 6] != 2, "A(2,0) -> 0,2");
    fail_if(counts_o2_aff[base_idx + 7] != 0 || counts_o2_unaff[base_idx + 7] != 0, "A(2,1) -> 0,0");
    fail_if(counts_o2_aff[base_idx + 8] != 0 || counts_o2_unaff[base_idx + 8] != 1, "A(2,2) -> 0,1");
    
    // Fold 4
    base_idx += 9;
    fail_if(counts_o2_aff[base_idx + 0] != 1 || counts_o2_unaff[base_idx + 0] != 1, "A(0,0) -> 1,1");
    fail_if(counts_o2_aff[base_idx + 1] != 1 || counts_o2_unaff[base_idx + 1] != 1, "A(0,1) -> 1,1");
    fail_if(counts_o2_aff[base_idx + 2] != 0 || counts_o2_unaff[base_idx + 2] != 1, "A(0,2) -> 0,1");
    fail_if(counts_o2_aff[base_idx + 3] != 0 || counts_o2_unaff[base_idx + 3] != 2, "A(1,0) -> 0,2");
    fail_if(counts_o2_aff[base_idx + 4] != 1 || counts_o2_unaff[base_idx + 4] != 0, "A(1,1) -> 1,0");
    fail_if(counts_o2_aff[base_idx + 5] != 0 || counts_o2_unaff[base_idx + 5] != 0, "A(1,2) -> 0,0");
    fail_if(counts_o2_aff[base_idx + 6] != 1 || counts_o2_unaff[base_idx + 6] != 1, "A(2,0) -> 1,1");
    fail_if(counts_o2_aff[base_idx + 7] != 0 || counts_o2_unaff[base_idx + 7] != 1, "A(2,1) -> 0,1");
    fail_if(counts_o2_aff[base_idx + 8] != 0 || counts_o2_unaff[base_idx + 8] != 1, "A(2,2) -> 0,1");
    
    free(combinations);
}
END_TEST
    
START_TEST(test_get_counts_all_folds_order_3) {
    int order = 3, num_folds = 5;
    int num_affected = 5, num_unaffected = 10;
    
    // Genotypes of all SNPs
    uint8_t masks_gt0[]   = { 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt1[]   = { 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0 };
    uint8_t masks_gt2[]   = { 1, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *masks_genotypes[] = { masks_gt0, masks_gt1, masks_gt2 };
    
    // Masks that associate each sample to its folds
    uint8_t masks_folds[] = { 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                              1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                              1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                              1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                              0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0 };
    
    masks_info info; masks_info_init(order, 1, num_affected, num_unaffected, &info);
    
    uint8_t *masks = _mm_malloc(info.num_combinations_in_a_row * info.num_masks * sizeof(uint8_t), 16);
    set_genotypes_masks(order, masks_genotypes, 1, masks, info);
    
    int num_combinations;
    uint8_t **combinations = get_genotype_combinations(order, &num_combinations);
    
    int num_counts = pow(NUM_GENOTYPES, order) * num_folds;
    int counts_o3_aff[num_counts], counts_o3_unaff[num_counts];
    combination_counts_all_folds(order, masks_folds, num_folds,
                                 combinations, masks, info, 
                                 counts_o3_aff, counts_o3_unaff);
    
    // Some elements from fold 0
    fail_if(counts_o3_aff[0] != 0 || counts_o3_unaff[0] != 0, "A(0,0,0) -> 0,0");
    fail_if(counts_o3_aff[1] != 2 || counts_o3_unaff[1] != 1, "A(0,0,1) -> 2,1");
    fail_if(counts_o3_aff[2] != 0 || counts_o3_unaff[2] != 0, "A(0,0,2) -> 0,0");
    fail_if(counts_o3_aff[4] != 0 || counts_o3_unaff[4] != 1, "A(0,1,1) -> 0,1");
    fail_if(counts_o3_aff[5] != 0 || counts_o3_unaff[5] != 0, "A(0,1,2) -> 0,0");
    fail_if(counts_o3_aff[8] != 0 || counts_o3_unaff[8] != 1, "A(0,2,2) -> 0,1");
    
    fail_if(counts_o3_aff[9]  != 0 || counts_o3_unaff[9]  != 0, "A(1,0,0) -> 0,0");
    fail_if(counts_o3_aff[11] != 0 || counts_o3_unaff[11] != 1, "A(1,0,2) -> 0,1");
    fail_if(counts_o3_aff[12] != 1 || counts_o3_unaff[12] != 0, "A(1,1,0) -> 1,0");
    fail_if(counts_o3_aff[15] != 0 || counts_o3_unaff[15] != 0, "A(1,2,0) -> 0,0");
    
    fail_if(counts_o3_aff[18] != 0 || counts_o3_unaff[18] != 0, "A(2,0,0) -> 0,0");
    fail_if(counts_o3_aff[19] != 1 || counts_o3_unaff[19] != 1, "A(2,0,1) -> 1,1");
    fail_if(counts_o3_aff[21] != 0 || counts_o3_unaff[21] != 1, "A(2,1,0) -> 0,1");
    fail_if(counts_o3_aff[24] != 0 || counts_o3_unaff[24] != 1, "A(2,2,0) -> 0,1");
    
    // Some elements from fold 1
    int base_idx = 27;
    fail_if(counts_o3_aff[base_idx + 0] != 0 || counts_o3_unaff[base_idx + 0] != 0, "A(0,0,0) -> 0,0");
    fail_if(counts_o3_aff[base_idx + 1] != 2 || counts_o3_unaff[base_idx + 1] != 0, "A(0,0,1) -> 2,0");
    fail_if(counts_o3_aff[base_idx + 2] != 0 || counts_o3_unaff[base_idx + 2] != 0, "A(0,0,2) -> 0,0");
    fail_if(counts_o3_aff[base_idx + 4] != 0 || counts_o3_unaff[base_idx + 4] != 1, "A(0,1,1) -> 0,1");
    fail_if(counts_o3_aff[base_idx + 5] != 1 || counts_o3_unaff[base_idx + 5] != 0, "A(0,1,2) -> 1,0");
    fail_if(counts_o3_aff[base_idx + 8] != 0 || counts_o3_unaff[base_idx + 8] != 1, "A(0,2,2) -> 0,1");
    
    fail_if(counts_o3_aff[base_idx + 9]  != 0 || counts_o3_unaff[base_idx + 9]  != 0, "A(1,0,0) -> 0,0");
    fail_if(counts_o3_aff[base_idx + 11] != 0 || counts_o3_unaff[base_idx + 11] != 2, "A(1,0,2) -> 0,2");
    fail_if(counts_o3_aff[base_idx + 12] != 0 || counts_o3_unaff[base_idx + 12] != 0, "A(1,1,0) -> 0,0");
    fail_if(counts_o3_aff[base_idx + 15] != 0 || counts_o3_unaff[base_idx + 15] != 0, "A(1,2,0) -> 0,0");
    
    fail_if(counts_o3_aff[base_idx + 18] != 0 || counts_o3_unaff[base_idx + 18] != 1, "A(2,0,0) -> 0,1");
    fail_if(counts_o3_aff[base_idx + 19] != 1 || counts_o3_unaff[base_idx + 19] != 1, "A(2,0,1) -> 1,1");
    fail_if(counts_o3_aff[base_idx + 21] != 0 || counts_o3_unaff[base_idx + 21] != 1, "A(2,1,0) -> 0,1");
    fail_if(counts_o3_aff[base_idx + 24] != 0 || counts_o3_unaff[base_idx + 24] != 1, "A(2,2,0) -> 0,1");
    
    free(combinations);
}
END_TEST


START_TEST(test_get_confusion_matrix) {
    int order = 2;
    int num_combinations;
    masks_info info;
    int matrix[4];
    
    // ---------- order 2 ------------
    
    uint8_t **possible_2d = get_genotype_combinations(order, &num_combinations);
    // Risky combinations: (1,0), (2,1), (2,2)
    masks_info_init(order, 1, 7, 5, &info);
    risky_combination *combination_2d = risky_combination_new(order, (int[2]){ 0, 1 }, possible_2d, 3, (int[3]){ 3, 7, 8 }, NULL, info);
    
    // 7 affected, 5 unaffected
    uint8_t gt2_0a[] = { 1, 1, 0, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t gt2_1a[] = { 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *genotypes_2da[2] = { gt2_0a, gt2_1a };
    uint8_t fold_masks_2da[] = { 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    confusion_matrix(order, combination_2d, genotypes_2da, fold_masks_2da, TRAINING, (int[2]) { 7, 5 }, (int[2]) { 0, 0 }, info, matrix);
    
    // printf("matrix = { %d, %d, %d, %d }\n", matrix[0], matrix[1], matrix[2], matrix[3]);
    fail_if(matrix[0] != 6, "(7 aff,5 unaff) TP = 6");
    fail_if(matrix[1] != 1, "(7 aff,5 unaff) FN = 1");
    fail_if(matrix[2] != 1, "(7 aff,5 unaff) FP = 1");
    fail_if(matrix[3] != 4, "(7 aff,5 unaff) TN = 4");
    
    // 4 affected, 8 unaffected
    masks_info_init(order, 1, 4, 8, &info);
    uint8_t gt2_0b[] = { 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t gt2_1b[] = { 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *genotypes_2db[2] = { gt2_0b, gt2_1b };
    uint8_t fold_masks_2db[] = { 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    confusion_matrix(order, combination_2d, genotypes_2db, fold_masks_2db, TRAINING, (int[2]) { 4, 8 }, (int[2]) { 0, 0 }, info, matrix);
    
    // printf("matrix = { %d, %d, %d, %d }\n", matrix[0], matrix[1], matrix[2], matrix[3]);
    fail_if(matrix[0] != 3, "(4 aff,8 unaff) TP = 3");
    fail_if(matrix[1] != 1, "(4 aff,8 unaff) FN = 1");
    fail_if(matrix[2] != 4, "(4 aff,8 unaff) FP = 4");
    fail_if(matrix[3] != 4, "(4 aff,8 unaff) TN = 4");
    
    risky_combination_free(combination_2d);
    
    // ---------- order 3 ------------
    order = 3;
    uint8_t gt3_0[] = { 1, 1, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t gt3_1[] = { 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t gt3_2[] = { 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *genotypes_3d[3] = { gt3_0, gt3_1, gt3_2 };
    uint8_t fold_masks_3d[] = { 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    uint8_t **possible_3d = get_genotype_combinations(order, &num_combinations);
    // Risky combinations: (1,0,1), (2,1,0), (2,2,1)
    masks_info_init(order, 1, 6, 6, &info);
    risky_combination *combination_3d = risky_combination_new(order, (int[3]){ 0, 1, 2 }, possible_3d, 4, (int[4]){ 4, 10, 21, 25 }, NULL, info);
    
    // 6 affected, 6 unaffected
    confusion_matrix(order, combination_3d, genotypes_3d, fold_masks_3d, TRAINING, (int[2]) { 6, 6 }, (int[2]) { 0, 0 }, info, matrix);
    
    // printf("matrix = { %d, %d, %d, %d }\n", matrix[0], matrix[1], matrix[2], matrix[3]);
    fail_if(matrix[0] != 6, "(6 aff,6 unaff) TP = 6");
    fail_if(matrix[1] != 0, "(6 aff,6 unaff) FN = 0");
    fail_if(matrix[2] != 3, "(6 aff,6 unaff) FP = 3");
    fail_if(matrix[3] != 3, "(6 aff,6 unaff) TN = 3");
    
    risky_combination_free(combination_3d);
}
END_TEST


START_TEST(test_get_confusion_matrix_excluding_samples) {
    int order = 2;
    int num_combinations;
    masks_info info;
    int matrix[4];
    
    // ---------- order 2 ------------
    
    uint8_t **possible_2d = get_genotype_combinations(order, &num_combinations);
    // Risky combinations: (1,0), (2,1), (2,2)
    masks_info_init(order, 1, 7, 5, &info);
    risky_combination *combination_2d = risky_combination_new(order, (int[2]){ 0, 1 }, possible_2d, 3, (int[3]){ 3, 7, 8 }, NULL, info);
    
    // 7 affected, 5 unaffected
    uint8_t gt2_0a[] = { 1, 1, 0, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t gt2_1a[] = { 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t *genotypes_2d[2] = { gt2_0a, gt2_1a };
    uint8_t fold_masks_2da[] = { 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    confusion_matrix(order, combination_2d, genotypes_2d, fold_masks_2da, TRAINING, (int[2]) { 4, 3 }, (int[2]) { 3, 2 }, info, matrix);
    
    printf("matrix = { %d, %d, %d, %d }\n", matrix[0], matrix[1], matrix[2], matrix[3]);
    fail_if(matrix[0] != 3, "(7 aff,5 unaff) TP = 3");
    fail_if(matrix[1] != 1, "(7 aff,5 unaff) FN = 1");
    fail_if(matrix[2] != 0, "(7 aff,5 unaff) FP = 0");
    fail_if(matrix[3] != 3, "(7 aff,5 unaff) TN = 3");
    
    confusion_matrix(order, combination_2d, genotypes_2d, fold_masks_2da, TESTING, (int[2]) { 4, 3 }, (int[2]) { 3, 2 }, info, matrix);
    
    printf("matrix = { %d, %d, %d, %d }\n", matrix[0], matrix[1], matrix[2], matrix[3]);
    fail_if(matrix[0] != 3, "(7 aff,5 unaff) TP = 3");
    fail_if(matrix[1] != 0, "(7 aff,5 unaff) FN = 0");
    fail_if(matrix[2] != 1, "(7 aff,5 unaff) FP = 1");
    fail_if(matrix[3] != 1, "(7 aff,5 unaff) TN = 1");
    
    uint8_t fold_masks_2db[] = { 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    confusion_matrix(order, combination_2d, genotypes_2d, fold_masks_2db, TRAINING, (int[2]) { 4, 2 }, (int[2]) { 3, 3 }, info, matrix);
    
    printf("matrix = { %d, %d, %d, %d }\n", matrix[0], matrix[1], matrix[2], matrix[3]);
    fail_if(matrix[0] != 3, "(7 aff,5 unaff) TP = 3");
    fail_if(matrix[1] != 1, "(7 aff,5 unaff) FN = 1");
    fail_if(matrix[2] != 0, "(7 aff,5 unaff) FP = 0");
    fail_if(matrix[3] != 2, "(7 aff,5 unaff) TN = 2");
    
    confusion_matrix(order, combination_2d, genotypes_2d, fold_masks_2db, TESTING, (int[2]) { 4, 2 }, (int[2]) { 3, 3 }, info, matrix);
    
    printf("matrix = { %d, %d, %d, %d }\n", matrix[0], matrix[1], matrix[2], matrix[3]);
    fail_if(matrix[0] != 3, "(7 aff,5 unaff) TP = 3");
    fail_if(matrix[1] != 0, "(7 aff,5 unaff) FN = 0");
    fail_if(matrix[2] != 1, "(7 aff,5 unaff) FP = 1");
    fail_if(matrix[3] != 2, "(7 aff,5 unaff) TN = 2");
    
    uint8_t fold_masks_2dc[] = { 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    confusion_matrix(order, combination_2d, genotypes_2d, fold_masks_2dc, TRAINING, (int[2]) { 6, 4 }, (int[2]) { 1, 1 }, info, matrix);
    
    printf("matrix = { %d, %d, %d, %d }\n", matrix[0], matrix[1], matrix[2], matrix[3]);
    fail_if(matrix[0] != 6, "(7 aff,5 unaff) TP = 6");
    fail_if(matrix[1] != 0, "(7 aff,5 unaff) FN = 0");
    fail_if(matrix[2] != 0, "(7 aff,5 unaff) FP = 0");
    fail_if(matrix[3] != 4, "(7 aff,5 unaff) TN = 4");
    
    confusion_matrix(order, combination_2d, genotypes_2d, fold_masks_2dc, TESTING, (int[2]) { 6, 4 }, (int[2]) { 1, 1 }, info, matrix);
    
    printf("matrix = { %d, %d, %d, %d }\n", matrix[0], matrix[1], matrix[2], matrix[3]);
    fail_if(matrix[0] != 0, "(7 aff,5 unaff) TP = 0");
    fail_if(matrix[1] != 1, "(7 aff,5 unaff) FN = 1");
    fail_if(matrix[2] != 1, "(7 aff,5 unaff) FP = 1");
    fail_if(matrix[3] != 0, "(7 aff,5 unaff) TN = 0");
    
    risky_combination_free(combination_2d);
}
END_TEST


START_TEST(test_model_evaluation_formulas) {
    unsigned int confusion_matrix_1[] = { 40, 2, 4, 10 };
    unsigned int confusion_matrix_2[] = { 20, 10, 10, 20 };
    
    // Accuracy
    fail_if(evaluate_model(confusion_matrix_1, CA) - 0.89285714 > 1e-6, "CA(1) = 0.89285714");
    fail_if(evaluate_model(confusion_matrix_2, CA) - 0.66666666 > 1e-6, "CA(2) = 0.6666...");
    
    // Balanced accuracy
    fail_if(evaluate_model(confusion_matrix_1, BA) - 0.83333333 > 1e-6, "BA(1) = 0.83333...");
    fail_if(evaluate_model(confusion_matrix_2, BA) - 0.66666666 > 1e-6, "BA(2) = 0.6666...");
    
    // Gamma
    fail_if(evaluate_model(confusion_matrix_1, GAMMA) - 0.96078431 > 1e-6, "GAMMA(1) = 0.96078431");
    fail_if(evaluate_model(confusion_matrix_2, GAMMA) - 0.6 > 1e-6, "GAMMA(2) = 0.6");
    
    // Kendall's Tau-B
    fail_if(evaluate_model(confusion_matrix_1, TAU_B) - 0.70352647 > 1e-6, "TAU_B(1) = 0.70352647");
    fail_if(evaluate_model(confusion_matrix_2, TAU_B) - 0.33333333 > 1e-6, "TAU_B(2) = 0.3333...");
    
}
END_TEST


/* ******************************
 *      Main entry point        *
 * ******************************/

int main (int argc, char *argv) {
    Suite *fs = create_test_suite();
    SRunner *fs_runner = srunner_create(fs);
    srunner_run_all(fs_runner, CK_NORMAL);
    int number_failed = srunner_ntests_failed (fs_runner);
    srunner_free (fs_runner);
    
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


Suite *create_test_suite(void) {
    TCase *tc_counts = tcase_create("Genotype counts");
    tcase_add_test(tc_counts, test_get_masks);
    tcase_add_test(tc_counts, test_get_counts);
    tcase_add_test(tc_counts, test_get_counts_all_folds_order_2);
    tcase_add_test(tc_counts, test_get_counts_all_folds_order_3);
    
    TCase *tc_ranking = tcase_create("Evaluation and ranking");
    tcase_add_test(tc_ranking, test_get_confusion_matrix);
    tcase_add_test(tc_ranking, test_get_confusion_matrix_excluding_samples);
    tcase_add_test(tc_ranking, test_model_evaluation_formulas);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Epistasis model");
    suite_add_tcase(fs, tc_counts);
    suite_add_tcase(fs, tc_ranking);
    
    return fs;
}
