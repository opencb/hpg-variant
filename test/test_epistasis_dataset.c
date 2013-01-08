#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <check.h>

#include <bioformats/vcf/vcf_file_structure.h>

#include "epistasis/dataset.h"
#include "epistasis/dataset_creator.c"


Suite *create_test_suite(void);

int num_samples = 20;
int num_affected = 10;
int num_unaffected = 10;
uint8_t phenotypes[20];

int num_records = 3;
char *format = "GT";
char *possible_gts[3] = { "0/0", "0/1", "1/1" };
vcf_record_t **records;

// epistasis_dataset *dataset;


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/

void setup_dataset(void) {
    // Create phenotypes
    for (int i = 0; i < num_samples; i++) {
        phenotypes[i] = i % 2;
    }
    
    // Create records
    records = malloc (num_records * sizeof(vcf_record_t*));
    
    // Record 0 gets samples { 0/0, 0/1, 1/1, 0/0, ... }
    records[0] = vcf_record_new();
    set_vcf_record_chromosome("1", 1, records[0]);
    set_vcf_record_position(1, records[0]);
    set_vcf_record_format(format, 2, records[0]);
    for (int i = 0; i < num_samples - 1; i++) {
        add_vcf_record_sample(possible_gts[i % 3], 3, records[0]);
    }
    add_vcf_record_sample("./.", 3, records[0]);
    
    // Record 1 gets samples { 0/1, 1/1, 0/0, 0/1, ... }
    records[1] = vcf_record_new();
    set_vcf_record_chromosome("1", 1, records[1]);
    set_vcf_record_position(2, records[1]);
    set_vcf_record_format(format, 2, records[1]);
    for (int i = 0; i < num_samples; i++) {
        add_vcf_record_sample(possible_gts[(i + 1) % 3], 3, records[1]);
    }
    // Record 2 gets samples { 1/1, 0/0, 0/1, 1/1, ... }
    records[2] = vcf_record_new();
    set_vcf_record_chromosome("1", 1, records[2]);
    set_vcf_record_position(3, records[2]);
    set_vcf_record_format(format, 2, records[2]);
    for (int i = 0; i < num_samples; i++) {
        add_vcf_record_sample(possible_gts[(i + 2) % 3], 3, records[2]);
    }
    
    // Create dataset
//     dataset = epistasis_dataset_new();

}

void teardown_dataset(void) {
//     epistasis_dataset_free(dataset);
    free(records);
}


/* ******************************
 *          Unit tests          *
 * ******************************/

START_TEST (test_destination) {
    int *destination = group_individuals_by_phenotype(phenotypes, num_affected, num_unaffected);
    
    for (int i = 0; i < num_samples; i++) {
//         printf("phenotype[%d] = %d\tdestination[%d] = %d\n", i, phenotypes[i], i, destination[i]);
        if (i % 2 > 0) {
            // Example: sample 1 must be destination[1] = (1 + 1) / 2 - 1 = 1 - 1 = 0
            // Example: sample 5 must be destination[5] = (5 + 1) / 2 - 1 = 3 - 1 = 2
            fail_unless(destination[i] == (i + 1) / 2 - 1, "Cases must be placed in the index = (i + 1) / 2 - 1");
        } else {
            // Example: sample 2 must be destination[2] = num_affected + ceil((2 + 1) / 2) - 1 = 10 + ceil(1.5) - 1 = 11
            // Example: sample 4 must be destination[4] = num_affected + ceil((4 + 1) / 2) - 1 = 10 + ceil(2.5) - 1 = 12
            fail_unless(destination[i] == num_affected + (i + 1) / 2, 
                        "Controls must be placed in the index = num_affected + (i + 1) / 2");
        }
    }
}
END_TEST

START_TEST (test_process_records) {
    int *destination = group_individuals_by_phenotype(phenotypes, num_affected, num_unaffected);
    uint8_t *genotypes = epistasis_dataset_process_records(records, num_records, destination, num_samples);
    
    fail_unless(genotypes[0] == 1, "Variant 1:1: Sample 0 must have genotype 0/1");
    fail_unless(genotypes[1] == 0, "Variant 1:1: Sample 1 must have genotype 0/0");
    fail_unless(genotypes[2] == 2, "Variant 1:1: Sample 2 must have genotype 1/1");
    fail_unless(genotypes[3] == 1, "Variant 1:1: Sample 3 must have genotype 0/1");
    fail_unless(genotypes[4] == 0, "Variant 1:1: Sample 4 must have genotype 0/0");
    fail_unless(genotypes[5] == 2, "Variant 1:1: Sample 5 must have genotype 1/1");
    fail_unless(genotypes[6] == 1, "Variant 1:1: Sample 6 must have genotype 0/1");
    fail_unless(genotypes[7] == 0, "Variant 1:1: Sample 7 must have genotype 0/0");
    fail_unless(genotypes[8] == 2, "Variant 1:1: Sample 8 must have genotype 1/1");
    fail_unless(genotypes[9] == 255, "Variant 1:1: Sample 9 must have genotype ./.");
    fail_unless(genotypes[10] == 0, "Variant 1:1: Sample 10 must have genotype 0/0");
    fail_unless(genotypes[11] == 2, "Variant 1:1: Sample 11 must have genotype 1/1");
    fail_unless(genotypes[18] == 1, "Variant 1:1: Sample 18 must have genotype 0/1");
    fail_unless(genotypes[19] == 0, "Variant 1:1: Sample 19 must have genotype 0/0");
    
    fail_unless(genotypes[num_samples + 0] == 2, "Variant 1:2: Sample 0 must have genotype 1/1");
    fail_unless(genotypes[num_samples + 1] == 1, "Variant 1:2: Sample 1 must have genotype 0/1");
    fail_unless(genotypes[num_samples + 2] == 0, "Variant 1:2: Sample 2 must have genotype 0/0");
    fail_unless(genotypes[num_samples + 3] == 2, "Variant 1:2: Sample 3 must have genotype 1/1");
    fail_unless(genotypes[num_samples + 4] == 1, "Variant 1:2: Sample 4 must have genotype 0/1");
    fail_unless(genotypes[num_samples + 5] == 0, "Variant 1:2: Sample 5 must have genotype 0/0");
    fail_unless(genotypes[num_samples + 6] == 2, "Variant 1:2: Sample 6 must have genotype 1/1");
    fail_unless(genotypes[num_samples + 7] == 1, "Variant 1:2: Sample 7 must have genotype 0/1");
    fail_unless(genotypes[num_samples + 8] == 0, "Variant 1:2: Sample 8 must have genotype 0/0");
    fail_unless(genotypes[num_samples + 9] == 2, "Variant 1:2: Sample 9 must have genotype 1/1");
    fail_unless(genotypes[num_samples + 10] == 1, "Variant 1:2: Sample 10 must have genotype 0/1");
    fail_unless(genotypes[num_samples + 11] == 0, "Variant 1:2: Sample 11 must have genotype 0/0");
    fail_unless(genotypes[num_samples + 18] == 2, "Variant 1:2: Sample 18 must have genotype 1/1");
    fail_unless(genotypes[num_samples + 19] == 1, "Variant 1:2: Sample 19 must have genotype 0/1");
    
    fail_unless(genotypes[2 * num_samples + 0] == 0, "Variant 1:3: Sample 0 must have genotype 0/0");
    fail_unless(genotypes[2 * num_samples + 1] == 2, "Variant 1:3: Sample 1 must have genotype 1/1");
    fail_unless(genotypes[2 * num_samples + 2] == 1, "Variant 1:3: Sample 2 must have genotype 0/1");
    fail_unless(genotypes[2 * num_samples + 3] == 0, "Variant 1:3: Sample 3 must have genotype 0/0");
    fail_unless(genotypes[2 * num_samples + 4] == 2, "Variant 1:3: Sample 4 must have genotype 1/1");
    fail_unless(genotypes[2 * num_samples + 5] == 1, "Variant 1:3: Sample 5 must have genotype 0/1");
    fail_unless(genotypes[2 * num_samples + 6] == 0, "Variant 1:1: Sample 6 must have genotype 0/0");
    fail_unless(genotypes[2 * num_samples + 7] == 2, "Variant 1:3: Sample 7 must have genotype 1/1");
    fail_unless(genotypes[2 * num_samples + 8] == 1, "Variant 1:3: Sample 8 must have genotype 0/1");
    fail_unless(genotypes[2 * num_samples + 9] == 0, "Variant 1:3: Sample 9 must have genotype 0/0");
    fail_unless(genotypes[2 * num_samples + 10] == 2, "Variant 1:3: Sample 10 must have genotype 1/1");
    fail_unless(genotypes[2 * num_samples + 11] == 1, "Variant 1:3: Sample 11 must have genotype 0/1");
    fail_unless(genotypes[2 * num_samples + 18] == 0, "Variant 1:3: Sample 18 must have genotype 0/0");
    fail_unless(genotypes[2 * num_samples + 19] == 2, "Variant 1:3: Sample 19 must have genotype 1/1");
}
END_TEST


START_TEST (test_get_block_stride) {
    fail_unless(get_block_stride(1024, 2) == 32, "1024 operations, order 2 -> stride 32");
    fail_unless(get_block_stride(10000000, 2) == 3163, "10M operations, order 2 -> stride 3163");
    
    fail_unless(get_block_stride(1000, 3) == 10, "1000 operations, order 3 -> stride 10");
    fail_unless(get_block_stride(1024, 3) == 11, "1024 operations, order 3 -> stride 11");
    fail_unless(get_block_stride(10000, 3) == 22, "10000 operations, order 3 -> stride 22");
    
    fail_unless(get_block_stride(10000, 4) == 10, "10000 operations, order 4 -> stride 10");
    fail_unless(get_block_stride(1000000, 5) == 16, "1M operations, order 5 -> stride 16");
}
END_TEST

START_TEST (test_get_next_block) {
    int num_blocks = 4;
    int order = 2;
    int block_2coords[] = { 0, 0 }, block_3coords[] = { 0, 0, 0 };
    int num_combinations = 1;
    
    // Order 2 combinations
    while (!get_next_block(num_blocks, order, block_2coords)) { num_combinations++; }
    fail_if(num_combinations != 10, "4 blocks, order 2 -> 10 combinations");
    
    block_2coords[0] = block_2coords[1] = 0;
    get_next_block(num_blocks, order, block_2coords); fail_if(block_2coords[0] != 0 || block_2coords[1] != 1, "Block: (0,0) -> (0,1)");
    get_next_block(num_blocks, order, block_2coords); fail_if(block_2coords[0] != 0 || block_2coords[1] != 2, "Block: (0,1) -> (0,2)");
    get_next_block(num_blocks, order, block_2coords); fail_if(block_2coords[0] != 0 || block_2coords[1] != 3, "Block: (0,2) -> (0,3)");
    get_next_block(num_blocks, order, block_2coords); fail_if(block_2coords[0] != 1 || block_2coords[1] != 1, "Block: (0,3) -> (1,1)");

    get_next_block(num_blocks, order, block_2coords); fail_if(block_2coords[0] != 1 || block_2coords[1] != 2, "Block: (1,1) -> (1,2)");
    get_next_block(num_blocks, order, block_2coords); fail_if(block_2coords[0] != 1 || block_2coords[1] != 3, "Block: (1,2) -> (1,3)");
    get_next_block(num_blocks, order, block_2coords); fail_if(block_2coords[0] != 2 || block_2coords[1] != 2, "Block: (1,3) -> (2,2)");
    
    get_next_block(num_blocks, order, block_2coords); fail_if(block_2coords[0] != 2 || block_2coords[1] != 3, "Block: (2,2) -> (2,3)");
    get_next_block(num_blocks, order, block_2coords); fail_if(block_2coords[0] != 3 || block_2coords[1] != 3, "Block: (2,3) -> (3,3)");
    
    
    // Order 3 combinations
    order = 3;
    num_combinations = 1;
    
    while (!get_next_block(num_blocks, order, block_3coords)) { num_combinations++; }
    fail_if(num_combinations != 20, "4 blocks, order 3 -> 20 combinations");
    
    block_3coords[0] = block_3coords[1] = block_3coords[2] = 0;
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 0 || block_3coords[1] != 0 || block_3coords[2] != 1, "Block: (0,0,0) -> (0,0,1)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 0 || block_3coords[1] != 0 || block_3coords[2] != 2, "Block: (0,0,1) -> (0,0,2)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 0 || block_3coords[1] != 0 || block_3coords[2] != 3, "Block: (0,0,2) -> (0,0,3)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 0 || block_3coords[1] != 1 || block_3coords[2] != 1, "Block: (0,0,3) -> (0,1,1)");

    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 0 || block_3coords[1] != 1 || block_3coords[2] != 2, "Block: (0,1,1) -> (0,1,2)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 0 || block_3coords[1] != 1 || block_3coords[2] != 3, "Block: (0,1,2) -> (0,1,3)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 0 || block_3coords[1] != 2 || block_3coords[2] != 2, "Block: (0,1,3) -> (0,2,2)");
    
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 0 || block_3coords[1] != 2 || block_3coords[2] != 3, "Block: (0,2,2) -> (0,2,3)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 0 || block_3coords[1] != 3 || block_3coords[2] != 3, "Block: (0,2,3) -> (0,3,3)");
    
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 1 || block_3coords[1] != 1 || block_3coords[2] != 1, "Block: (0,3,3) -> (1,1,1)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 1 || block_3coords[1] != 1 || block_3coords[2] != 2, "Block: (1,1,1) -> (1,1,2)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 1 || block_3coords[1] != 1 || block_3coords[2] != 3, "Block: (1,1,2) -> (1,1,3)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 1 || block_3coords[1] != 2 || block_3coords[2] != 2, "Block: (1,1,3) -> (1,2,2)");

    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 1 || block_3coords[1] != 2 || block_3coords[2] != 3, "Block: (1,2,2) -> (1,2,3)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 1 || block_3coords[1] != 3 || block_3coords[2] != 3, "Block: (1,2,3) -> (1,3,3)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 2 || block_3coords[1] != 2 || block_3coords[2] != 2, "Block: (1,3,3) -> (2,2,2)");
    
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 2 || block_3coords[1] != 2 || block_3coords[2] != 3, "Block: (2,2,2) -> (2,2,3)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 2 || block_3coords[1] != 3 || block_3coords[2] != 3, "Block: (2,2,3) -> (2,3,3)");
    get_next_block(num_blocks, order, block_3coords); fail_if(block_3coords[0] != 3 || block_3coords[1] != 3 || block_3coords[2] != 3, "Block: (2,3,3) -> (3,3,3)");
}
END_TEST

START_TEST (test_get_first_combination_in_block) {
    int stride = 10;
    int order = 2;
    int block_2coords[] = { 0, 0 }, block_3coords[] = { 0, 0, 0 };
    int *first_coords;
    
    // Order 2 combinations
    first_coords = get_first_combination_in_block(order, block_2coords, stride);
    fail_if(first_coords[0] != 0 || first_coords[1] != 1, "B(0,0) -> Pos(0,1)");
    
    block_2coords[0] = 0, block_2coords[1] = 1;
    first_coords = get_first_combination_in_block(order, block_2coords, stride);
    fail_if(first_coords[0] != 0 || first_coords[1] != 10, "B(0,1) -> Pos(0,10)");
    
    block_2coords[0] = 1, block_2coords[1] = 1;
    first_coords = get_first_combination_in_block(order, block_2coords, stride);
    fail_if(first_coords[0] != 10 || first_coords[1] != 11, "B(1,1) -> Pos(10,11)");
    
    block_2coords[0] = 1, block_2coords[1] = 3;
    first_coords = get_first_combination_in_block(order, block_2coords, stride);
    fail_if(first_coords[0] != 10 || first_coords[1] != 30, "B(1,3) -> Pos(10,30)");
    
    block_2coords[0] = 2, block_2coords[1] = 2;
    first_coords = get_first_combination_in_block(order, block_2coords, stride);
    fail_if(first_coords[0] != 20 || first_coords[1] != 21, "B(2,2) -> Pos(20,21)");
    
    block_2coords[0] = 2, block_2coords[1] = 1;
    first_coords = get_first_combination_in_block(order, block_2coords, stride);
    fail_if(first_coords[0] != 20 || first_coords[1] != 21, "B(2,1) -> Pos(20,21)");
    
    // Order 3 combinations
    order = 3;
    
    first_coords = get_first_combination_in_block(order, block_3coords, stride);
    fail_if(first_coords[0] != 0 || first_coords[1] != 1 || first_coords[2] != 2, "B(0,0) -> Pos(0,1,2)");
    
    block_3coords[0] = 0, block_3coords[1] = 1, block_3coords[1] = 1;
    first_coords = get_first_combination_in_block(order, block_3coords, stride);
    fail_if(first_coords[0] != 0 || first_coords[1] != 10 || first_coords[2] != 11, "B(0,1,1) -> Pos(0,10,11)");
    
    block_3coords[0] = 0, block_3coords[1] = 1, block_3coords[1] = 1;
    first_coords = get_first_combination_in_block(order, block_3coords, stride);
    fail_if(first_coords[0] != 0 || first_coords[1] != 10 || first_coords[2] != 11, "B(0,1,1) -> Pos(0,10,11)");
    
    block_3coords[0] = 1, block_3coords[1] = 1, block_3coords[1] = 1;
    first_coords = get_first_combination_in_block(order, block_3coords, stride);
    fail_if(first_coords[0] != 10 || first_coords[1] != 11 || first_coords[2] != 12, "B(1,1,1) -> Pos(10,11,12)");
    
    block_3coords[0] = 2, block_3coords[1] = 3, block_3coords[1] = 3;
    first_coords = get_first_combination_in_block(order, block_3coords, stride);
    fail_if(first_coords[0] != 20 || first_coords[1] != 30 || first_coords[2] != 31, "B(2,3,3) -> Pos(20,30,31)");
    
    block_3coords[0] = 3, block_3coords[1] = 3, block_3coords[1] = 2;
    first_coords = get_first_combination_in_block(order, block_3coords, stride);
    fail_if(first_coords[0] != 30 || first_coords[1] != 31 || first_coords[2] != 32, "B(3,3,2) -> Pos(30,31,32)");
    
    block_3coords[0] = 3, block_3coords[1] = 3, block_3coords[1] = 3;
    first_coords = get_first_combination_in_block(order, block_3coords, stride);
    fail_if(first_coords[0] != 30 || first_coords[1] != 31 || first_coords[2] != 32, "B(3,3,3) -> Pos(30,31,32)");
}
END_TEST

START_TEST (test_get_next_combination_in_block) {
    
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
    TCase *tc_creation = tcase_create("Dataset creation");
    tcase_add_unchecked_fixture(tc_creation, setup_dataset, teardown_dataset);
    tcase_add_test(tc_creation, test_destination);
    tcase_add_test(tc_creation, test_process_records);
    
    TCase *tc_balancing = tcase_create("Work distribution and load balancing");
    tcase_add_test(tc_balancing, test_get_block_stride);
    tcase_add_test(tc_balancing, test_get_next_block);
    tcase_add_test(tc_balancing, test_get_first_combination_in_block);
    
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Epistasis dataset");
    suite_add_tcase(fs, tc_creation);
    suite_add_tcase(fs, tc_balancing);
    
    return fs;
}
