#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "gwas/epistasis/dataset.h"

int main(int argc, char *argv[]) {
    size_t order, num_variants, stride, num_blocks_per_dim, num_block_coords, max_num_block_coords;
    int num_mpi_ranks = 1, mpi_rank;
    
    if (argc < 3) {
        printf("Usage: mpirun --np N mpi <order> <num_variants> <stride>\nExample: mpirun --np 4 mpi 2 100000 20\n");
        exit(1);
    }
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    
    order = atoi(argv[1]);
    num_variants = atol(argv[2]);
    stride = atoi(argv[3]);
    num_blocks_per_dim = ceil((double) num_variants / stride);
    num_block_coords = 0;
    max_num_block_coords = (size_t) pow(num_blocks_per_dim, order);
    
    int *block_coords, *my_block_coords;
    int curr_idx, next_idx;
    
    if (mpi_rank == 0) {
        printf("Blocks per dim = %zu\tMax block coords = %zu\n", num_blocks_per_dim, max_num_block_coords);
        block_coords = calloc(max_num_block_coords * order, sizeof(int));
        
        // Calculate all blocks coordinates
        do {
            curr_idx = num_block_coords * order;
            next_idx = curr_idx + order;
            memcpy(block_coords + next_idx, block_coords + curr_idx, order * sizeof(int));
            curr_idx = next_idx;
            num_block_coords++;
/*
            printf("(%d %d) ", block_coords[curr_idx], block_coords[curr_idx+1]);
*/
        } while (get_next_block(num_blocks_per_dim, order, block_coords + curr_idx));
    }
    
    // Send total number of blocks to all processes
    MPI_Bcast(&num_block_coords, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Prepare arguments for scatterv
    int block_counts[num_mpi_ranks];
    int block_offsets[num_mpi_ranks];

    if (mpi_rank == 0) {
        printf("nbc = %zu\tmod = %zu\n", num_block_coords, num_block_coords % num_mpi_ranks);
    }
    
    for (int p = 0; p < num_mpi_ranks; p++) {
        if (p < num_block_coords % num_mpi_ranks) {
            block_counts[p] = order * (num_block_coords / num_mpi_ranks + 1);
        } else {
            block_counts[p] = order * (num_block_coords / num_mpi_ranks);
        }
    }

    block_offsets[0] = 0;
    for (int p = 1; p < num_mpi_ranks; p++) {
        block_offsets[p] = block_offsets[p-1] + block_counts[p-1];
    }
    
    my_block_coords = calloc(block_counts[mpi_rank] * order, sizeof(int));
    
    if (mpi_rank == 0) {
        for (int p = 0; p < num_mpi_ranks; p++) {
            printf("(%d, off %d) ", block_counts[p], block_offsets[p]);
        }
        printf("\n");
    }
    
    // MPI_Scatterv sends block coordinates to all processes
    MPI_Scatterv(block_coords, block_counts, block_offsets, MPI_INT, 
                 my_block_coords, block_counts[mpi_rank], MPI_INT,
                 0, MPI_COMM_WORLD);
    
    printf("I'm MPI process %d of %d and have to work with %zu blocks\n", mpi_rank, num_mpi_ranks, block_counts[mpi_rank] / order);
    
    for (int i = 0; i < block_counts[mpi_rank] / order; i++) {
        printf("*%d* (%d %d)\n", mpi_rank, my_block_coords[i * order], my_block_coords[i * order + 1]);
    }
    
    
    if (mpi_rank == 0) {
        free(block_coords);
    }
    free(my_block_coords);
    
    MPI_Finalize();
    
    return 0;
}
