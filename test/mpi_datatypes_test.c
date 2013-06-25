#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#define NUM_GENOTYPES 3

struct risky_combination {
    double accuracy;
    int order;
    int num_risky_genotypes;
/*
    uint8_t *genotypes;
    int *combination;
*/
};

int main(int argc, char *argv[]) {
    int num_mpi_ranks = 1, mpi_rank;
    int order = 2;
    int num_member_blocks = 2;
    struct risky_combination sent, received;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    
    MPI_Datatype mpi_risky_combination_type;
    
    // Length of the struct members
/*
    int lengths[num_member_blocks] = { 1, 2, pow(NUM_GENOTYPES, order), order };
*/
    int lengths[] = { 1, 2 };
    // Datatype of the struct members
/*
    MPI_Datatype types[] = { MPI_DOUBLE, MPI_INT, MPI_BYTE, MPI_INT };
*/
    MPI_Datatype types[] = { MPI_DOUBLE, MPI_INT };
    // Address and offset of the struct members
    MPI_Aint addr_struct, addr_accuracy, addr_order, addr_genotypes, addr_combination;
    MPI_Aint d_extent;
/*
    MPI_Aint displs[];
*/
    MPI_Aint offsets[num_member_blocks];
    MPI_Get_address(&sent, &addr_struct);
    MPI_Get_address(&sent.accuracy, &addr_accuracy);
    MPI_Get_address(&sent.order, &addr_order);
    MPI_Type_extent(MPI_DOUBLE, &d_extent);
    
    offsets[0] = 0;
    offsets[1] = d_extent;
    
    MPI_Type_struct(num_member_blocks, lengths, offsets, types, &mpi_risky_combination_type);
    MPI_Type_commit(&mpi_risky_combination_type);
    
    srand(133);
    
    if (mpi_rank == 0) {
        for (int i = 1; i < num_mpi_ranks; i++) {
            // Create and send a risky_combination
            sent.accuracy = (double) (rand() % 1000) / 1000;
            sent.order = rand() % 4 + 1;
            sent.num_risky_genotypes = rand() % 9 + 1;
        
            printf("%d -> %d send (%.3f, %d, %d)\n", mpi_rank, i, sent.accuracy, sent.order, sent.num_risky_genotypes);
            MPI_Send(&sent, 1, mpi_risky_combination_type, i, 0, MPI_COMM_WORLD); 
        } 
    } else {
        MPI_Status stat;
        MPI_Recv(&received, 1, mpi_risky_combination_type, 0, 0, MPI_COMM_WORLD, &stat);
        printf("-> %d recv (%.3f, %d, %d)\n", mpi_rank, received.accuracy, received.order, received.num_risky_genotypes);
    }
    
    MPI_Finalize();
    
    return 0;
}
