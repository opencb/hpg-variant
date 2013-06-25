#include <stdint.h>
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
    uint8_t *genotypes;
    int *combination;
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
    MPI_Status stat;
    
    // Length of the struct members
    int lengths[] = { 1, 2 };
    // Datatype of the struct members
    MPI_Datatype types[] = { MPI_DOUBLE, MPI_INT };
    // Offset of the struct members
    MPI_Aint offsets[num_member_blocks];
    MPI_Aint d_extent;
    
    MPI_Type_extent(MPI_DOUBLE, &d_extent);
    offsets[0] = 0;
    offsets[1] = d_extent;
    
    MPI_Type_struct(num_member_blocks, lengths, offsets, types, &mpi_risky_combination_type);
    MPI_Type_commit(&mpi_risky_combination_type);
    
    srand(133);
    
    if (mpi_rank == 0) {
        sent.combination = malloc (order * sizeof(int));
        sent.genotypes = malloc (pow(NUM_GENOTYPES, order) * sizeof(uint8_t));
        for (int i = 1; i < num_mpi_ranks; i++) {
            // Create and send a risky_combination
            sent.accuracy = (double) (rand() % 1000) / 1000;
            sent.order = order;
            sent.num_risky_genotypes = rand() % 9 + 1;
            for (int j = 0; j < order; j++) {
                sent.combination[j] = rand() % 100;
            }
            for (int j = 0; j < pow(NUM_GENOTYPES, order); j++) {
                sent.genotypes[j] = rand() % NUM_GENOTYPES;
            }
        
            // Send members of built-in types
            MPI_Send(&sent, 1, mpi_risky_combination_type, i, 0, MPI_COMM_WORLD);
            // Send combination (int pointer)
            MPI_Send(sent.combination, order, MPI_INT, i, 0, MPI_COMM_WORLD);
            // Send genotypes (uint8_t pointer)
            MPI_Send(sent.genotypes, pow(NUM_GENOTYPES, order), MPI_BYTE, i, 0, MPI_COMM_WORLD);
            
            // Show sent data
            printf("%d -> %d send (%.3f, %d, %d, %d, %d)\n", mpi_rank, i, 
                    sent.accuracy, sent.order, sent.num_risky_genotypes, sent.combination[0], sent.combination[1]);
            printf("%d -> %d send GTS (%d %d %d %d %d %d %d %d)\n", mpi_rank, i,
                    sent.genotypes[0], sent.genotypes[1], sent.genotypes[2], sent.genotypes[3], sent.genotypes[4],
                    sent.genotypes[5], sent.genotypes[6], sent.genotypes[7], sent.genotypes[8]);
        } 
        free(sent.combination);
        free(sent.genotypes);
    } else {
        received.combination = malloc (order * sizeof(int));
        received.genotypes = malloc (pow(NUM_GENOTYPES, order) * sizeof(uint8_t));
        // Receive members of built-in types
        MPI_Recv(&received, 1, mpi_risky_combination_type, 0, 0, MPI_COMM_WORLD, &stat);
        // Receive combination (int pointer)
        MPI_Recv(received.combination, order, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        // Receive genotypes (uint8_t pointer)
        MPI_Recv(received.genotypes, pow(NUM_GENOTYPES, order), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &stat);
        
        // Show received data
        printf("-> %d recv (%.3f, %d, %d, %d, %d)\n", mpi_rank, 
                received.accuracy, received.order, received.num_risky_genotypes, received.combination[0], received.combination[1]);
        printf("-> %d recv GTS (%d %d %d %d %d %d %d %d)\n", mpi_rank, 
                received.genotypes[0], received.genotypes[1], received.genotypes[2], received.genotypes[3], received.genotypes[4],
                received.genotypes[5], received.genotypes[6], received.genotypes[7], received.genotypes[8]);
        
        free(received.combination);
        free(received.genotypes);
    }
    
    MPI_Finalize();
    
    return 0;
}
