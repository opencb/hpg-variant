#include "mpi_epistasis_helper.h"


void risky_combination_mpi_init(risky_combination_mpi_t type) {
    // Length of each block of struct members
    type.lengths[0] = 1;
    type.lengths[1] = 2;
    
    // Offset of each block of struct members
    MPI_Aint d_extent;
    MPI_Type_extent(MPI_DOUBLE, &d_extent);
    type.offsets[0] = 0;
    type.offsets[1] = d_extent;
    
    // Type of each block of struct members
    type.types[0] = MPI_DOUBLE;
    type.types[1] = MPI_INT;
    
    MPI_Type_struct(2, type.lengths, type.offsets, type.types, &(type.datatype));
    MPI_Type_commit(&(type.datatype));
}

void risky_combination_mpi_free(risky_combination_mpi_t type) {
    MPI_Type_free(&(type.datatype));
}

void send_ranking_risky_size_mpi(struct heap *ranking, int dest, MPI_Comm comm) {
    long size_cast_for_mpi = (long) ranking->size;
    MPI_Send(&size_cast_for_mpi, 1, MPI_LONG, dest, TAG_RANKING_RISKY_SIZE, comm);
}

size_t receive_ranking_risky_size_mpi(int src, MPI_Comm comm, MPI_Status stat) {
    long received_size;
    MPI_Recv(&received_size, 1, MPI_LONG, src, TAG_RANKING_RISKY_SIZE, comm, &stat);
    return (size_t) received_size;
}

void send_risky_combination_mpi(risky_combination *risky, risky_combination_mpi_t type, int dest, MPI_Comm comm) {
    MPI_Datatype mpi_risky_combination_type;
    // Length of the struct members
    int lengths[] = { 1, 2 };
    // Datatype of the struct members
    MPI_Datatype types[] = { MPI_DOUBLE, MPI_INT };
    // Offset of the struct members
    MPI_Aint d_extent;
    MPI_Type_extent(MPI_DOUBLE, &d_extent);
    MPI_Aint offsets[] = { 0, d_extent };
    
    MPI_Type_struct(2, lengths, offsets, types, &mpi_risky_combination_type);
    MPI_Type_commit(&mpi_risky_combination_type);
    
    // Send members of built-in types
    MPI_Send(risky, 1, mpi_risky_combination_type, dest, TAG_RANKING_RISKY_ELEM, comm);
    // Send combination (int pointer)
    MPI_Send(risky->combination, risky->order, MPI_INT, dest, TAG_RANKING_RISKY_ELEM, comm);
    // Send genotypes (uint8_t pointer)
    MPI_Send(risky->genotypes, risky->num_risky_genotypes * risky->order, MPI_BYTE, dest, TAG_RANKING_RISKY_ELEM, comm);
}

risky_combination* receive_risky_combination_mpi(int order, risky_combination_mpi_t type, masks_info info, int src, MPI_Comm comm, MPI_Status stat) {
    int comb[order]; memset(comb, 0, order * sizeof(int));
    risky_combination *received = risky_combination_new(order, comb, NULL, 0, NULL, NULL, info);
    
    MPI_Datatype mpi_risky_combination_type;
    // Length of the struct members
    int lengths[] = { 1, 2 };
    // Datatype of the struct members
    MPI_Datatype types[] = { MPI_DOUBLE, MPI_INT };
    // Offset of the struct members
    MPI_Aint d_extent;
    MPI_Type_extent(MPI_DOUBLE, &d_extent);
    MPI_Aint offsets[] = { 0, d_extent };
    
    MPI_Type_struct(2, lengths, offsets, types, &mpi_risky_combination_type);
    MPI_Type_commit(&mpi_risky_combination_type);
    
    // Receive members of built-in types
    MPI_Recv(received, 1, mpi_risky_combination_type, src, TAG_RANKING_RISKY_ELEM, comm, &stat);
    // Receive combination (int pointer)
    MPI_Recv(received->combination, received->order, MPI_INT, src, TAG_RANKING_RISKY_ELEM, comm, &stat);
    // Receive genotypes (uint8_t pointer)
    MPI_Recv(received->genotypes, received->num_risky_genotypes * received->order, MPI_BYTE, src, TAG_RANKING_RISKY_ELEM, comm, &stat);
    
    return received;
}
