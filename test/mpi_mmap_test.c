#include <assert.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <mpi.h>

int main (int argc, char *argv[]) {
    int myrank;
    MPI_File fd;
    MPI_Offset len;
    MPI_Status status;
    
    if (argc < 2) {
        printf("Usage: program <filename>\n");
        exit(1);
    }
    
    char *filename = argv[1];
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fd);
    MPI_File_get_size(fd, &len);
    
    size_t file_len = (size_t) len;
    printf("File %s length = %llu bytes / %zu bytes\n", filename, len, file_len);
    
    uint8_t *map = malloc(len * sizeof(uint8_t));
    
    MPI_File_seek(fd, 0, MPI_SEEK_SET);
    MPI_Barrier(MPI_COMM_WORLD);

    size_t num_variants;
    int num_affected, num_unaffected;
    size_t genotypes_offset = sizeof(size_t) + sizeof(uint32_t) + sizeof(uint32_t);
    
    MPI_File_read(fd, map, len, MPI_BYTE, &status);
    
    printf("%d %d %d %d\n", map[0], map[1], map[2], map[3]);
    size_t *map_size_t = (size_t*) map;
    uint32_t *map_uint32 = (uint32_t*) (map + sizeof(size_t));
    num_variants = map_size_t[0];
    num_affected = map_uint32[0];
    num_unaffected = map_uint32[1];
    printf("num variants = %zu, aff = %d, unaff = %d\n", num_variants, num_affected, num_unaffected);
    
    map = map + genotypes_offset;
    
    for (int i = 0; i < 100; i+=5) {
        printf("%d %d %d %d %d ", map[i], map[i+1], map[i+2], map[i+3], map[i+4]);
    }
    
    printf("\n");
    
    MPI_File_close(&fd);
    
    MPI_Finalize();
    
    map = map - genotypes_offset;
    free(map);
    
    return 0;
}
