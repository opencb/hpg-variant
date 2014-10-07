#include <assert.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>


void *mmap_file(size_t *len, const char *filename) {
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        perror("Error opening file");
    }
    
    struct stat st[1];
    if (fstat(fd, st)) {
        perror("Error while getting file information");
    }
    *len = (size_t) st->st_size;

    if (!*len) {
        close(fd);
        return NULL;
    }

    void *map = mmap(NULL, *len, PROT_READ, MAP_PRIVATE, fd, 0);
    if (MAP_FAILED == map) {
        perror("mmap failed");
    }
    close(fd);
    
    return map;
}

uint8_t *epistasis_dataset_load(int *num_affected, int *num_unaffected, size_t *num_variants, size_t *file_len, size_t *genotypes_offset, char *filename) {
    FILE *fp = fopen(filename, "rb");
    fread(num_variants, 1, sizeof(size_t), fp);
    fread(num_affected, 1, sizeof(uint32_t), fp);
    fread(num_unaffected, 1, sizeof(uint32_t), fp);
    fclose(fp);
    
    *genotypes_offset = sizeof(size_t) + sizeof(uint32_t) + sizeof(uint32_t);
    
    return mmap_file(file_len, filename);
}

int epistasis_dataset_close(uint8_t *contents, size_t file_len) {
    assert(contents);
    assert(file_len > 0);
    return munmap((void*) contents, file_len);
}

int main(int argc, char *argv[]) {
    int num_affected, num_unaffected, num_samples;
    size_t num_variants, file_len, genotypes_offset;
    char *filename, *output_filename;

    if (argc < 2) {
        printf("Usage: hpgvariant2mdrjava <filename>\n");
        exit(1);
    }

    filename = argv[1];
    output_filename = malloc ((strlen(filename) + 5) * sizeof(char));
    sprintf(output_filename, "%s.txt", filename);

    printf("output file is %s\n", output_filename);

    uint8_t *input_file = epistasis_dataset_load(&num_affected, &num_unaffected, &num_variants, &file_len, &genotypes_offset, filename);
    uint8_t *genotypes = input_file + genotypes_offset;
    num_samples = num_affected + num_unaffected;
    FILE *outfile = fopen(output_filename, "w");

    printf("genotypes offset is %zu\n", genotypes_offset);

    int i,j;
    for (i = 0; i < num_variants; i++) {
        fprintf(outfile, "S%d\t", i);
    }
    fprintf(outfile, "Class\n");

    for (j = 0; j < num_affected; j++) {
        for (i = 0; i < num_variants; i++) {
            uint8_t gt = genotypes[i * num_samples + j];
            if (gt == 255) gt = 3;
            fprintf(outfile, "%d\t", gt); 
        }
        fprintf(outfile, "0\n");
    }

    for (j = num_affected; j < num_samples; j++) {
        for (i = 0; i < num_variants; i++) {
            uint8_t gt = genotypes[i * num_samples + j];
            if (gt == 255) gt = 3;
            fprintf(outfile, "%d\t", gt); 
        }
        fprintf(outfile, "1\n");
    }

    fclose(outfile);
    epistasis_dataset_close(input_file, file_len);

    return 0;
}

