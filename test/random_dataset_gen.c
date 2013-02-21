#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


int main(int argc, char *argv[]) {
    if (argc < 4) {
        printf("Usage: %s <num_variants> <num_affected> <num_unaffected>", argv[0]);
        return 1;
    }
    
    size_t num_variants = atoi(argv[1]);
    uint32_t num_affected = atoi(argv[2]);
    uint32_t num_unaffected = atoi(argv[3]);
    int num_samples = num_affected + num_unaffected;
    
    uint8_t genotypes[num_variants * num_samples];
    
    struct timeval tv;
    gettimeofday(&tv, NULL);
    int usec = tv.tv_usec;
    srand(usec);
    
    for (int i = 0; i < num_variants; i++) {
        for (int j = 0; j < num_samples; j++) {
            genotypes[i * num_samples + j] = rand() % 3;
        }
    }
    
    printf("Dataset = {\n");
    for (int i = 0; i < num_variants; i++) {
        for (int j = 0; j < num_samples; j++) {
            printf("%d ", genotypes[i * num_samples + j]);
        }
        printf("\n");
    }
    printf("}\n\n");
    
    // Write my dataset
    char myfilename[128];
    sprintf(myfilename, "mydataset_%zu_%d_%d.bin", num_variants, num_affected, num_unaffected);
    FILE *mine = fopen(myfilename, "wb");
    
    fwrite(&num_variants, sizeof(size_t), 1, mine);
    fwrite(&num_affected, sizeof(uint32_t), 1, mine);
    fwrite(&num_unaffected, sizeof(uint32_t), 1, mine);
    fwrite(genotypes, sizeof(uint8_t), num_variants * num_samples, mine);
    
    fclose(mine);
    
    // Write original MDR dataset
    char theirfilename[128];
    sprintf(theirfilename, "theirdataset_%zu_%d_%d.txt", num_variants, num_affected, num_unaffected);
    FILE *theirs = fopen(theirfilename, "w");
    
    for (int j = 0; j < num_variants; j++) {
        fprintf(theirs, "S%d\t", j);
    }
    fprintf(theirs, "Class\n");
    
    printf("Their dataset = {\n");
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_variants; j++) {
            fprintf(theirs, "%d\t", genotypes[j * num_samples + i]);
            printf("%d ", genotypes[j * num_samples + i]);
        }
        if (i < num_affected) {
            fprintf(theirs, "1\t");
            printf("1");
        } else {
            fprintf(theirs, "0\t");
            printf("0");
        }
        
        fprintf(theirs, "\n");
        printf("\n");
    }
    printf("}\n\n");
    
    fclose(theirs);
    
    return 0;
}