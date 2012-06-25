CC = gcc
CFLAGS = -std=c99 -O3 -march=native
CFLAGS_DEBUG = -std=c99 -g

# Source code folders
INC_DIR = $(PWD)/include
LIBS_DIR = $(PWD)/libs
CONTAINERS_DIR = $(LIBS_DIR)/containers
COMMONS_DIR = $(LIBS_DIR)/commons
BIOFORMATS_DIR = $(LIBS_DIR)/bioformats
MATH_DIR = $(LIBS_DIR)/math

# Include and libs folders
INCLUDES = -I . -I $(LIBS_DIR) -I $(INC_DIR) -I /usr/include/libxml2 -I/usr/local/include
LIBS = -L/usr/lib/x86_64-linux-gnu -lcurl -Wl,-Bsymbolic-functions -lconfig -lcprops -fopenmp -lm -lxml2 -lgsl -lgslcblas -largtable2
LIBS_TEST = -lcheck

INCLUDES_STATIC = -I . -I $(LIBS_DIR) -I $(INC_DIR) -I /usr/include/libxml2 -I/usr/local/include
LIBS_STATIC = -L$(LIBS_DIR) -L/usr/lib/x86_64-linux-gnu -lcurl -Wl,-Bsymbolic-functions -lconfig -lcprops -fopenmp -lm -lxml2 -lgsl -lgslcblas -largtable2

# Source files dependencies
VCF_OBJS = $(BIOFORMATS_DIR)/vcf/vcf_*.o
GFF_OBJS = $(BIOFORMATS_DIR)/gff/gff_*.o
PED_OBJS = $(BIOFORMATS_DIR)/ped/ped_*.o
REGION_TABLE_OBJS = $(BIOFORMATS_DIR)/features/region/region.o $(CONTAINERS_DIR)/region_table.o $(CONTAINERS_DIR)/region_table_utils.o
MISC_OBJS = $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/result.o \
	$(CONTAINERS_DIR)/list.o $(BIOFORMATS_DIR)/family.o $(MATH_DIR)/fisher.o

# Project source files
EFFECT_FILES = effect/main_effect.c effect/effect_options_parsing.c effect/effect_runner.c effect/auxiliary_files_writer.c
GWAS_FILES = gwas/main_gwas.c gwas/gwas_options_parsing.c gwas/tdt.c gwas/tdt_runner.c gwas/assoc.c gwas/assoc_basic_test.c gwas/assoc_fisher_test.c gwas/assoc_runner.c gwas/checks_family.c
HPG_VARIANT_FILES = global_options.c hpg_variant_utils.c

# Project object files
EFFECT_OBJS = main_effect.o effect_options_parsing.o effect_runner.o auxiliary_files_writer.o
GWAS_OBJS = main_gwas.o gwas_options_parsing.o tdt.o tdt_runner.o assoc.o assoc_basic_test.o assoc_fisher_test.o assoc_runner.o checks_family.o
HPG_VARIANT_OBJS = global_options.o hpg_variant_utils.o

ALL_FILES = $(HPG_VARIANT_OBJS) $(EFFECT_OBJS) $(GWAS_OBJS) $(VCF_OBJS) $(GFF_OBJS) $(PED_OBJS) $(REGION_TABLE_OBJS) $(MISC_OBJS)

# Targets
all: compile-dependencies hpg-variant

deploy: compile-dependencies-static $(HPG_VARIANT_FILES)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -c main.c $(HPG_VARIANT_FILES) $(EFFECT_FILES) $(GWAS_FILES) $(INCLUDES) $(LIBS_STATIC)
	test -d bin || mkdir bin
	cp hpg-variant.cfg bin
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o bin/hpg-variant main.o $(ALL_FILES) $(INCLUDES_STATIC) $(LIBS_STATIC)
	
hpg-variant: compile-dependencies $(HPG_VARIANT_FILES)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -c main.c $(HPG_VARIANT_FILES) $(EFFECT_FILES) $(GWAS_FILES) $(INCLUDES) $(LIBS)
	test -d bin || mkdir bin
	cp hpg-variant.cfg bin
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o bin/$@ main.o $(ALL_FILES) $(INCLUDES) $(LIBS)

testing: test/test_effect_runner.c test/test_tdt_runner.c $(ALL_FILES)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o test/effect.test test/test_effect_runner.c $(ALL_FILES) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o test/tdt.test test/test_tdt_runner.c $(ALL_FILES) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o test/checks_family.test test/test_checks_family.c $(ALL_FILES) $(INCLUDES) $(LIBS) $(LIBS_TEST)

compile-dependencies:
	make family.o && \
	cd $(MATH_DIR) && make compile &&  \
	cd $(COMMONS_DIR) && make compile &&  \
	cd $(CONTAINERS_DIR) && make compile &&  \
	cd $(BIOFORMATS_DIR)/features/region && make region.o &&  \
	cd $(BIOFORMATS_DIR)/gff && make compile &&  \
	cd $(BIOFORMATS_DIR)/ped && make compile &&  \
	cd $(BIOFORMATS_DIR)/vcf && make compile

compile-dependencies-static:
	cd $(MATH_DIR) && make compile &&  \
	cd $(COMMONS_DIR) && make compile &&  \
	cd $(CONTAINERS_DIR) && make compile-static &&  \
	cd $(BIOFORMATS_DIR)/features/region && make region.o &&  \
	cd $(BIOFORMATS_DIR)/gff && make compile-static &&  \
	cd $(BIOFORMATS_DIR)/ped && make compile-static &&  \
	cd $(BIOFORMATS_DIR)/vcf && make compile-static

family.o:
	cd $(BIOFORMATS_DIR) && \
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c -o $(BIOFORMATS_DIR)/$@ $(BIOFORMATS_DIR)/family.c $(INCLUDES) $(LIBS)

clean:
	rm -f *.o
	rm -f $(COMMONS_DIR)/*.o
	rm -f $(CONTAINERS_DIR)/*.o
	rm -f $(BIOFORMATS_DIR)/vcf/*.o
	rm -f $(BIOFORMATS_DIR)/gff/*.o
	rm -f $(BIOFORMATS_DIR)/ped/*.o
	rm -f $(BIOFORMATS_DIR)/features/region/*.o
	rm -rf bin

