CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

# Source code folders
LIBS_ROOT = $(PWD)/libs
CONTAINERS_DIR = $(LIBS_ROOT)/containers
COMMONS_DIR = $(LIBS_ROOT)/commons
BIOINFO_DATA_DIR = $(LIBS_ROOT)/bioformats
REGION_DIR = $(BIOINFO_DATA_DIR)/features/region

# Include and libs folders
INCLUDES = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR) -I $(BIOINFO_DATA_DIR) -I $(BIOINFO_DATA_DIR)/vcf/ -I $(BIOINFO_DATA_DIR)/gff/ -I $(BIOINFO_DATA_DIR)/ped/ -I . -I ./effect/ -I ./gwas -I /usr/include/libxml2
LIBS = -L/usr/lib/x86_64-linux-gnu -lcurl -Wl,-Bsymbolic-functions -lconfig -lcprops -fopenmp -lm -lxml2
LIBS_TEST = -lcheck

INCLUDES_STATIC = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR) -I $(BIOINFO_DATA_DIR) -I $(BIOINFO_DATA_DIR)/vcf/ -I $(BIOINFO_DATA_DIR)/gff/ -I $(BIOINFO_DATA_DIR)/ped/ -I . -I ./effect/ -I ./gwas -I /usr/include/libxml2 -I ./include
LIBS_STATIC = -Llibs -L/usr/lib/x86_64-linux-gnu -lcurl -Wl,-Bsymbolic-functions -lconfig -lcprops -fopenmp -lm -lxml2

# Source files dependencies
VCF_OBJS = $(BIOINFO_DATA_DIR)/vcf/vcf_*.o
GFF_OBJS = $(BIOINFO_DATA_DIR)/gff/gff_*.o
PED_OBJS = $(BIOINFO_DATA_DIR)/ped/ped_*.o $(BIOINFO_DATA_DIR)/family.o
REGION_TABLE_OBJS = $(REGION_DIR)/region.o $(CONTAINERS_DIR)/region_table.o $(CONTAINERS_DIR)/region_table_utils.o
MISC_OBJS = $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/result.o $(CONTAINERS_DIR)/list.o

# Project source files
EFFECT_FILES = effect/main_effect.c effect/effect_options_parsing.c effect/effect_runner.c effect/auxiliary_files_writer.c
GWAS_FILES = gwas/main_gwas.c gwas/gwas_options_parsing.c gwas/tdt.c gwas/tdt_runner.c gwas/checks_family.c
HPG_VARIANT_FILES = global_options.c hpg_variant_utils.c

# Project object files
EFFECT_OBJS = main_effect.o effect_options_parsing.o effect_runner.o auxiliary_files_writer.o
GWAS_OBJS = main_gwas.o gwas_options_parsing.o tdt.o tdt_runner.o checks_family.o
HPG_VARIANT_OBJS = global_options.o hpg_variant_utils.o

ALL_FILES = $(HPG_VARIANT_OBJS) $(EFFECT_OBJS) $(GWAS_OBJS) $(VCF_OBJS) $(GFF_OBJS) $(PED_OBJS) $(REGION_TABLE_OBJS) $(MISC_OBJS)

# Targets
all: compile-dependencies hpg-variant

deploy: compile-dependencies $(HPG_VARIANT_FILES)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c main.c $(HPG_VARIANT_FILES) $(EFFECT_FILES) $(GWAS_FILES) $(INCLUDES) $(LIBS)
	test -d bin || mkdir bin
	cp hpg-variant.cfg bin
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o bin/hpg-variant main.o $(ALL_FILES) $(INCLUDES_STATIC) $(LIBS_STATIC)
	
hpg-variant: compile-dependencies $(HPG_VARIANT_FILES)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -c main.c $(HPG_VARIANT_FILES) $(EFFECT_FILES) $(GWAS_FILES) $(INCLUDES) $(LIBS)
	test -d bin || mkdir bin
	cp hpg-variant.cfg bin
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o bin/$@ main.o $(ALL_FILES) $(INCLUDES) $(LIBS)

testing: test/test_effect_runner.c test/test_tdt_runner.c $(ALL_FILES)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o test/effect.test test/test_effect_runner.c $(ALL_FILES) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o test/tdt.test test/test_tdt_runner.c $(ALL_FILES) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o test/checks_family.test test/test_checks_family.c $(ALL_FILES) $(INCLUDES) $(LIBS) $(LIBS_TEST)

compile-dependencies:
	make family.o && \
	cd $(COMMONS_DIR) && make file_utils.o http_utils.o string_utils.o log.o result.o &&  \
	cd $(REGION_DIR) && make region.o &&  \
	cd $(CONTAINERS_DIR) && make list.o region_table.o region_table_utils.o &&  \
	cd $(BIOINFO_DATA_DIR)/gff && make compile &&  \
	cd $(BIOINFO_DATA_DIR)/ped && make compile &&  \
	cd $(BIOINFO_DATA_DIR)/vcf && make compile

family.o:
	cd $(BIOINFO_DATA_DIR) && \
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c -o $(BIOINFO_DATA_DIR)/$@ $(BIOINFO_DATA_DIR)/family.c $(INCLUDES) $(LIBS)

clean:
	rm -f *.o
	rm -f $(CONTAINERS_DIR)/*.o
	rm -f $(COMMONS_DIR)/*.o
	rm -f $(BIOINFO_DATA_DIR)/vcf/*.o
	rm -f $(BIOINFO_DATA_DIR)/gff/*.o
	rm -f $(BIOINFO_DATA_DIR)/ped/*.o
	rm -f $(REGION_DIR)/*.o
	rm -rf bin

