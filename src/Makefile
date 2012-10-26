# Source code folders
ROOT_DIR = $(PWD)/..
INC_DIR = $(ROOT_DIR)/include
LIBS_DIR = $(ROOT_DIR)/libs
SRC_DIR = $(ROOT_DIR)/src
BIN_DIR = $(ROOT_DIR)/bin
TEST_DIR = $(ROOT_DIR)/test

# Libraries folders
BIOINFO_LIBS_DIR = $(LIBS_DIR)/bioinfo-libs
COMMON_LIBS_DIR = $(LIBS_DIR)/common-libs
CONTAINERS_DIR = $(COMMON_LIBS_DIR)/containers
COMMONS_DIR = $(COMMON_LIBS_DIR)/commons
BIOFORMATS_DIR = $(BIOINFO_LIBS_DIR)/bioformats
MATH_DIR = $(LIBS_DIR)/math

# -I (includes) and -L (libraries) paths
INCLUDES = -I $(SRC_DIR) -I $(LIBS_DIR) -I $(BIOINFO_LIBS_DIR) -I $(COMMON_LIBS_DIR) -I $(INC_DIR) -I /usr/include/libxml2 -I/usr/local/include
LIBS = -L/usr/lib/x86_64-linux-gnu -lcurl -Wl,-Bsymbolic-functions -lconfig -lcprops -fopenmp -lm -lxml2 -lgsl -lgslcblas -largtable2
LIBS_TEST = -lcheck

INCLUDES_STATIC = -I $(SRC_DIR) -I $(LIBS_DIR) -I $(BIOINFO_LIBS_DIR) -I $(COMMON_LIBS_DIR) -I $(INC_DIR) -I /usr/include/libxml2 -I/usr/local/include
LIBS_STATIC = -L$(LIBS_DIR) -L/usr/lib/x86_64-linux-gnu -lcurl -Wl,-Bsymbolic-functions -lconfig -lcprops -fopenmp -lm -lxml2 -lgsl -lgslcblas -largtable2


# Project dependencies
VCF_OBJS = $(BIOFORMATS_DIR)/vcf/vcf_*.o
GFF_OBJS = $(BIOFORMATS_DIR)/gff/gff_*.o
PED_OBJS = $(BIOFORMATS_DIR)/ped/ped_*.o
FAMILY_OBJS = $(BIOFORMATS_DIR)/family/family.o $(BIOFORMATS_DIR)/family/checks_family.o
REGION_TABLE_OBJS = $(BIOFORMATS_DIR)/features/region/region.o $(BIOFORMATS_DIR)/features/region/region_table.o $(BIOFORMATS_DIR)/features/region/region_table_utils.o
MISC_OBJS = $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/result.o \
	$(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(BIOFORMATS_DIR)/family.o $(MATH_DIR)/fisher.o


# Global targets
dist: compile-dependencies effect-dist gwas-dist vcf-dist
debug: compile-dependencies effect-debug gwas-debug vcf-debug


# hpg-var-effect
effect-dist: compile-dependencies
	cd effect && $(MAKE) effect-dist

effect-debug: compile-dependencies
	cd effect && $(MAKE) effect-debug


# hpg-var-gwas
gwas-dist: compile-dependencies
	cd gwas && $(MAKE) gwas-dist

gwas-debug: compile-dependencies
	cd gwas && $(MAKE) gwas-debug


# hpg-var-vcf
vcf-dist: compile-dependencies
	cd vcf-tools && $(MAKE) vcf-dist

vcf-debug: compile-dependencies
	cd vcf-tools && $(MAKE) vcf-debug


# Build dependencies
compile-dependencies:
	cd $(MATH_DIR) && make compile &&  \
	cd $(COMMONS_DIR) && make compile &&  \
	cd $(CONTAINERS_DIR) && make &&  \
	cd $(BIOFORMATS_DIR)/family && make &&  \
	cd $(BIOFORMATS_DIR)/features/region && make &&  \
	cd $(BIOFORMATS_DIR)/gff && make compile &&  \
	cd $(BIOFORMATS_DIR)/ped && make compile &&  \
	cd $(BIOFORMATS_DIR)/vcf && make compile

compile-dependencies-static:
	cd $(MATH_DIR) && make compile &&  \
	cd $(COMMONS_DIR) && make compile &&  \
	cd $(CONTAINERS_DIR) && make &&  \
	cd $(BIOFORMATS_DIR)/family && make &&  \
	cd $(BIOFORMATS_DIR)/features/region && make &&  \
	cd $(BIOFORMATS_DIR)/gff && make compile-static &&  \
	cd $(BIOFORMATS_DIR)/ped && make compile-static &&  \
	cd $(BIOFORMATS_DIR)/vcf && make compile-static

family.o:
	cd $(BIOFORMATS_DIR) && \
	$(CC) $(CFLAGS) -c -o $(BIOFORMATS_DIR)/$@ $(BIOFORMATS_DIR)/family.c $(INCLUDES) $(LIBS)
