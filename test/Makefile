CC = gcc
CFLAGS = -std=c99 -O3 -march=native -D_XOPEN_SOURCE=600 -D_GNU_SOURCE
CFLAGS_DEBUG = -std=c99 -g -O0 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE

# Source code folders
ROOT_DIR = $(PWD)/../
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
REGION_TABLE_OBJS = $(BIOFORMATS_DIR)/features/region/region.o $(BIOFORMATS_DIR)/features/region/region_table.o $(BIOFORMATS_DIR)/features/region/region_table_utils.o
MISC_OBJS = $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/result.o \
	$(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(BIOFORMATS_DIR)/family/family.o $(BIOFORMATS_DIR)/family/checks_family.o $(MATH_DIR)/stats/fisher.o
DEPEND_OBJS = $(VCF_OBJS) $(GFF_OBJS) $(PED_OBJS) $(REGION_TABLE_OBJS) $(MISC_OBJS)

# Project files
# EFFECT_OBJS = $(SRC_DIR)/effect/*.o $(SRC_DIR)/*.o
# GWAS_OBJS = $(SRC_DIR)/gwas/*.o $(SRC_DIR)/gwas/assoc/*.o $(SRC_DIR)/gwas/tdt/*.o $(SRC_DIR)/*.o
EFFECT_OBJS = $(SRC_DIR)/effect/auxiliary_files_writer.o $(SRC_DIR)/effect/effect_options_parsing.o $(SRC_DIR)/effect/effect_runner.o $(SRC_DIR)/*.o
GWAS_OBJS = $(SRC_DIR)/gwas/assoc/*.o $(SRC_DIR)/gwas/tdt/*.o $(SRC_DIR)/hpg_variant_utils.o $(SRC_DIR)/shared_options.o
VCF_TOOLS_OBJS = $(SRC_DIR)/vcf-tools/*.o $(SRC_DIR)/vcf-tools/filter/*.o $(SRC_DIR)/vcf-tools/merge/*.o $(SRC_DIR)/vcf-tools/split/*.o $(SRC_DIR)/vcf-tools/stats/*.o  $(SRC_DIR)/*.o


all: build

build: $(TEST_DIR)/test_checks_family.c $(TEST_DIR)/test_effect_runner.c $(TEST_DIR)/test_merge.c  $(TEST_DIR)/test_tdt_runner.c
	$(CC) $(CFLAGS_DEBUG) -o $(TEST_DIR)/checks_family.test $(TEST_DIR)/test_checks_family.c $(GWAS_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS_DEBUG) -o $(TEST_DIR)/effect.test $(TEST_DIR)/test_effect_runner.c $(EFFECT_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS_DEBUG) -o $(TEST_DIR)/merge.test $(TEST_DIR)/test_merge.c $(SRC_DIR)/vcf-tools/filter/*.o $(SRC_DIR)/vcf-tools/merge/*.o $(SRC_DIR)/vcf-tools/split/*.o $(SRC_DIR)/vcf-tools/stats/*.o $(SRC_DIR)/*.o $(DEPEND_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS_DEBUG) -o $(TEST_DIR)/tdt.test $(TEST_DIR)/test_tdt_runner.c $(GWAS_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
