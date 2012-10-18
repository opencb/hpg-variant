CC = gcc
CFLAGS = -std=c99 -O3 -march=native -D_XOPEN_SOURCE=600 -D_GNU_SOURCE
CFLAGS_DEBUG = -std=c99 -g -O0 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE

# Source code folders
INC_DIR = $(PWD)/include
LIBS_DIR = $(PWD)/libs
SRC_DIR = $(PWD)/src
BIN_DIR = $(PWD)/bin

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
	$(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(BIOFORMATS_DIR)/family.o $(MATH_DIR)/fisher.o

# Project source files
GLOBAL_FILES = $(SRC_DIR)/shared_options.c $(SRC_DIR)/hpg_variant_utils.c
EFFECT_FILES = $(SRC_DIR)/effect/*.c $(GLOBAL_FILES)
GWAS_FILES = $(SRC_DIR)/gwas/*.c $(SRC_DIR)/gwas/assoc/*.c $(SRC_DIR)/gwas/tdt/*.c $(GLOBAL_FILES)
# VARIANT_FILES = $(SRC_DIR)/effect/*.c $(SRC_DIR)/gwas/*.c
VCF_TOOLS_FILES = $(SRC_DIR)/vcf-tools/*.c $(SRC_DIR)/vcf-tools/filter/*.c $(SRC_DIR)/vcf-tools/merge/*.c $(SRC_DIR)/vcf-tools/split/*.c $(SRC_DIR)/vcf-tools/stats/*.c $(GLOBAL_FILES)
# ALL_FILES = $(SRC_DIR)/shared_options.c $(SRC_DIR)/hpg_variant_utils.c $(VARIANT_FILES) $(VCF_TOOLS_FILES)

# Project object files
HPG_VARIANT_OBJS = *.o
ALL_OBJS = $(HPG_VARIANT_OBJS) $(VCF_OBJS) $(GFF_OBJS) $(PED_OBJS) $(REGION_TABLE_OBJS) $(MISC_OBJS)

# Targets
all: clean clean-bin compile-dependencies effect-debug clean gwas-debug clean vcf-tools-debug

effect-deploy: clean compile-dependencies-static $(EFFECT_FILES)
	$(CC) $(CFLAGS) -c $(EFFECT_FILES) $(INCLUDES) $(LIBS_STATIC)
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/hpg-var-effect $(ALL_OBJS) $(INCLUDES_STATIC) $(LIBS_STATIC)

effect-debug: clean compile-dependencies $(EFFECT_FILES)
	$(CC) $(CFLAGS_DEBUG) -c $(EFFECT_FILES) $(INCLUDES) $(LIBS)
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS_DEBUG) -o $(BIN_DIR)/hpg-var-effect $(ALL_OBJS) $(INCLUDES) $(LIBS)

gwas-deploy: clean compile-dependencies-static $(GWAS_FILES)
	$(CC) $(CFLAGS) -c $(GWAS_FILES) $(INCLUDES) $(LIBS_STATIC)
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/hpg-var-gwas $(ALL_OBJS) $(INCLUDES_STATIC) $(LIBS_STATIC)

gwas-debug: clean compile-dependencies $(GWAS_FILES)
	$(CC) $(CFLAGS_DEBUG) -c $(GWAS_FILES) $(INCLUDES) $(LIBS)
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS_DEBUG) -o $(BIN_DIR)/hpg-var-gwas $(ALL_OBJS) $(INCLUDES) $(LIBS)

vcf-tools-deploy: clean compile-dependencies-static $(VCF_TOOLS_FILES)
	$(CC) $(CFLAGS) -c $(VCF_TOOLS_FILES) $(INCLUDES) $(LIBS_STATIC)
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/hpg-var-vcf $(ALL_OBJS) $(INCLUDES_STATIC) $(LIBS_STATIC)

vcf-tools-debug: clean compile-dependencies $(VCF_TOOLS_FILES)
	$(CC) $(CFLAGS_DEBUG) -c $(VCF_TOOLS_FILES) $(INCLUDES) $(LIBS)
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS_DEBUG) -o $(BIN_DIR)/hpg-var-vcf $(ALL_OBJS) $(INCLUDES) $(LIBS)


# deploy: compile-dependencies-static $(ALL_FILES)
# 	$(CC) $(CFLAGS) -c $(SRC_DIR)/main.c $(ALL_FILES) $(INCLUDES) $(LIBS_STATIC)
# 	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
# 	cp hpg-variant.cfg $(BIN_DIR)
# 	cp vcf-info-fields.cfg $(BIN_DIR)
# 	$(CC) $(CFLAGS) -o $(BIN_DIR)/hpg-variant $(ALL_OBJS) $(INCLUDES_STATIC) $(LIBS_STATIC)
# 	
# debug: compile-dependencies $(ALL_FILES)
# 	$(CC) $(CFLAGS_DEBUG) -c $(SRC_DIR)/main.c $(ALL_FILES) $(INCLUDES) $(LIBS)
# 	test -d $(BIN_DIR) || mkdir $(BIN_DIR)
# 	cp hpg-variant.cfg $(BIN_DIR)
# 	cp vcf-info-fields.cfg $(BIN_DIR)
# 	$(CC) $(CFLAGS_DEBUG) -o $(BIN_DIR)/hpg-variant $(ALL_OBJS) $(INCLUDES) $(LIBS)

testing: test/test_effect_runner.c test/test_tdt_runner.c $(ALL_OBJS)
	$(CC) $(CFLAGS_DEBUG) -o test/effect.test test/test_effect_runner.c $(ALL_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS_DEBUG) -o test/tdt.test test/test_tdt_runner.c $(ALL_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS_DEBUG) -o test/checks_family.test test/test_checks_family.c $(ALL_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	cp hpg-variant.cfg test
	cp vcf-info-fields.cfg test


compile-dependencies:
	make family.o && \
	cd $(MATH_DIR) && make compile &&  \
	cd $(COMMONS_DIR) && make compile &&  \
	cd $(CONTAINERS_DIR) && make &&  \
	cd $(BIOFORMATS_DIR)/features/region && make &&  \
	cd $(BIOFORMATS_DIR)/gff && make compile &&  \
	cd $(BIOFORMATS_DIR)/ped && make compile &&  \
	cd $(BIOFORMATS_DIR)/vcf && make compile

compile-dependencies-static:
	cd $(MATH_DIR) && make compile &&  \
	cd $(COMMONS_DIR) && make compile &&  \
	cd $(CONTAINERS_DIR) && make &&  \
	cd $(BIOFORMATS_DIR)/features/region && make &&  \
	cd $(BIOFORMATS_DIR)/gff && make compile-static &&  \
	cd $(BIOFORMATS_DIR)/ped && make compile-static &&  \
	cd $(BIOFORMATS_DIR)/vcf && make compile-static

family.o:
	cd $(BIOFORMATS_DIR) && \
	$(CC) $(CFLAGS) -c -o $(BIOFORMATS_DIR)/$@ $(BIOFORMATS_DIR)/family.c $(INCLUDES) $(LIBS)


clean:
	rm -f *.o
	rm -f $(COMMONS_DIR)/*.o
	rm -f $(CONTAINERS_DIR)/*.o
	rm -f $(BIOFORMATS_DIR)/vcf/*.o
	rm -f $(BIOFORMATS_DIR)/gff/*.o
	rm -f $(BIOFORMATS_DIR)/ped/*.o
	rm -f $(BIOFORMATS_DIR)/features/region/*.o

clean-bin:
	rm -rf bin
