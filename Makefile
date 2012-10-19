CC = gcc
CFLAGS = -std=c99 -O3 -march=native -D_XOPEN_SOURCE=600 -D_GNU_SOURCE
CFLAGS_DEBUG = -std=c99 -g -O0 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE

# Source code folders
INC_DIR = $(PWD)/include
LIBS_DIR = $(PWD)/libs
SRC_DIR = $(PWD)/src
BIN_DIR = $(PWD)/bin
TEST_DIR = $(PWD)/test

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
VCF_TOOLS_FILES = $(SRC_DIR)/vcf-tools/*.c $(SRC_DIR)/vcf-tools/filter/*.c $(SRC_DIR)/vcf-tools/merge/*.c $(SRC_DIR)/vcf-tools/split/*.c $(SRC_DIR)/vcf-tools/stats/*.c $(GLOBAL_FILES)

# Project object files
HPG_VARIANT_OBJS = *.o
EFFECT_OBJS = $(SRC_DIR)/effect/*.o $(SRC_DIR)/*.o
GWAS_OBJS = $(SRC_DIR)/gwas/*.o $(SRC_DIR)/gwas/assoc/*.o $(SRC_DIR)/gwas/tdt/*.o $(SRC_DIR)/*.o
VCF_TOOLS_OBJS = $(SRC_DIR)/vcf-tools/*.o $(SRC_DIR)/vcf-tools/filter/*.o $(SRC_DIR)/vcf-tools/merge/*.o $(SRC_DIR)/vcf-tools/split/*.o $(SRC_DIR)/vcf-tools/stats/*.o $(SRC_DIR)/*.o
DEPEND_OBJS = $(VCF_OBJS) $(GFF_OBJS) $(PED_OBJS) $(REGION_TABLE_OBJS) $(MISC_OBJS)
ALL_OBJS = $(HPG_VARIANT_OBJS) $(VCF_OBJS) $(GFF_OBJS) $(PED_OBJS) $(REGION_TABLE_OBJS) $(MISC_OBJS)

# Global target
all: clean compile-dependencies effect-debug gwas-debug vcf-tools-debug

# hpg-var-effect targets
effect-deploy: compile-dependencies $(SRC_DIR)/effect/%.o
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/hpg-var-effect $(EFFECT_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS)

effect-obj-deploy: $(EFFECT_FILES)
	for F in $(EFFECT_FILES); do \
		$(CC) $(CFLAGS) -c $$F -o $$F.o $(INCLUDES) $(LIBS); \
	done;

effect-debug: compile-dependencies effect-obj-debug
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS_DEBUG) -o $(BIN_DIR)/hpg-var-effect $(EFFECT_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS)

effect-obj-debug: $(EFFECT_FILES)
	for F in $(EFFECT_FILES); do \
		$(CC) $(CFLAGS_DEBUG) -c $$F -o $$F.o $(INCLUDES) $(LIBS); \
	done;

# hpg-var-gwas targets
gwas-deploy: compile-dependencies $(SRC_DIR)/gwas/%.o
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/hpg-var-gwas $(GWAS_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS)

gwas-obj-deploy: $(GWAS_FILES)
	for F in $(GWAS_FILES); do \
		$(CC) $(CFLAGS) -c $$F -o $$F.o $(INCLUDES) $(LIBS); \
	done;

gwas-debug: compile-dependencies gwas-obj-debug
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS_DEBUG) -o $(BIN_DIR)/hpg-var-gwas $(GWAS_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS)

gwas-obj-debug: $(GWAS_FILES)
	for F in $(GWAS_FILES); do \
		$(CC) $(CFLAGS_DEBUG) -c $$F -o $$F.o $(INCLUDES) $(LIBS); \
	done;

# hpg-var-vcf targets
vcf-deploy: compile-dependencies $(SRC_DIR)/vcf-tools/%.o
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS) -o $(BIN_DIR)/hpg-var-vcf $(VCF_TOOLS_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS)

vcf-obj-deploy: $(VCF_TOOLS_FILES)
	for F in $(VCF_TOOLS_FILES); do \
		$(CC) $(CFLAGS) -c $$F -o $$F.o $(INCLUDES) $(LIBS); \
	done;

vcf-debug: compile-dependencies vcf-obj-debug
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	$(CC) $(CFLAGS_DEBUG) -o $(BIN_DIR)/hpg-var-vcf $(VCF_TOOLS_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS)

vcf-obj-debug: $(VCF_TOOLS_FILES)
	for F in $(VCF_TOOLS_FILES); do \
		$(CC) $(CFLAGS_DEBUG) -c $$F -o $$F.o $(INCLUDES) $(LIBS); \
	done;

testing: $(TEST_DIR)/test_effect_runner.c $(TEST_DIR)/test_tdt_runner.c effect-debug gwas-debug vcf-debug
	rm $(SRC_DIR)/effect/main_effect.c.o $(SRC_DIR)/gwas/main_gwas.c.o $(SRC_DIR)/vcf-tools/main_vcf_tools.c.o
	$(CC) $(CFLAGS_DEBUG) -o $(TEST_DIR)/checks_family.test $(TEST_DIR)/test_checks_family.c $(GWAS_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS_DEBUG) -o $(TEST_DIR)/effect.test $(TEST_DIR)/test_effect_runner.c $(EFFECT_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS_DEBUG) -o $(TEST_DIR)/merge.test $(TEST_DIR)/test_merge.c $(SRC_DIR)/vcf-tools/filter/*.o $(SRC_DIR)/vcf-tools/merge/*.o $(SRC_DIR)/vcf-tools/split/*.o $(SRC_DIR)/vcf-tools/stats/*.o $(SRC_DIR)/*.o $(DEPEND_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	$(CC) $(CFLAGS_DEBUG) -o $(TEST_DIR)/tdt.test $(TEST_DIR)/test_tdt_runner.c $(GWAS_OBJS) $(DEPEND_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
	cp hpg-variant.cfg $(TEST_DIR)
	cp vcf-info-fields.cfg $(TEST_DIR)


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
	rm -rf *.o
	rm -rf bin
