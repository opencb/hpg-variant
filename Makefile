CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

# Source code folders
INC_DIR = $(PWD)/include
LIBS_DIR = $(PWD)/libs
SRC_DIR = $(PWD)/src
BIOINFO_LIBS_DIR = $(LIBS_DIR)/bioinfo-libs
COMMON_LIBS_DIR = $(LIBS_DIR)/common-libs

CONTAINERS_DIR = $(COMMON_LIBS_DIR)/containers
COMMONS_DIR = $(COMMON_LIBS_DIR)/commons
BIOFORMATS_DIR = $(BIOINFO_LIBS_DIR)/bioformats

# Include and libs folders
INCLUDES = -I $(SRC_DIR) -I $(LIBS_DIR) -I $(BIOINFO_LIBS_DIR) -I $(COMMON_LIBS_DIR)  -I $(INC_DIR)
LIBS = -L/usr/lib/x86_64-linux-gnu -lconfig -lcprops -fopenmp -lm -lcurl -Wl,-Bsymbolic-functions -largtable2
LIBS_TEST = -lcheck

INCLUDES_STATIC = -I $(SRC_DIR) -I $(LIBS_DIR) -I $(BIOINFO_LIBS_DIR) -I $(COMMON_LIBS_DIR)  -I $(INC_DIR)
LIBS_STATIC = -L$(LIBS_DIR) -L/usr/lib/x86_64-linux-gnu -lconfig -lcprops -fopenmp -lm -lcurl -Wl,-Bsymbolic-functions -largtable2

# Source files dependencies
VCF_OBJS = $(BIOFORMATS_DIR)/vcf/vcf_*.o
GFF_OBJS = $(BIOFORMATS_DIR)/gff/gff_*.o
MISC_OBJS = $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/string_utils.o \
	$(CONTAINERS_DIR)/list.o $(CONTAINERS_DIR)/array_list.o \
	$(BIOFORMATS_DIR)/features/region/region.o $(BIOFORMATS_DIR)/features/region/region_table.o $(BIOFORMATS_DIR)/features/region/region_table_utils.o

# Project source files
VCF_TOOLS_FILES = $(SRC_DIR)/global_options.c $(SRC_DIR)/hpg_vcf_tools_utils.c $(SRC_DIR)/filter/*.c $(SRC_DIR)/merge/merge.c $(SRC_DIR)/split/*.c $(SRC_DIR)/stats/*.c

# Project object files
VCF_TOOLS_OBJS = *.o
ALL_OBJS = $(VCF_TOOLS_OBJS) $(VCF_OBJS) $(GFF_OBJS) $(MISC_OBJS)

# Targets
all: compile-dependencies hpg-vcf

hpg-vcf: compile-dependencies $(VCF_TOOLS_FILES)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -c $(SRC_DIR)/main.c $(VCF_TOOLS_FILES) $(INCLUDES) $(LIBS)
	test -d bin || mkdir bin
	cp hpg-vcf-tools.cfg bin
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o bin/$@ $(ALL_OBJS) $(INCLUDES) $(LIBS)

deploy: compile-dependencies-static $(VCF_TOOLS_FILES)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(SRC_DIR)/main.c $(VCF_TOOLS_FILES) $(INCLUDES_STATIC) $(LIBS_STATIC)
	test -d bin || mkdir bin
	cp hpg-vcf-tools.cfg bin
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o bin/hpg-vcf $(ALL_OBJS) $(INCLUDES_STATIC) $(LIBS_STATIC)

test: hpg-vcf
# 	TEST_OBJS = `echo $(ALL_OBJS) | sed s/main.o//g`
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o test/merge.test test/test_merge.c $(VCF_TOOLS_FILES) $(VCF_OBJS) $(GFF_OBJS) $(MISC_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)
# 	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o test/merge.test test/test_merge.o $(TEST_OBJS) $(INCLUDES) $(LIBS) $(LIBS_TEST)

compile-dependencies:
	make family.o && \
	cd $(COMMONS_DIR) && make compile &&  \
	cd $(CONTAINERS_DIR) && make compile &&  \
	cd $(BIOFORMATS_DIR)/features/region && make &&  \
	cd $(BIOFORMATS_DIR)/gff && make compile &&  \
	cd $(BIOFORMATS_DIR)/vcf && make compile

compile-dependencies-static:
	cd $(COMMONS_DIR) && make compile &&  \
	cd $(CONTAINERS_DIR) && make &&  \
	cd $(BIOFORMATS_DIR)/features/region && make &&  \
	cd $(BIOFORMATS_DIR)/gff && make compile-static &&  \
	cd $(BIOFORMATS_DIR)/vcf && make compile-static && \
	cd $(LIBS_DIR)/sunrisedd && make all

family.o:
	cd $(BIOFORMATS_DIR) && \
	$(CC) $(CFLAGS) -c -o $(BIOFORMATS_DIR)/$@ $(BIOFORMATS_DIR)/family.c $(INCLUDES) $(LIBS)

clean:
	rm -f *.o
	rm -f $(CONTAINERS_DIR)/*.o
	rm -f $(COMMONS_DIR)/*.o
	rm -f $(BIOFORMATS_DIR)/vcf/*.o
	rm -f $(BIOFORMATS_DIR)/gff/*.o
	rm -f $(BIOFORMATS_DIR)/ped/*.o
	rm -f $(BIOFORMATS_DIR)/features/region/*.o
	rm -f $(LIBS_DIR)/sunrisedd/*.o
	rm -f $(LIBS_DIR)/sunrisedd/sunrisedd.a
	rm -rf bin
