CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

# Source code folders
INC_DIR = $(PWD)/include
LIBS_DIR = $(PWD)/libs
CONTAINERS_DIR = $(LIBS_DIR)/containers
COMMONS_DIR = $(LIBS_DIR)/commons
BIOFORMATS_DIR = $(LIBS_DIR)/bioformats

# Include and libs folders
INCLUDES = -I . -I $(LIBS_DIR) -I $(INC_DIR)
LIBS = -L/usr/lib/x86_64-linux-gnu -lconfig -lcprops -fopenmp -lm -lcurl -Wl,-Bsymbolic-functions -largtable2
LIBS_TEST = -lcheck

INCLUDES_STATIC = -I . -I $(LIBS_DIR) -I $(INC_DIR)
LIBS_STATIC = -Llibs -L/usr/lib/x86_64-linux-gnu -lconfig -lcprops -fopenmp -lm -lcurl -Wl,-Bsymbolic-functions -largtable2

# Source files dependencies
VCF_OBJS = $(BIOFORMATS_DIR)/vcf/vcf_*.o
GFF_OBJS = $(BIOFORMATS_DIR)/gff/gff_*.o
MISC_OBJS = $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/string_utils.o \
	$(CONTAINERS_DIR)/list.o $(CONTAINERS_DIR)/array_list.o \
	$(CONTAINERS_DIR)/region_table.o $(CONTAINERS_DIR)/region_table_utils.o $(BIOFORMATS_DIR)/features/region/region.o

# Project source files
FILTER_FILES = filter/main_filter.c filter/filter_runner.c filter/filter_options_parsing.c
SPLIT_FILES = split/main_split.c split/split.c split/split_runner.c split/split_options_parsing.c
STATS_FILES = stats/main_stats.c stats/stats_runner.c stats/stats_options_parsing.c
VCF_TOOLS_FILES = global_options.c hpg_vcf_tools_utils.c $(FILTER_FILES) $(SPLIT_FILES) $(STATS_FILES)

# Project object files
FILTER_OBJS = main_filter.o filter_runner.o filter_options_parsing.o
SPLIT_OBJS = main_split.o split.o split_runner.o split_options_parsing.o
STATS_OBJS = main_stats.o stats_runner.o stats_options_parsing.o
VCF_TOOLS_OBJS = global_options.o hpg_vcf_tools_utils.o $(FILTER_OBJS) $(SPLIT_OBJS) $(STATS_OBJS) $(VCF_OBJS) $(GFF_OBJS) $(MISC_OBJS)

# Targets
all: compile-dependencies hpg-vcf

hpg-vcf: compile-dependencies $(VCF_TOOLS_FILES)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -c main.c $(VCF_TOOLS_FILES) $(INCLUDES) $(LIBS)
	test -d bin || mkdir bin
	cp hpg-vcf-tools.cfg bin
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o bin/$@ main.o $(VCF_TOOLS_OBJS) $(INCLUDES) $(LIBS)

deploy: compile-dependencies-static $(VCF_TOOLS_FILES)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -c main.c $(VCF_TOOLS_FILES) $(INCLUDES_STATIC) $(LIBS_STATIC)
	test -d bin || mkdir bin
	cp hpg-vcf-tools.cfg bin
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o bin/hpg-vcf main.o $(VCF_TOOLS_OBJS) $(INCLUDES_STATIC) $(LIBS_STATIC)

compile-dependencies:
	cd $(COMMONS_DIR) && make compile &&  \
	cd $(CONTAINERS_DIR) && make compile &&  \
	cd $(BIOFORMATS_DIR)/features/region && make region.o &&  \
	cd $(BIOFORMATS_DIR)/gff && make compile &&  \
	cd $(BIOFORMATS_DIR)/vcf && make compile

compile-dependencies-static:
	cd $(COMMONS_DIR) && make compile &&  \
	cd $(CONTAINERS_DIR) && make compile-static &&  \
	cd $(BIOFORMATS_DIR)/features/region && make region.o &&  \
	cd $(BIOFORMATS_DIR)/gff && make compile-static &&  \
	cd $(BIOFORMATS_DIR)/vcf && make compile-static

clean:
	rm -f *.o
	rm -f $(CONTAINERS_DIR)/*.o
	rm -f $(COMMONS_DIR)/*.o
	rm -f $(BIOFORMATS_DIR)/features/region/*.o
	rm -rf bin
