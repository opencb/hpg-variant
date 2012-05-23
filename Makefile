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
INCLUDES = -I . -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR) -I $(BIOINFO_DATA_DIR)/vcf/ -I $(BIOINFO_DATA_DIR)/gff/  
LIBS = -L/usr/lib/x86_64-linux-gnu -lconfig -lcprops -fopenmp -lm
LIBS_TEST = -lcheck

# Source files dependencies
VCF_FILES = $(BIOINFO_DATA_DIR)/vcf/*.o
GFF_FILES = $(BIOINFO_DATA_DIR)/gff/*.o
REGION_TABLE_FILES = $(REGION_DIR)/region.o $(CONTAINERS_DIR)/region_table.o $(CONTAINERS_DIR)/region_table_utils.o
MISC_FILES = $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.c $(COMMONS_DIR)/string_utils.o $(CONTAINERS_DIR)/list.o

# Project source files
FILTER_FILES = filter/main_filter.c filter/filter_runner.c filter/filter_options_parsing.c
SPLIT_FILES = split/main_split.c split/split.c split/split_runner.c split/split_options_parsing.c
STATS_FILES = stats/main_stats.c stats/stats.c stats/stats_runner.c stats/stats_options_parsing.c
VCF_TOOLS_FILES = global_options.c hpg_vcf_tools_utils.c $(FILTER_FILES) $(SPLIT_FILES) $(STATS_FILES) $(VCF_FILES) $(GFF_FILES) $(MISC_FILES) $(REGION_TABLE_FILES)


# Targets
all: compile-dependencies hpg-vcf-tools

hpg-vcf-tools: compile-dependencies $(VCF_TOOLS_FILES) 
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ main.c $(VCF_TOOLS_FILES) $(INCLUDES) $(LIBS)
#	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ main.c $(VCF_TOOLS_FILES) $(INCLUDES) $(LIBS)

testing: test/test_stats.c $(VCF_TOOLS_FILES)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o test/stats.test test/test_stats.c $(VCF_TOOLS_FILES) $(INCLUDES) $(LIBS) $(LIBS_TEST)
# 	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o test/stats.test test/test_stats.c $(VCF_TOOLS_FILES) $(INCLUDES) $(LIBS) $(LIBS_TEST)

compile-dependencies:
	cd $(COMMONS_DIR) && make file_utils.o string_utils.o log.o &&  \
	cd $(REGION_DIR) && make region.o &&  \
	cd $(CONTAINERS_DIR) && make list.o region_table.o region_table_utils.o &&  \
	cd $(BIOINFO_DATA_DIR)/gff && make compile &&  \
	cd $(BIOINFO_DATA_DIR)/vcf && make compile

clean:
	rm -f *.o
	rm -f $(CONTAINERS_DIR)/*.o
	rm -f $(COMMONS_DIR)/*.o
	rm -f $(REGION_DIR)/*.o
	rm hpg-vcf-tools 

