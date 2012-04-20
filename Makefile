CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

# Source code folders
LIBS_ROOT = libs
CONTAINERS_DIR = $(LIBS_ROOT)/containers
COMMONS_DIR = $(LIBS_ROOT)/commons
BIOINFO_DATA_DIR = $(LIBS_ROOT)/bioformats
REGION_DIR = $(BIOINFO_DATA_DIR)/features/region

# Include and libs folders
INCLUDES = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR) -I $(BIOINFO_DATA_DIR)/vcf/ -I $(BIOINFO_DATA_DIR)/gff/  
LIBS = -L/usr/lib/x86_64-linux-gnu -lcurl -Wl,-Bsymbolic-functions -lconfig -lcprops -fopenmp -lm

# Source files dependencies
#VCF_FILES = $(BIOINFO_DATA_DIR)/vcf/vcf_file.c $(BIOINFO_DATA_DIR)/vcf/vcf_reader.c $(BIOINFO_DATA_DIR)/vcf/vcf_batch.c \
	$(BIOINFO_DATA_DIR)/vcf/vcf_read.c $(BIOINFO_DATA_DIR)/vcf/vcf_write.c \
	$(BIOINFO_DATA_DIR)/vcf/util.c $(BIOINFO_DATA_DIR)/vcf/vcf_filters.c
VCF_FILES = $(BIOINFO_DATA_DIR)/vcf/*.o
GFF_FILES = $(BIOINFO_DATA_DIR)/gff/*.o
REGION_TABLE_FILES = $(REGION_DIR)/region.o $(CONTAINERS_DIR)/region_table.o $(CONTAINERS_DIR)/region_table_utils.o
MISC_FILES = $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.c $(COMMONS_DIR)/string_utils.o $(CONTAINERS_DIR)/list.o

# Project source files
#VCF_TOOLS_FILES = main.c global_options.c vcf_tools_options_parsing.c vcf_tools_runner.c $(MISC_FILES)
VCF_TOOLS_FILES = main.c global_options.c vcf_tools_options_parsing.c vcf_tools_runner.c $(VCF_FILES) $(GFF_FILES) $(MISC_FILES) $(REGION_TABLE_FILES)


all: log.o list.o file_utils.o string_utils.o region.o region_table.o region_table_utils.o gff-reader vcf-reader vcf-tools

vcf-tools: $(VCF_TOOLS_FILES) gff-reader vcf-reader
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(VCF_TOOLS_FILES) $(INCLUDES) $(LIBS)
#	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ $(VCF_TOOLS_FILES) $(INCLUDES) $(LIBS)

gff-reader: 
	cd $(BIOINFO_DATA_DIR)/gff && make compile

vcf-reader: gff-reader
	cd $(BIOINFO_DATA_DIR)/vcf && make compile

file_utils.o:
	cd $(COMMONS_DIR) && make file_utils.o

string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o

log.o:
	cd $(COMMONS_DIR) && make log.o

list.o:
	cd $(CONTAINERS_DIR) && make list.o

region.o:
	cd $(REGION_DIR) && make region.o

region_table.o:
	cd $(CONTAINERS_DIR) && make region_table.o

region_table_utils.o:
	cd $(CONTAINERS_DIR) && make region_table_utils.o

clean:
	rm vcf-tools 

