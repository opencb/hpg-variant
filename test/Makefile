CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

ROOT = ~/appl/bioinfo-c

SRC_DIR = ..
BIOINFO_DATA_DIR = $(ROOT)/bioinfo-data
REGION_DIR = $(BIOINFO_DATA_DIR)/features/region
CONTAINERS_DIR = $(ROOT)/containers
COMMONS_DIR = $(ROOT)/commons

INCLUDES = -I $(SRC_DIR) -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(BIOINFO_DATA_DIR)/vcf -I $(REGION_DIR)
LIBS = -L/usr/lib/x86_64-linux-gnu -lcurl -Wl,-Bsymbolic-functions -lconfig -lcprops -fopenmp -lcheck -lm



VCF_FILES = $(BIOINFO_DATA_DIR)/vcf/vcf_file.c $(BIOINFO_DATA_DIR)/vcf/vcf_reader.c $(BIOINFO_DATA_DIR)/vcf/vcf_batch.c \
	$(BIOINFO_DATA_DIR)/vcf/vcf_read.c $(BIOINFO_DATA_DIR)/vcf/vcf_write.c \
	$(BIOINFO_DATA_DIR)/vcf/util.c $(BIOINFO_DATA_DIR)/vcf/vcf_filters.c
MISC_FILES = $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/http_utils.o $(CONTAINERS_DIR)/list.o $(REGION_DIR)/region.o $(CONTAINERS_DIR)/region_table.o 


EFFECT_FILES = $(SRC_DIR)/effect/main_effect.c $(SRC_DIR)/effect/effect_options_parsing.c $(SRC_DIR)/effect/effect_runner.c
VCF_TOOLS_FILES = $(SRC_DIR)/vcf-tools/main_vcf_tools.c $(SRC_DIR)/vcf-tools/vcf_tools_options_parsing.c $(SRC_DIR)/vcf-tools/vcf_tools_runner.c
HPG_VARIANT_FILES = $(SRC_DIR)/global_options.c $(EFFECT_FILES) $(VCF_TOOLS_FILES) $(VCF_FILES) $(MISC_FILES)

TEST_EFFECT_FILES = test_effect_runner.c $(HPG_VARIANT_FILES) 


all: list.o http_utils.o region.o region_table.o test_effect

test_effect: $(TEST_EFFECT_FILES)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(TEST_EFFECT_FILES) $(INCLUDES) $(LIBS) \
	&& ./test_effect

list.o:
	cd $(CONTAINERS_DIR) && make list.o

http_utils.o:
	cd $(COMMONS_DIR) && make http_utils.o

region.o:
	cd $(REGION_DIR) && make region.o

region_table.o:
	cd $(CONTAINERS_DIR) && make region_table.o

clean:
	rm test_effect

