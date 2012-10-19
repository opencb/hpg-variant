SRC_DIR = $(PWD)/src
BIN_DIR = $(PWD)/bin
TEST_DIR = $(PWD)/test

all: clean build

# Deployment
build: build-effect build-gwas build-vcf

build-effect:
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	cd $(SRC_DIR) && $(MAKE) effect-deploy

build-gwas:
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	cd $(SRC_DIR) && $(MAKE) gwas-deploy

build-vcf:
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	cd $(SRC_DIR) && $(MAKE) vcf-deploy

# Debugging
debug: debug-effect debug-gwas debug-vcf

debug-effect:
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	cd $(SRC_DIR) && $(MAKE) effect-debug

debug-gwas:
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	cd $(SRC_DIR) && $(MAKE) gwas-debug

debug-vcf:
	mkdir -p $(BIN_DIR)
	cp hpg-variant.cfg $(BIN_DIR)
	cp vcf-info-fields.cfg $(BIN_DIR)
	cd $(SRC_DIR) && $(MAKE) vcf-debug

# Testing
test: $(TEST_DIR)/test_effect_runner.c $(TEST_DIR)/test_tdt_runner.c debug
	cp hpg-variant.cfg $(TEST_DIR)
	cp vcf-info-fields.cfg $(TEST_DIR)
	cd $(TEST_DIR) && $(MAKE)

# Clean-up
clean:
	rm -rf *.o
	rm -rf $(BIN_DIR)
