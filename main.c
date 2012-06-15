#include "main.h"

int main(int argc, char *argv[])
{
	init_log_custom(2, 1, "hpg-variant.log");

    const char *config = find_configuration_file(argc, argv);
    
    mmap_vcf = 0;
	// Get tool
	char *tool = argv[1];
	int exit_code = 0;
	
	// Parse tool args and run tool
	if (strcmp(tool, "call") == 0) {
		fprintf(stderr, "%s tool not yet implemented\n", tool);
		exit_code = NOT_IMPLEMENTED_TOOL;
        
	} else if (strcmp(tool, "effect") == 0) {
		exit_code = effect(argc - 1, argv + 1, config);
        
	} else if (strcmp(tool, "functional-analysis") == 0) {
		fprintf(stderr, "%s tool not yet implemented\n", tool);
        exit_code = NOT_IMPLEMENTED_TOOL;
        
	} else if (strcmp(tool, "gwas") == 0) {
        exit_code = genomic_analysis(argc - 1, argv + 1, config);
        
	} else if (strcmp(tool, "pathway-analysis") == 0) {
		fprintf(stderr, "%s tool not yet implemented\n", tool);
        exit_code = NOT_IMPLEMENTED_TOOL;
        
	} else {
		fprintf(stderr, "The requested tool does not exist! (%s)\n", tool);
        exit_code = NOT_IMPLEMENTED_TOOL;
	}
	
	if (exit_code > 0)
	{
		fprintf(stderr, "Tool %s terminated with failure (exit code = %d)\n", tool, exit_code);
	}
	return exit_code;
}
