#include "vcf_tools.h"

int vcf_handling(int argc, char *argv[], const char *config) {
    // Get sub-tool
    char *tool = argv[1];
    int exit_code = 0;
    
    // Parse tool args and run tool
    if (strcmp(tool, "filter") == 0) {
        exit_code = vcf_tool_filter(argc - 1, argv + 1, config);
        
    } else if (strcmp(tool, "merge") == 0) {
        exit_code = vcf_tool_merge(argc - 1, argv + 1, config);
        
    } else if (strcmp(tool, "split") == 0) {
        exit_code = vcf_tool_split(argc - 1, argv + 1, config);
        
    } else if (strcmp(tool, "stats") == 0) {
        exit_code = vcf_tool_stats(argc - 1, argv + 1, config);
        
    } else {
        fprintf(stderr, "The requested VCF handling tool does not exist! (%s)\n", tool);
        exit_code = NOT_IMPLEMENTED_TOOL;
    }
    
    if (exit_code > 0) {
        fprintf(stderr, "Tool %s terminated with failure (exit code = %d)\n", tool, exit_code);
    }
    return exit_code;
}

