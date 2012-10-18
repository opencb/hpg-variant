/*
 * Copyright (c) 2012 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * hpg-variant is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-variant. If not, see <http://www.gnu.org/licenses/>.
 */

#include "main.h"

int main(int argc, char *argv[]) {
    if (argc == 1 || !strcmp(argv[1], "--help")) {
        printf("Usage: hpg-variant < effect | gwas | vcf > < tool-options >\nFor more information about a certain tool, type hpg-variant tool-name --help\n");
        return 0;
    }
    
    init_log_custom(2, 1, "hpg-variant.log");

    const char *config = find_configuration_file(argc, argv);
    
    // Get tool
    char *tool = argv[1];
    int exit_code = 0;

    // Parse tool args and run tool
    if (!strcmp(tool, "call")) {
        fprintf(stderr, "%s tool not yet implemented\n", tool);
        exit_code = NOT_IMPLEMENTED_TOOL;
        
    } else if (!strcmp(tool, "effect")) {
        exit_code = effect(argc - 1, argv + 1, config);
        
    } else if (!strcmp(tool, "functional-analysis")) {
        fprintf(stderr, "%s tool not yet implemented\n", tool);
        exit_code = NOT_IMPLEMENTED_TOOL;
        
    } else if (!strcmp(tool, "gwas")) {
        exit_code = genomic_analysis(argc - 1, argv + 1, config);
        
    } else if (!strcmp(tool, "pathway-analysis")) {
        fprintf(stderr, "%s tool not yet implemented\n", tool);
        exit_code = NOT_IMPLEMENTED_TOOL;
        
    } else if (!strcmp(tool, "vcf")) {
        exit_code = vcf_handling(argc - 1, argv + 1, config);
        
    } else {
        fprintf(stderr, "The requested tool does not exist! (%s)\n", tool);
        exit_code = NOT_IMPLEMENTED_TOOL;
    }

    if (exit_code > 0) {
        fprintf(stderr, "Tool %s terminated with failure (exit code = %d)\n", tool, exit_code);
    }
    return exit_code;
}
