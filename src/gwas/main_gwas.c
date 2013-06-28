/*
 * Copyright (c) 2012-2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
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

#include "main_gwas.h"

int main(int argc, char *argv[]) {
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        printf("Usage: %s < assoc | epi | tdt > < tool-options >\nFor more information about a certain tool, type %s tool-name --help\n", argv[0], argv[0]);
        return 0;
    }
    
    init_log_custom(2, 1, "hpg-var-gwas.log", "w");

    array_list_t *config_search_paths = NULL;
    const char *config = NULL;
    int config_len;

#ifdef _USE_MPI
    int mpi_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (mpi_rank == 0) {
#endif
        
    config_search_paths = get_configuration_search_paths(argc, argv);
    config = retrieve_config_file("hpg-variant.conf", config_search_paths);
    config_len = strlen(config);
        
#ifdef _USE_MPI
    }
    MPI_Bcast(&config_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mpi_rank != 0) {
        config = malloc (config_len * sizeof(char));
    }

    MPI_Bcast(config, config_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    LOG_INFO_F("Rank %d got config file: %.*s\n", mpi_rank, config_len, config);
#endif

    // Get sub-tool
    char *tool = argv[1];
    int exit_code = 0;
    
    // Parse tool args and run tool
    if (strcmp(tool, "assoc") == 0) {
        exit_code = association(argc - 1, argv + 1, config);
        
    } else if (strcmp(tool, "epi") == 0) {
        exit_code = epistasis(argc - 1, argv + 1, config);
        
    } else if (strcmp(tool, "tdt") == 0) {
        exit_code = tdt(argc - 1, argv + 1, config);
 
    } else {
        fprintf(stderr, "The requested genome-wide analysis tool does not exist! (%s)\n", tool);
        exit_code = NOT_IMPLEMENTED_TOOL;
    }
    
    if (exit_code > 0) {
        fprintf(stderr, "Tool %s terminated with failure (exit code = %d)\n", tool, exit_code);
    }
  
#ifdef _USE_MPI
    if (mpi_rank == 0) { 
        if (config_search_paths) { array_list_free(config_search_paths, free); }
    }
    
    MPI_Finalize();
#else
    array_list_free(config_search_paths, free);
    
#endif
    
    stop_log();
    return exit_code;
}
