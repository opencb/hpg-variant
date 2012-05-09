#include "gwas.h"
#include "tdt_runner.h"

int genomic_analysis(int argc, char *argv[]) {
    LOG_DEBUG_F("gwas called with %d args\n", argc);

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    global_options_data_t *global_options_data = init_global_options_data();
    gwas_options_data_t *options_data = init_options_data();


    /* ******************************
     *       Execution steps        *
     * ******************************/

    // Step 1: read options from configuration file
    int config_errors = read_global_configuration("hpg-variant.cfg", global_options_data);
    config_errors &= read_gwas_configuration("hpg-variant.cfg", options_data);
    LOG_INFO_F("Config read with errors = %d\n", config_errors);
    
    if (config_errors) {
        return CANT_READ_CONFIG_FILE;
    }
    
    // Step 2: parse command-line options
    parse_gwas_options(argc, argv, options_data, global_options_data);

    // Step 3: check that all options are set with valid values
    // Mandatory options that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_global_opts, check_gwas_opts;

    check_global_opts = verify_global_options(global_options_data);
    if (check_global_opts > 0)
    {
        return check_global_opts;
    }

    check_gwas_opts = verify_gwas_options(global_options_data, options_data);
    if (check_gwas_opts > 0)
    {
        return check_gwas_opts;
    }

    // Step 4: Perform the operations related to the selected GWAS sub-tool
//     switch (options_data->task) {
//         case TDT:
    if (options_data->task == TDT) {
            run_tdt_test(global_options_data, options_data);
    }
//             break;
//     }
    
    free_options_data(options_data);
    free_global_options_data(global_options_data);

    return 0;
}

gwas_options_data_t *init_options_data(void) {
    gwas_options_data_t *options_data = (gwas_options_data_t*) calloc (1, sizeof(gwas_options_data_t));
    
    options_data->task = NONE;
    options_data->max_batches = 10;
    options_data->batch_size = 2000;
    options_data->num_threads = 4;

    return options_data;
}

void free_options_data(gwas_options_data_t *options_data) {
    if (options_data->chain)    { cp_heap_destroy(options_data->chain); }
    free(options_data);
}
