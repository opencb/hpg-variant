#include "vcf_tools.h"
#include "vcf_tools_runner.h"

int main(int argc, char *argv[])
{
    init_log_custom(2, 1, "hpg-vcf-tools.log");

    /* ******************************
        * 	Modifiable options	*
        * ******************************/

    global_options_data_t *global_options_data = init_global_options_data();
    vcf_tools_options_data_t *options_data = init_options_data();


    /* ******************************
        * 	 Execution steps	*
        * ******************************/

    // Step 1: read options from configuration file
    int config_read = read_global_configuration("hpg-vcf-tools.cfg", global_options_data);
    config_read &= read_vcf_tools_configuration("hpg-vcf-tools.cfg", options_data);
    LOG_DEBUG_F("Config read successfully = %d\n", config_read);

    // Step 2: parse command-line options
    parse_vcf_tools_options(argc, argv, options_data, global_options_data);

    // Step 3: check that all options are set with valid values
    // Mandatory that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_global_opts, check_vcf_tools_opts;

    check_global_opts = verify_global_options(global_options_data);
    if (check_global_opts > 0)
    {
        return check_global_opts;
    }

    check_vcf_tools_opts = verify_vcf_tools_options(global_options_data, options_data);
    if (check_vcf_tools_opts > 0)
    {
        return check_vcf_tools_opts;
    }

    // Step 4: Perform the requested tasks
    execute_vcf_tools(global_options_data, options_data);

    free_options_data(options_data);
    free_global_options_data(global_options_data);

    return 0;
}

vcf_tools_options_data_t *init_options_data()
{
    vcf_tools_options_data_t *options_data = (vcf_tools_options_data_t*) malloc (sizeof(vcf_tools_options_data_t));

    options_data->num_threads = 1;
    options_data->max_batches = 10;
    options_data->batch_size = 20000;

    return options_data;
}

void free_options_data(vcf_tools_options_data_t *options_data)
{
    free(options_data);
}

