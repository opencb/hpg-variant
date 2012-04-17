#include "vcf_tools.h"


int read_vcf_tools_configuration(const char *filename, vcf_tools_options_data_t *options_data) {
    config_t *config = (config_t*) malloc (sizeof(config_t));
    
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE)
    {
        fprintf(stderr, "config file error: %s\n", config_error_text(config));
        return ret_code;
    }
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "vcf-tools.num-threads", &(options_data->num_threads));
    if (ret_code == CONFIG_FALSE)
    {
        fprintf(stderr, "num-threads config error: %s\n", config_error_text(config));
        return ret_code;
    }
    dprintf("num-threads = %ld\n", options_data->num_threads);
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "vcf-tools.max-batches", &(options_data->max_batches));
    if (ret_code == CONFIG_FALSE)
    {
        fprintf(stderr, "max-batches config error: %s\n", config_error_text(config));
        return ret_code;
    }
    dprintf("max-batches = %ld\n", options_data->max_batches);
    
    // Read number of threads to perform the operations
    ret_code = config_lookup_int(config, "vcf-tools.batch-size", &(options_data->batch_size));
    if (ret_code == CONFIG_FALSE)
    {
        fprintf(stderr, "batch-size config error: %s\n", config_error_text(config));
        return ret_code;
    }
    dprintf("batch-size = %ld\n", options_data->batch_size);
    
    return ret_code;
}


void parse_vcf_tools_options(int argc, char *argv[], vcf_tools_options_data_t *options_data, global_options_data_t *global_options_data) {
	const struct option *options = merge_options(vcf_tools_options, NUM_VCF_TOOLS_OPTIONS);

	char *tmp_string_field;
	int tmp_int_field;
    filter_t *filter;
//     char filter_input[];
    
	int c;
	// Last option read, for when the global options parser is invoked
	int previous_opt_index = optind;
	
	int debug = 1;
	while ((c = getopt_long (argc, argv, ":A:O:f:r:s:", options, &optind)) != -1) {
		dprintf("<main> c = %c, opt_idx = %d\n", c, optind);
		switch (c) {
			case 'A':
			case 'O':
				optind = parse_global_options(argc, argv, global_options_data, previous_opt_index);
				break;
            case 'f':
                tmp_string_field = (char*) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(tmp_string_field, optarg);
//                 filter_input = optarg;
                filter = create_region_filter(tmp_string_field, 1);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                dprintf("regions file = %s\n", optarg);
                break;
			case 'r':
                tmp_string_field = (char*) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(tmp_string_field, optarg);
//                 filter_input = optarg;
                filter = create_region_filter(tmp_string_field, 0);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                dprintf("regions = %s\n", optarg);
				break;
			case 's':
                filter = create_snp_filter(optarg);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
				dprintf("snp filter to %s SNPs\n", (optarg == NULL)? "include" : optarg);
				break;
			case '?':
			default:
				dprintf("Option unknown\n");
				break;
		}
		
		previous_opt_index = optind;
	}
}


int verify_vcf_tools_options(global_options_data_t *global_options_data, vcf_tools_options_data_t *options_data)
{
    // Check whether the input VCF file is defined
    if (global_options_data->vcf_filename == NULL || strlen(global_options_data->vcf_filename) == 0) {
        fprintf(stderr, "Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
	
	return 0;
}

