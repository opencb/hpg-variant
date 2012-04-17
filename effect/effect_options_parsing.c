#include "effect.h"


int read_effect_configuration(const char *filename, effect_options_data_t *options_data)
{
    if (filename == NULL || options_data == NULL) {
        return -1;
    }
    
	config_t *config = (config_t*) malloc (sizeof(config_t));
	int ret_code = config_read_file(config, filename);
	if (ret_code == CONFIG_FALSE)
	{
		fprintf(stderr, "config file error: %s\n", config_error_text(config));
		return ret_code;
	}
	
    char *tmp_url, *tmp_version, *tmp_species;
    
	// Read URL
    ret_code = config_lookup_string(config, "effect.url", &tmp_url);
	if (ret_code == CONFIG_FALSE)
	{
		fprintf(stderr, "url config error: %s\n", config_error_text(config));
		return ret_code;
	} else {
        free(options_data->host_url);
        options_data->host_url = (char*) calloc (strlen(tmp_url)+1, sizeof(char));
        strncat(options_data->host_url, tmp_url, strlen(tmp_url));
        dprintf("URL = %s (%zu chars)\n", options_data->host_url, strlen(options_data->host_url));
    }
	
	// Read version
    ret_code = config_lookup_string(config, "effect.version", &tmp_version);
	if (ret_code == CONFIG_FALSE)
	{
		fprintf(stderr, "version config error: %s\n", config_error_text(config));
		return ret_code;
	} else {
        free(options_data->version);
        options_data->version = (char*) calloc (strlen(tmp_version)+1, sizeof(char));
        strncat(options_data->version, tmp_version, strlen(tmp_version));
        dprintf("version = %s (%zu chars)\n", options_data->version, strlen(options_data->version));
    }
	
	// Read species
	ret_code = config_lookup_string(config, "effect.species", &tmp_species);
	if (ret_code == CONFIG_FALSE)
	{
		fprintf(stderr, "species config error: %s\n", config_error_text(config));
		return ret_code;
	} else {
        free(options_data->species);
        options_data->species = (char*) calloc (strlen(tmp_species)+1, sizeof(char));
        strncat(options_data->species, tmp_species, strlen(tmp_species));
        dprintf("species = %s (%zu chars)\n", options_data->species, strlen(options_data->species));
    }
	
	// Read number of threads that will make request to the web service
	ret_code = config_lookup_int(config, "effect.num-threads", &(options_data->num_threads));
	if (ret_code == CONFIG_FALSE)
	{
		fprintf(stderr, "num-threads config error: %s\n", config_error_text(config));
		return ret_code;
	}
	dprintf("num-threads = %ld\n", options_data->num_threads);
	
    // Read maximum number of batches that can be stored at certain moment
    ret_code = config_lookup_int(config, "effect.max-batches", &(options_data->max_batches));
    if (ret_code == CONFIG_FALSE)
    {
        fprintf(stderr, "max-batches config error: %s\n", config_error_text(config));
        return ret_code;
    }
    dprintf("max-batches = %ld\n", options_data->max_batches);
    
    // Read size of every batch read
    ret_code = config_lookup_int(config, "effect.batch-size", &(options_data->batch_size));
    if (ret_code == CONFIG_FALSE)
    {
        fprintf(stderr, "batch-size config error: %s\n", config_error_text(config));
        return ret_code;
    }
    dprintf("batch-size = %ld\n", options_data->batch_size);
    
	// Read number of variants per request to the web service
	ret_code = config_lookup_int(config, "effect.variants-per-request", &(options_data->variants_per_request));
	if (ret_code == CONFIG_FALSE)
	{
		fprintf(stderr, "variants-per-request config error: %s\n", config_error_text(config));
		return ret_code;
	}
	dprintf("variants-per-request = %ld\n", options_data->variants_per_request);
	
    config_destroy(config);
    free(config);
    
	return ret_code;
}


void parse_effect_options(int argc, char *argv[], effect_options_data_t *options_data, global_options_data_t *global_options_data)
{
	struct option *options = merge_options(effect_options, NUM_EFFECT_OPTIONS);

    char *tmp_string_field;
    int tmp_int_field;
    filter_t *filter;

	int c;
	// Last option read, for when the global options parser is invoked
	int previous_opt_index = optind;
	
	int debug = 1;
	while ((c = getopt_long (argc, argv, ":A:O:e:n:s:u:", options, &optind)) != -1)
	{
		dprintf("<main> c = %c, opt_idx = %d\n", c, optind);
		switch (c)
		{
			case 'A':
			case 'O':
				optind = parse_global_options(argc, argv, global_options_data, previous_opt_index);
				break;
			case 'e':
				options_data->excludes = (char*) malloc((strlen(optarg)+1) * sizeof(char));
				strcpy(options_data->excludes, optarg);
				dprintf("exclude = %s\n", options_data->excludes);
				break;
            case 'f':
                tmp_string_field = (char*) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(tmp_string_field, optarg);
                filter = create_region_filter(tmp_string_field, 1);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                dprintf("regions file = %s\n", optarg);
                break;
			case 'n':
				tmp_int_field = atoi(optarg);
				if (!tmp_int_field)
				{
					fprintf(stderr, "The requested number of threads is not valid (%s)\n", optarg);
				} else
				{
					options_data->num_threads = tmp_int_field;
				}
				dprintf("num-threads = %ld\n", options_data->num_threads);
				break;
			case 'p':
				tmp_int_field = atoi(optarg);
				if (!tmp_int_field)
				{
					fprintf(stderr, "The requested number of variants per request is not valid (%s)\n", optarg);
				} else
				{
					options_data->variants_per_request = tmp_int_field;
				}
				dprintf("variants-per-request = %ld\n", options_data->variants_per_request);
				break;
            case 'r':
                tmp_string_field = (char*) malloc((strlen(optarg)+1) * sizeof(char));
                strcpy(tmp_string_field, optarg);
                filter = create_region_filter(tmp_string_field, 0);
                options_data->chain = add_to_filter_chain(filter, options_data->chain);
                dprintf("regions = %s\n", optarg);
                break;
			case 's':
				// options_data->species is const char*, so it must be freed and reassigned
				tmp_string_field = (char*) malloc((strlen(optarg)+1) * sizeof(char));
				strcpy(tmp_string_field, optarg);
				free((void*) options_data->species);
				options_data->species = tmp_string_field;
				dprintf("species = %s\n", options_data->species);
				break;
			case 'u':
				// options_data->url is const char*, so it must be freed and reassigned
				tmp_string_field = (char*) malloc((strlen(optarg)+1) * sizeof(char));
				strcpy(tmp_string_field, optarg);
				free((void*) options_data->host_url);
				options_data->host_url = tmp_string_field;
				dprintf("host url = %s\n", options_data->host_url);
				break;
			case 'v':
				// options_data->version is const char*, so it must be freed and reassigned
				tmp_string_field = (char*) malloc((strlen(optarg)+1) * sizeof(char));
				strcpy(tmp_string_field, optarg);
				free((void*) options_data->version);
				options_data->version = tmp_string_field;
				dprintf("version = %s\n", options_data->version);
				break;
			case '?':
			default:
				printf("Option unknown\n");
				break;
		}
		
		previous_opt_index = optind;
	}
	
	free(options);
}


int verify_effect_options(global_options_data_t *global_options_data, effect_options_data_t *options_data)
{
    // Check whether the input VCF file is defined
    if (global_options_data->vcf_filename == NULL || strlen(global_options_data->vcf_filename) == 0) {
        fprintf(stderr, "Please specify the input VCF file.\n");
        return VCF_FILE_NOT_SPECIFIED;
    }
    
	// Check whether the host URL is defined
	if (options_data->host_url == NULL || strlen(options_data->host_url) == 0)
	{
		fprintf(stderr, "Please specify the host URL to the web service.\n");
		return EFFECT_HOST_URL_NOT_SPECIFIED;
	}
	
	// Check whether the version is defined
	if (options_data->version == NULL || strlen(options_data->version) == 0)
	{
		fprintf(stderr, "Please specify the version.\n");
		return EFFECT_VERSION_NOT_SPECIFIED;
	}
	
	// Check whether the species is defined
	if (options_data->species == NULL || strlen(options_data->species) == 0)
	{
		fprintf(stderr, "Please specify the species to take as reference.\n");
		return EFFECT_SPECIES_NOT_SPECIFIED;
	}
	
	return 0;
}
