#include "effect_runner.h"

// Output file descriptors (must be kept between function calls)
static cp_hashtable *output_files = NULL;
static FILE *all_variants_file = NULL;
static FILE *summary_file = NULL;
static FILE *snp_phenotype_file = NULL;
static FILE *mutation_phenotype_file = NULL;

// Lines of the output data in the main .txt files
static list_t *output_list;
// Consequence type counters (for summary, must be kept between web service calls)
static cp_hashtable *summary_count = NULL;
// Gene list (for genes-with-variants, must be kept between web service calls)
static cp_hashtable *gene_list = NULL;

// Line buffers and their maximum size (one per thread)
static char **line;
static char **output_line;
static int *max_line_size;

static char **snp_line;
static char **snp_output_line;
static int *snp_max_line_size;

static char **mutation_line;
static char **mutation_output_line;
static int *mutation_max_line_size;

// Output directory (non-accessible directly from CURL callback function)
static char *output_directory;
static size_t output_directory_len;

static int batch_num;


int run_effect(char **urls, shared_options_data_t *shared_options, effect_options_data_t *options_data) {
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("batches", 1, shared_options->max_batches, read_list);

    int ret_code = 0;
    double start, stop, total;
    vcf_file_t *file = vcf_open(shared_options->vcf_filename);
    
    if (!file) {
        LOG_FATAL("VCF file does not exist!\n");
    }
    
    output_directory = shared_options->output_directory;
    output_directory_len = strlen(output_directory);
    
    ret_code = create_directory(output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", output_directory);
    }
    
    // Remove all .txt files in folder
    ret_code = delete_files_by_extension(output_directory, "txt");
    if (ret_code != 0) {
        return ret_code;
    }
    
    // Initialize environment for connecting to the web service
    ret_code = init_http_environment(0);
    if (ret_code != 0) {
        return ret_code;
    }
    
    // Initialize collections of file descriptors and summary counters
    ret_code = initialize_ws_output(shared_options, options_data);
    if (ret_code != 0) {
        return ret_code;
    }
    
    // Create job.status file
    char job_status_filename[output_directory_len + 10];
    sprintf(job_status_filename, "%s/job.status", output_directory);
    FILE *job_status = new_job_status_file(job_status_filename);
    if (!job_status) {
        LOG_FATAL("Can't create job status file");
    } else {
        update_job_status_file(0, job_status);
    }
    
 
#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            ret_code = vcf_parse_batches(read_list, shared_options->batch_size, file, 0);

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) {
                LOG_ERROR_F("Error %d while reading the file %s\n", ret_code, file->filename);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            list_decr_writers(read_list);
        }
        
#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            omp_set_num_threads(shared_options->num_threads);
            
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            filter_t **filters = NULL;
            int num_filters = 0;
            if (shared_options->chain != NULL) {
                filters = sort_filter_chain(shared_options->chain, &num_filters);
            }
            FILE *passed_file = NULL, *failed_file = NULL;
            get_output_files(shared_options, &passed_file, &failed_file);
    
            start = omp_get_wtime();

            int i = 0;
            list_item_t* item = NULL;
            
            int ret_ws_0, ret_ws_1, ret_ws_2;
            while ((item = list_remove_item(read_list)) != NULL) {
                if (i == 0) {
                    // Write file format, header entries and delimiter
                    if (passed_file != NULL) { vcf_write_to_file(file, passed_file); }
                    if (failed_file != NULL) { vcf_write_to_file(file, failed_file); }
                    
                    LOG_DEBUG("VCF header written\n");
                }
                vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
                array_list_t *input_records = batch;
                array_list_t *passed_records = NULL, *failed_records = NULL;

                if (i % 20 == 0) {
                    LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->size, batch->capacity);
                }

                if (filters == NULL) {
                    passed_records = input_records;
                } else {
                    failed_records = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
                    passed_records = run_filter_chain(input_records, failed_records, filters, num_filters);
                }

                // Write records that passed to a separate file, and query the WS with them as args
                if (passed_records->size > 0) {
                    // Divide the list of passed records in ranges of size defined in config file
                    int num_chunks;
                    int *chunk_sizes;
                    int *chunk_starts = create_chunks(passed_records->size, shared_options->entries_per_thread, &num_chunks, &chunk_sizes);
                    
                    // OpenMP: Launch a thread for each range
                    #pragma omp parallel for
                    for (int j = 0; j < num_chunks; j++) {
                        LOG_DEBUG_F("[%d] WS invocation\n", omp_get_thread_num());
                        LOG_DEBUG_F("[%d] -- effect WS\n", omp_get_thread_num());
                        ret_ws_0 = invoke_effect_ws(urls[0], (vcf_record_t**) (passed_records->items + chunk_starts[j]), chunk_sizes[j], options_data->excludes);
                        if (!options_data->no_phenotypes) {
                            LOG_DEBUG_F("[%d] -- snp WS\n", omp_get_thread_num());
                            ret_ws_1 = invoke_snp_phenotype_ws(urls[1], (vcf_record_t**) (passed_records->items + chunk_starts[j]), chunk_sizes[j]);
                            LOG_DEBUG_F("[%d] -- mutation WS\n", omp_get_thread_num());
                            ret_ws_2 = invoke_mutation_phenotype_ws(urls[2], (vcf_record_t**) (passed_records->items + chunk_starts[j]), chunk_sizes[j]);
                        }
                    }
                    
                    free(chunk_starts);
                    free(chunk_sizes);
                    
                    LOG_INFO_F("*** %dth web services invocation finished\n", i);
                    
                    if (ret_ws_0 || ret_ws_1 || ret_ws_2) {
                        if (ret_ws_0) {
                            LOG_ERROR_F("Effect web service error: %s\n", get_last_http_error(ret_ws_0));
                        }
                        if (ret_ws_1) {
                            LOG_ERROR_F("SNP phenotype web service error: %s\n", get_last_http_error(ret_ws_1));
                        }
                        if (ret_ws_2) {
                            LOG_ERROR_F("Mutations phenotype web service error: %s\n", get_last_http_error(ret_ws_2));
                        }
                        
                        LOG_FATAL("Can not continue execution after a web service error occurred");
                        break;
                    }
                }
                
                // Write records that passed and failed to separate files
                if (passed_file != NULL && failed_file != NULL) {
                    if (passed_records != NULL && passed_records->size > 0) {
                        write_batch(passed_records, passed_file);
                    }
                    if (failed_records != NULL && failed_records->size > 0) {
                        write_batch(failed_records, failed_file);
                    }
                }
                
                // Free items in both lists (not their internal data)
                if (passed_records != input_records) {
                    LOG_DEBUG_F("[Batch %d] %zu passed records\n", i, passed_records->size);
                    array_list_free(passed_records, NULL);
                }
                if (failed_records) {
                    LOG_DEBUG_F("[Batch %d] %zu failed records\n", i, failed_records->size);
                    array_list_free(failed_records, NULL);
                }
                // Free batch and its contents
                vcf_batch_free(item->data_p);
                list_item_free(item);
                
                i++;
            }

            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Free resources
            if (passed_file) { fclose(passed_file); }
            if (failed_file) { fclose(failed_file); }
            
            // Free filters
            for (i = 0; i < num_filters; i++) {
                filter_t *filter = filters[i];
                filter->free_func(filter);
            }
            free(filters);
            
            // Decrease list writers count
            for (i = 0; i < shared_options->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }
        
#pragma omp section
        {
            // Thread which writes the results to all_variants, summary and one file per consequence type
            int ret = 0;
            char *line;
            list_item_t* item = NULL;
            FILE *fd = NULL;
            while ((item = list_remove_item(output_list)) != NULL) {
                line = item->data_p;
                
                // Type greater than 0: consequence type identified by its SO code
                // Type equals to -1: SNP phenotype
                // Type equals to -2: mutation phenotype
                if (item->type > 0) {
                    // Write entry in the consequence type file
                    fd = cp_hashtable_get(output_files, &(item->type));
                    int ret = fprintf(fd, "%s\n", line);
                    if (ret < 0) {
                        LOG_ERROR_F("Error writing to file: '%s'\n", line);
                    }
                    
                    // Write in all_variants
                    ret = fprintf(all_variants_file, "%s\n", line);
                    if (ret < 0) {
                        LOG_ERROR_F("Error writing to all_variants: '%s'\n", line);
                    }
                    
                } else if (item->type == SNP_PHENOTYPE) {
                    ret = fprintf(snp_phenotype_file, "%s\n", line);
                    if (ret < 0) {
                        LOG_ERROR_F("Error writing to snp_phenotypes: '%s'\n", line);
                    }
                    
                } else if (item->type == MUTATION_PHENOTYPE) {
                    ret = fprintf(mutation_phenotype_file, "%s\n", line);
                    if (ret < 0) {
                        LOG_ERROR_F("Error writing to mutation_phenotypes: '%s'\n", line);
                    }
                }
                
                free(line);
                list_item_free(item);
            }
            
        }
    }

    write_summary_file(summary_count, summary_file);
    write_genes_with_variants_file(gene_list, output_directory);
    write_result_file(shared_options, options_data, summary_count, output_directory);

    ret_code = free_ws_output(shared_options->num_threads);
    free(read_list);
    free(output_list);
    vcf_close(file);
    
    update_job_status_file(100, job_status);
    close_job_status_file(job_status);
    
    return ret_code;
}


char *compose_effect_ws_request(const char *category, const char *method, shared_options_data_t *options_data) {
    if (options_data->host_url == NULL || options_data->version == NULL || options_data->species == NULL) {
        return NULL;
    }
    
    // URL Constants
    const char *ws_root_url = "cellbase/rest/";
//     const char *ws_name_url = "genomic/variant/";//consequence_type";
    const char *ws_extra_params = "?header=false";
    
    // Length of URL parts
    const int host_url_len = strlen(options_data->host_url);
    const int ws_root_len = strlen(ws_root_url);
    const int version_len = strlen(options_data->version);
    const int species_len = strlen(options_data->species);
    const int ws_name_len = strlen(category);
    const int method_len = strlen(method);
    const int ws_extra_len = strlen(ws_extra_params);
    const int result_len = host_url_len + ws_root_len + version_len + species_len + ws_name_len + method_len + ws_extra_len + 5;
    
    char *result_url = (char*) calloc (result_len, sizeof(char));
    
    // Host URL
    strncat(result_url, options_data->host_url, host_url_len);
    if (result_url[host_url_len - 1] != '/') {
        strncat(result_url, "/", 1);
    }
    
    // Root of the web service
    strncat(result_url, ws_root_url, ws_root_len);
    
    // Version
    strncat(result_url, options_data->version, version_len);
    if (result_url[strlen(result_url) - 1] != '/') {
        strncat(result_url, "/", 1);
    }
    
    // Species
    strncat(result_url, options_data->species, species_len);
    if (result_url[strlen(result_url) - 1] != '/') {
        strncat(result_url, "/", 1);
    }
    
    // Name of the web service
    strncat(result_url, category, ws_name_len);
    strncat(result_url, "/", 1);
    strncat(result_url, method, method_len);
    
    // Extra arguments of the web service
    strncat(result_url, ws_extra_params, ws_extra_len);
    
    return result_url;
}

int invoke_effect_ws(const char *url, vcf_record_t **records, int num_records, char *excludes) {
    CURL *curl;
    CURLcode ret_code = CURLE_OK;

    struct curl_httppost *formpost = NULL;
    struct curl_httppost *lastptr = NULL;
    
    const char *output_format = "txt";
    
    int variants_len = 512, current_index = 0;
    char *variants = (char*) calloc (variants_len, sizeof(char));
    
    int chr_len, position_len, reference_len, alternate_len;
    int new_len_range;

    LOG_DEBUG_F("[%d] WS for batch #%d\n", omp_get_thread_num(), batch_num);
    batch_num++;
    
    for (int i = 0; i < num_records; i++) {
        vcf_record_t *record = records[i];
        chr_len = strlen(record->chromosome);
        reference_len = strlen(record->reference);
        alternate_len = strlen(record->alternate);
        new_len_range = current_index + chr_len + reference_len + alternate_len + 32;
        
        LOG_DEBUG_F("%s:%lu:%s:%s\n", record->chromosome, record->position, record->reference, record->alternate);
        
        // Reallocate memory if next record won't fit
        if (variants_len < (current_index + new_len_range + 1)) {
            char *aux = (char*) realloc(variants, (variants_len + new_len_range + 1) * sizeof(char));
            if (aux) { 
                variants = aux; 
                variants_len += new_len_range;
            }
        }
        
        // Append region info to buffer
        strncat(variants, record->chromosome, chr_len);
        strncat(variants, ":", 1);
        current_index += chr_len + 1;
        sprintf(variants + current_index, "%lu:", record->position);
        strncat(variants, record->reference, reference_len);
        strncat(variants, ":", 1);
        strncat(variants, record->alternate, alternate_len);
        strncat(variants, ",", 1);
        current_index = strlen(variants);
    }
    
    LOG_DEBUG_F("variants = %s\n", variants);
    LOG_DEBUG_F("excludes = %s\n", excludes);
    
    char *params[CONSEQUENCE_TYPE_WS_NUM_PARAMS] = { "of", "variants", "exclude" };
    char *params_values[CONSEQUENCE_TYPE_WS_NUM_PARAMS] = { output_format, variants, excludes };
    
    ret_code = http_post(url, params, params_values, CONSEQUENCE_TYPE_WS_NUM_PARAMS, write_effect_ws_results);
    
    free(variants);
    
    return ret_code;
}

static size_t write_effect_ws_results(char *contents, size_t size, size_t nmemb, void *userdata) {
    int tid = omp_get_thread_num();
    
    int i = 0;
    int data_read_len = 0, next_line_len = 0;
    // Whether the SO code field (previous to the consequence type name) has been found
    int *SO_found = (int*) malloc (sizeof(int));
    // Whether the buffer was consumed with a line read just partially
    int premature_end = 0;
    
    size_t realsize = size * nmemb;
    
    int *count;
    
    char *data = contents;
    char tmp_consequence_type[128];
    char *aux_buffer;
    char *output_text;
    
    
    LOG_DEBUG_F("Effect WS invoked, response size = %zu bytes\n", realsize);
    
    while (data_read_len < realsize) {
        assert((line + tid) != NULL);
        assert((max_line_size + tid) != NULL);
        
        LOG_DEBUG_F("[%d] loop iteration #%d\n", tid, i);
        // Get length of data to copy
        next_line_len = strcspn(data, "\n");
        
        // If the line[tid] is too long for the current buffers, reallocate a little more than the needed memory
        if (strlen(line[tid]) + next_line_len + 1 > max_line_size[tid]) {
            LOG_DEBUG_F("Line too long (%d elements, but %zu needed) in batch #%d\n", 
                        max_line_size[tid], strlen(line[tid]) + next_line_len, batch_num);
//             char *out_buf = (char*) calloc (next_line_len+1, sizeof(char));
//             snprintf(out_buf, next_line_len, "%s", data);
//             LOG_INFO_F("[%d] too big data is: '%s'\n", tid, out_buf);
            char *aux_1 = (char*) realloc (line[tid], (max_line_size[tid] + next_line_len + 1) * sizeof(char));
            char *aux_2 = (char*) realloc (output_line[tid], (max_line_size[tid] + next_line_len + 1) * sizeof(char));
            
            if (!aux_1 || !aux_2) {
                LOG_ERROR("Can't resize buffers\n");
                // Can't resize buffers -> can't keep reading the file
                if (!aux_1) { free(line[tid]); }
                if (!aux_2) { free(output_line[tid]); }
                return data_read_len;
            }
            
            line[tid] = aux_1;
            output_line[tid] = aux_2;
            max_line_size[tid] += next_line_len + 1;
            LOG_DEBUG_F("[%d] buffers realloc'd (%d)\n", tid, max_line_size[tid]);
        }
        
        LOG_DEBUG_F("[%d] position = %d, read = %d, max_size = %zu\n", 
                    i, next_line_len, data_read_len, realsize);
        
        if (data_read_len + next_line_len >= realsize) {
            // Save current state (line[tid] partially read)
            strncat(line[tid], data, next_line_len);
            chomp(line[tid]);
            line[tid][strlen(line[tid])] = '\0';
            premature_end = 1;
            LOG_DEBUG_F("widow line[tid] = '%s'\n", line[tid]);
            data_read_len = realsize;
            break;
        }
        
        strncat(line[tid], data, next_line_len);
        strncat(output_line[tid], line[tid], strlen(line[tid]));
     
        LOG_DEBUG_F("[%d] copy to buffer (%zu)\n", tid, strlen(line[tid]));
    
        int num_substrings;
        char *copy_buf = strdup(line[tid]);
        char **split_result = split(copy_buf, "\t", &num_substrings);
        free(copy_buf);
        
        // Find consequence type name (always after SO field)
        *SO_found = 0;
        if (num_substrings == 25) {
            LOG_DEBUG_F("gene = %s\tSO = %d\tCT = %s\n", split_result[17], atoi(split_result[18] + 3), split_result[19]);
            if (!cp_hashtable_contains(gene_list, split_result[17])) {
                cp_hashtable_put(gene_list, strdup(split_result[17]), NULL);
            }
            *SO_found = atoi(split_result[18] + 3);
           memset(tmp_consequence_type, 0, 128 * sizeof(char));
           strncat(tmp_consequence_type, split_result[19], strlen(split_result[19]));
        } else {
            LOG_INFO_F("[%d] Non-valid line found: '%s'\n", tid, line[tid]);
            memset(line[tid], 0, strlen(line[tid]));
            memset(output_line[tid], 0, strlen(output_line[tid]));
            
            for (int s = 0; s < num_substrings; s++) {
                free(split_result[s]);
            }
            free(split_result);
            continue;
        }
        
        for (int s = 0; s < num_substrings; s++) {
            free(split_result[s]);
        }
        free(split_result);
        
        if (!*SO_found) { // SO:000000 is not valid
            LOG_DEBUG_F("[%d] Non-valid SO found (0)\n", tid);
            memset(line[tid], 0, strlen(line[tid]));
            memset(output_line[tid], 0, strlen(output_line[tid]));
            continue;
        }

        LOG_DEBUG_F("[%d] SO found = %d\n", tid, *SO_found);
        size_t consequence_type_len = strlen(tmp_consequence_type);
     
        // If file does not exist, create its descriptor and summary counter
        FILE *aux_file = cp_hashtable_get(output_files, SO_found);
        if (!aux_file) {
#pragma omp critical
            {
                // This construction avoids 2 threads trying to insert the same CT
                aux_file = cp_hashtable_get(output_files, SO_found);
                if (!aux_file) {
                    char filename[output_directory_len + consequence_type_len + 6];
                    memset(filename, 0, (output_directory_len + consequence_type_len + 6) * sizeof(char));
                    strncat(filename, output_directory, output_directory_len);
                    strncat(filename, "/", 1);
                    strncat(filename, tmp_consequence_type, consequence_type_len);
                    strncat(filename, ".txt", 4);
                    aux_file = fopen(filename, "a");
                    
                    // Add to hashtables (file descriptors and summary counters)
                    int *SO_stored = (int*) malloc (sizeof(int));
                    *SO_stored = *SO_found;
                    cp_hashtable_put(output_files, SO_stored, aux_file);

                    LOG_INFO_F("[%d] new CT = %s\n", tid, tmp_consequence_type);
                }
            }
        }
        
        // Write line[tid] to file corresponding to its consequence type
        if (aux_file) { 
#pragma omp critical
            {
                // TODO move critical one level below?
                count = (int*) cp_hashtable_get(summary_count, tmp_consequence_type);
                if (count == NULL) {
                    char *consequence_type = (char*) calloc (consequence_type_len+1, sizeof(char));
                    strncat(consequence_type, tmp_consequence_type, consequence_type_len);
                    assert(!strcmp(consequence_type, tmp_consequence_type));
                    count = (int*) malloc (sizeof(int));
                    *count = 0;
                    cp_hashtable_put(summary_count, consequence_type, count);
                    LOG_DEBUG_F("[%d] Initialized summary count for %s\n", tmp_consequence_type);
                }
                // Increment counter for summary
                (*count)++;
            }
            
            LOG_DEBUG_F("[%d] before writing %s\n", tid, tmp_consequence_type);
            output_text = strdup(output_line[tid]);
            list_item_t *output_item = list_item_new(tid, *SO_found, output_text);
            list_insert_item(output_item, output_list);
            LOG_DEBUG_F("[%d] after writing %s\n", tid, tmp_consequence_type);
        }
        
        data += next_line_len+1;
        data_read_len += next_line_len+1;
        
        memset(line[tid], 0, strlen(line[tid]));
        memset(output_line[tid], 0, strlen(output_line[tid]));
        
        i++;
    }
 
    // Empty buffer for next callback invocation
    if (!premature_end) {
        memset(line[tid], 0, strlen(line[tid]));
        memset(output_line[tid], 0, strlen(line[tid]));
    }
    free(SO_found);

    return data_read_len;
}

int invoke_snp_phenotype_ws(const char *url, vcf_record_t **records, int num_records) {
    CURL *curl;
    CURLcode ret_code = CURLE_OK;

    struct curl_httppost *formpost = NULL;
    struct curl_httppost *lastptr = NULL;
    
    const char *output_format = "txt";
    
    int variants_len = 512, current_index = 0;
    char *variants = (char*) calloc (variants_len, sizeof(char));
    
    int id_len, new_len_range;

    LOG_DEBUG_F("[%d] WS for batch #%d\n", omp_get_thread_num(), batch_num);
    batch_num++;
    
    for (int i = 0; i < num_records; i++) {
        vcf_record_t *record = records[i];
        if (!strcmp(".", record->id)) {
            continue;
        }
        
        id_len = strlen(record->id);
        new_len_range = current_index + id_len + 32;
        
        LOG_DEBUG_F("%s:%lu:%s:%s\n", record->chromosome, record->position, record->reference, record->alternate);
        
        // Reallocate memory if next record won't fit
        if (variants_len < (current_index + new_len_range + 1)) {
            char *aux = (char*) realloc(variants, (variants_len + new_len_range + 1) * sizeof(char));
            if (aux) { 
                variants = aux; 
                variants_len += new_len_range;
            }
        }
        
        // Append region info to buffer
        strncat(variants, record->id, id_len);
        strncat(variants, ",", 1);
        current_index += id_len + 2;
    }
    
    LOG_DEBUG_F("snps = %s\n", variants);
    
    char *params[CONSEQUENCE_TYPE_WS_NUM_PARAMS-1] = { "of", "snps" };
    char *params_values[CONSEQUENCE_TYPE_WS_NUM_PARAMS-1] = { output_format, variants };
    
    ret_code = http_post(url, params, params_values, CONSEQUENCE_TYPE_WS_NUM_PARAMS-1, write_snp_phenotype_ws_results);
    
    free(variants);
    
    return ret_code;
}

int invoke_mutation_phenotype_ws(const char *url, vcf_record_t **records, int num_records) {
    CURL *curl;
    CURLcode ret_code = CURLE_OK;

    struct curl_httppost *formpost = NULL;
    struct curl_httppost *lastptr = NULL;
    
    const char *output_format = "txt";
    
    int variants_len = 512, current_index = 0;
    char *variants = (char*) calloc (variants_len, sizeof(char));
    
    int chr_len, position_len, reference_len, alternate_len;
    int new_len_range;

    LOG_DEBUG_F("[%d] WS for batch #%d\n", omp_get_thread_num(), batch_num);
    batch_num++;
    
    for (int i = 0; i < num_records; i++) {
        vcf_record_t *record = records[i];
        chr_len = strlen(record->chromosome);
        reference_len = strlen(record->reference);
        alternate_len = strlen(record->alternate);
        new_len_range = current_index + chr_len + reference_len + alternate_len + 32;
        
        LOG_DEBUG_F("mutation phenotype of %s:%lu:%s:%s\n", record->chromosome, record->position, record->reference, record->alternate);
        
        // Reallocate memory if next record won't fit
        if (variants_len < (current_index + new_len_range + 1)) {
            char *aux = (char*) realloc(variants, (variants_len + new_len_range + 1) * sizeof(char));
            if (aux) { 
                variants = aux; 
                variants_len += new_len_range;
            }
        }
        
        // Append region info to buffer
        strncat(variants, record->chromosome, chr_len);
        strncat(variants, ":", 1);
        current_index += chr_len + 1;
        sprintf(variants + current_index, "%lu:", record->position);
        strncat(variants, record->reference, reference_len);
        strncat(variants, ":", 1);
        strncat(variants, record->alternate, alternate_len);
        strncat(variants, ",", 1);
        current_index = strlen(variants);
    }
    
    LOG_DEBUG_F("variants = %s\n", variants);
    
    char *params[CONSEQUENCE_TYPE_WS_NUM_PARAMS-1] = { "of", "variants" };
    char *params_values[CONSEQUENCE_TYPE_WS_NUM_PARAMS-1] = { output_format, variants };
    
    ret_code = http_post(url, params, params_values, CONSEQUENCE_TYPE_WS_NUM_PARAMS-1, write_mutation_phenotype_ws_results);
    
    free(variants);
    
    return ret_code;
}

static size_t write_snp_phenotype_ws_results(char *contents, size_t size, size_t nmemb, void *userdata) {
    int tid = omp_get_thread_num();
    
    int i = 0;
    int data_read_len = 0, next_line_len = 0;
    // Whether the buffer was consumed with a line read just partially
    int premature_end = 0;
    
    size_t realsize = size * nmemb;
    
    char *data = contents;
    char *output_text;
    
    
    LOG_DEBUG_F("SNP phenotype WS invoked, response size = %zu bytes\n", realsize);
    
    while (data_read_len < realsize) {
        assert((snp_line + tid) != NULL);
        assert((snp_max_line_size + tid) != NULL);
        
        LOG_DEBUG_F("[%d] loop iteration #%d\n", tid, i);
        // Get length of data to copy
        next_line_len = strcspn(data, "\n");
        
        // If the snp_line[tid] is too long for the current buffers, reallocate a little more than the needed memory
        if (strlen(snp_line[tid]) + next_line_len + 1 > snp_max_line_size[tid]) {
            LOG_DEBUG_F("Line too long (%d elements, but %zu needed) in batch #%d\n", 
                        snp_max_line_size[tid], strlen(snp_line[tid]) + next_line_len, batch_num);
//             char *out_buf = (char*) calloc (next_line_len+1, sizeof(char));
//             snprintf(out_buf, next_line_len, "%s", data);
//             LOG_INFO_F("[%d] too big data is: '%s'\n", tid, out_buf);
            char *aux_1 = (char*) realloc (snp_line[tid], (snp_max_line_size[tid] + next_line_len + 1) * sizeof(char));
            char *aux_2 = (char*) realloc (snp_output_line[tid], (snp_max_line_size[tid] + next_line_len + 1) * sizeof(char));
            
            if (!aux_1 || !aux_2) {
                LOG_ERROR("Can't resize buffers\n");
                // Can't resize buffers -> can't keep reading the file
                if (!aux_1) { free(snp_line[tid]); }
                if (!aux_2) { free(snp_output_line[tid]); }
                return data_read_len;
            }
            
            snp_line[tid] = aux_1;
            snp_output_line[tid] = aux_2;
            snp_max_line_size[tid] += next_line_len + 1;
            LOG_DEBUG_F("[%d] buffers realloc'd (%d)\n", tid, snp_max_line_size[tid]);
        }
        
        LOG_DEBUG_F("[%d] position = %d, read = %d, max_size = %zu\n", 
                    i, next_line_len, data_read_len, realsize);
        
        if (data_read_len + next_line_len >= realsize) {
            // Save current state (snp_line[tid] partially read)
            strncat(snp_line[tid], data, next_line_len);
            chomp(snp_line[tid]);
            snp_line[tid][strlen(snp_line[tid])] = '\0';
            premature_end = 1;
            LOG_DEBUG_F("widow snp_line[tid] = '%s'\n", snp_line[tid]);
            data_read_len = realsize;
            break;
        }
        
        strncat(snp_line[tid], data, next_line_len);
        strncat(snp_output_line[tid], snp_line[tid], strlen(snp_line[tid]));
     
        LOG_DEBUG_F("[%d] copy to buffer (%zu)\n", tid, strlen(snp_line[tid]));
    
        LOG_DEBUG_F("[%d] before writing snp phenotype\n", tid);
        output_text = strdup(snp_output_line[tid]);
        list_item_t *output_item = list_item_new(tid, SNP_PHENOTYPE, output_text);
        list_insert_item(output_item, output_list);
        LOG_DEBUG_F("[%d] after writing snp phenotype\n", tid);
            
        data += next_line_len+1;
        data_read_len += next_line_len+1;
        
        memset(snp_line[tid], 0, strlen(snp_line[tid]));
        memset(snp_output_line[tid], 0, strlen(snp_output_line[tid]));
        
        i++;
    }
 
    // Empty buffer for next callback invocation
    if (!premature_end) {
        memset(snp_line[tid], 0, strlen(snp_line[tid]));
        memset(snp_output_line[tid], 0, strlen(snp_line[tid]));
    }

    return data_read_len;
}

static size_t write_mutation_phenotype_ws_results(char *contents, size_t size, size_t nmemb, void *userdata) {
    int tid = omp_get_thread_num();
    
    int i = 0;
    int data_read_len = 0, next_line_len = 0;
    // Whether the buffer was consumed with a line read just partially
    int premature_end = 0;
    
    size_t realsize = size * nmemb;
    
    char *data = contents;
    char *output_text;
    
    
    LOG_DEBUG_F("Mutation phenotype WS invoked, response size = %zu bytes -> %s\n", realsize, data);
    
    while (data_read_len < realsize) {
        assert((mutation_line + tid) != NULL);
        assert((mutation_max_line_size + tid) != NULL);
        
        LOG_DEBUG_F("[%d] loop iteration #%d\n", tid, i);
        // Get length of data to copy
        next_line_len = strcspn(data, "\n");
        
        // If the mutation_line[tid] is too long for the current buffers, reallocate a little more than the needed memory
        if (strlen(mutation_line[tid]) + next_line_len + 1 > mutation_max_line_size[tid]) {
            LOG_DEBUG_F("Line too long (%d elements, but %zu needed) in batch #%d\n", 
                        mutation_max_line_size[tid], strlen(mutation_line[tid]) + next_line_len, batch_num);
//             char *out_buf = (char*) calloc (next_line_len+1, sizeof(char));
//             snprintf(out_buf, next_line_len, "%s", data);
//             LOG_INFO_F("[%d] too big data is: '%s'\n", tid, out_buf);
            char *aux_1 = (char*) realloc (mutation_line[tid], (mutation_max_line_size[tid] + next_line_len + 1) * sizeof(char));
            char *aux_2 = (char*) realloc (mutation_output_line[tid], (mutation_max_line_size[tid] + next_line_len + 1) * sizeof(char));
            
            if (!aux_1 || !aux_2) {
                LOG_ERROR("Can't resize buffers\n");
                // Can't resize buffers -> can't keep reading the file
                if (!aux_1) { free(mutation_line[tid]); }
                if (!aux_2) { free(mutation_output_line[tid]); }
                return data_read_len;
            }
            
            mutation_line[tid] = aux_1;
            mutation_output_line[tid] = aux_2;
            mutation_max_line_size[tid] += next_line_len + 1;
            LOG_DEBUG_F("[%d] buffers realloc'd (%d)\n", tid, mutation_max_line_size[tid]);
        }
        
        LOG_DEBUG_F("[%d] position = %d, read = %d, max_size = %zu\n", 
                    i, next_line_len, data_read_len, realsize);
        
        if (data_read_len + next_line_len >= realsize) {
            // Save current state (mutation_line[tid] partially read)
            strncat(mutation_line[tid], data, next_line_len);
            chomp(mutation_line[tid]);
            mutation_line[tid][strlen(mutation_line[tid])] = '\0';
            premature_end = 1;
            LOG_DEBUG_F("widow mutation_line[tid] = '%s'\n", mutation_line[tid]);
            data_read_len = realsize;
            break;
        }
        
        strncat(mutation_line[tid], data, next_line_len);
        strncat(mutation_output_line[tid], mutation_line[tid], strlen(mutation_line[tid]));
     
        LOG_DEBUG_F("[%d] copy to buffer (%zu)\n", tid, strlen(mutation_line[tid]));
    
        LOG_DEBUG_F("[%d] before writing mutation phenotype\n", tid);
        output_text = strdup(mutation_output_line[tid]);
        list_item_t *output_item = list_item_new(tid, MUTATION_PHENOTYPE, output_text);
        list_insert_item(output_item, output_list);
        LOG_DEBUG_F("[%d] after writing mutation phenotype\n", tid);
            
        data += next_line_len+1;
        data_read_len += next_line_len+1;
        
        memset(mutation_line[tid], 0, strlen(mutation_line[tid]));
        memset(mutation_output_line[tid], 0, strlen(mutation_output_line[tid]));
        
        i++;
    }
 
    // Empty buffer for next callback invocation
    if (!premature_end) {
        memset(mutation_line[tid], 0, strlen(mutation_line[tid]));
        memset(mutation_output_line[tid], 0, strlen(mutation_line[tid]));
    }

    return data_read_len;
}

int initialize_ws_output(shared_options_data_t *shared_options, effect_options_data_t *options_data){
    int num_threads = shared_options->num_threads;
    char *outdir = shared_options->output_directory;
    
    // Initialize output text list
    output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", num_threads, shared_options->max_batches * shared_options->batch_size, output_list);
    
    // Initialize collections of file descriptors
    output_files = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                                 50,
                                                 cp_hash_int,
                                                 (cp_compare_fn) int_cmp,
                                                 NULL,
                                                 (cp_destructor_fn) free_file_key1,
                                                 NULL,
                                                 (cp_destructor_fn) free_file_descriptor
                                                );
    
    char all_variants_filename[strlen(outdir) + 18];
    sprintf(all_variants_filename, "%s/all_variants.txt", outdir);
    all_variants_file = fopen(all_variants_filename, "a");
    if (!all_variants_file) {   // Can't store results
        return 1;
    }
    char *key = (char*) calloc (13, sizeof(char));
    strncat(key, "all_variants", 12);
    cp_hashtable_put(output_files, key, all_variants_file);
    
    char summary_filename[strlen(outdir) + 13];
    sprintf(summary_filename, "%s/summary.txt", outdir);
    summary_file = fopen(summary_filename, "a");
    if (!summary_file) {   // Can't store results
        return 2;
    }
    key = (char*) calloc (8, sizeof(char));
    strncat(key, "summary", 7);
    cp_hashtable_put(output_files, key, summary_file);
    
    char snp_phenotype_filename[strlen(outdir) + 20];
    sprintf(snp_phenotype_filename, "%s/snp_phenotypes.txt", outdir);
    snp_phenotype_file = fopen(snp_phenotype_filename, "a");
    if (!snp_phenotype_filename) {
        return 3;
    }
    key = (char*) calloc (15, sizeof(char));
    strncat(key, "snp_phenotypes", 14);
    cp_hashtable_put(output_files, key, snp_phenotype_file);
    
    char mutation_phenotype_filename[strlen(outdir) + 25];
    sprintf(mutation_phenotype_filename, "%s/mutation_phenotypes.txt", outdir);
    mutation_phenotype_file = fopen(mutation_phenotype_filename, "a");
    if (!mutation_phenotype_filename) {
        return 3;
    }
    key = (char*) calloc (20, sizeof(char));
    strncat(key, "mutation_phenotypes", 19);
    cp_hashtable_put(output_files, key, mutation_phenotype_file);
    
    
    // Initialize summary counters and genes list
    summary_count = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                                 64,
                                                 cp_hash_istring,
                                                 (cp_compare_fn) strcasecmp,
                                                 NULL,
                                                 (cp_destructor_fn) free_file_key2,
                                                 NULL,
                                                 (cp_destructor_fn) free_summary_counter
                                                 );

    gene_list = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                              64,
                                              cp_hash_istring,
                                              (cp_compare_fn) strcasecmp,
                                              NULL,
                                              (cp_destructor_fn) free_file_key2,
                                              NULL,
                                              NULL
                                             );
    
    // Create a buffer for each thread
    line = (char**) calloc (num_threads, sizeof(char*));
    output_line = (char**) calloc (num_threads, sizeof(char*));
    max_line_size = (int*) calloc (num_threads, sizeof(int));
    
    snp_line = (char**) calloc (num_threads, sizeof(char*));
    snp_output_line = (char**) calloc (num_threads, sizeof(char*));
    snp_max_line_size = (int*) calloc (num_threads, sizeof(int));
    
    mutation_line = (char**) calloc (num_threads, sizeof(char*));
    mutation_output_line = (char**) calloc (num_threads, sizeof(char*));
    mutation_max_line_size = (int*) calloc (num_threads, sizeof(int));
    
    for (int i = 0; i < num_threads; i++) {
        max_line_size[i] = snp_max_line_size[i] = mutation_max_line_size[i] = 512;
        
        line[i] = (char*) calloc (max_line_size[i], sizeof(char));
        output_line[i] = (char*) calloc (max_line_size[i], sizeof(char));
        snp_line[i] = (char*) calloc (snp_max_line_size[i], sizeof(char));
        snp_output_line[i] = (char*) calloc (snp_max_line_size[i], sizeof(char));
        mutation_line[i] = (char*) calloc (mutation_max_line_size[i], sizeof(char));
        mutation_output_line[i] = (char*) calloc (mutation_max_line_size[i], sizeof(char));
    }
         
    return 0;
}

int free_ws_output(int num_threads) {
    // Free file descriptors, summary counters and gene list
    cp_hashtable_destroy(output_files);
    cp_hashtable_destroy(summary_count);
    cp_hashtable_destroy(gene_list);
    
    // Free line buffers
    for (int i = 0; i < num_threads; i++) {
        free(line[i]);
        free(output_line[i]);
        free(snp_line[i]);
        free(snp_output_line[i]);
        free(mutation_line[i]);
        free(mutation_output_line[i]);
    }
    
    free(max_line_size);
    free(line);
    free(output_line);
        
    free(snp_max_line_size);
    free(snp_line);
    free(snp_output_line);
        
    free(mutation_max_line_size);
    free(mutation_line);
    free(mutation_output_line);
    
    return 0;
}

static void free_file_key1(int *key) {
    LOG_DEBUG_F("Free file key 1: %d\n", *key);
    free(key);
}

static void free_file_descriptor(FILE *fd) {
    LOG_DEBUG("Free file descriptor\n");
    fclose(fd);
}

static void free_file_key2(char *key) {
    LOG_DEBUG_F("Free file key 2: %s\n", key);
    free(key);
}

static void free_summary_counter(int *count) {
    LOG_DEBUG_F("Free summary counter %d\n", *count);
    free(count);
}

static int int_cmp(int *a, int *b) {
    return *a - *b;
}
