#include "effect_runner.h"

void write_summary_file(cp_hashtable *summary_count, FILE *summary_file) {
     char *consequence_type;
     int *count;
     
     char **keys = (char**) cp_hashtable_get_keys(summary_count);
     int num_keys = cp_hashtable_count(summary_count);
     for (int i = 0; i < num_keys; i++) {
         consequence_type = keys[i];
         count = (int*) cp_hashtable_get(summary_count, consequence_type);
         fprintf(summary_file, "%s\t%d\n", consequence_type, *count);
     }
     free(keys);
}

void write_genes_with_variants_file(cp_hashtable *gene_list, char *output_directory) {
    char filename[] = "genes_with_variants.txt";
    FILE *file = NULL;
    char *aux_buffer;
    
    // Create genes file
    aux_buffer = (char*) calloc (strlen(output_directory) + strlen(filename) + 2, sizeof(char));
    sprintf(aux_buffer, "%s/%s", output_directory, filename);
    file = fopen(aux_buffer, "w");
    free(aux_buffer);
    
    char **keys = (char**) cp_hashtable_get_keys(gene_list);
    int num_keys = cp_hashtable_count(gene_list);
    for (int i = 0; i < num_keys; i++) {
        fprintf(file, "%s\n", keys[i]);
    }
    
    fclose(file);
}


void write_result_file(global_options_data_t *global_options_data, effect_options_data_t *options_data, cp_hashtable *summary_count, char *output_directory) {
    char *aux_buffer, *input_vcf_buffer;
    result_file_t *result_file = NULL;
    
    // Create result file
    aux_buffer = (char*) calloc (strlen(output_directory) + strlen("result.xml") + 2, sizeof(char));
    sprintf(aux_buffer, "%s/result.xml", output_directory);
    result_file = result_file_new("v0.7", aux_buffer);
    
    // Add meta data
    time_t now;
    char day_buffer[10];
    time(&now);
    strftime(day_buffer, 15, "%Y-%m-%d", localtime(&now));
    
    result_item_t *meta_item_version = result_item_new("version", result_file->version, "SDK version", "MESSAGE", "", "", "");
    result_item_t *meta_item_date = result_item_new("date", day_buffer, "Job date", "MESSAGE", "", "", "");
    result_item_t *meta_item_tool = result_item_new("tool", "consequence-type", "Tool executed", "MESSAGE", "", "", "");
    
    result_add_meta_item(meta_item_version, result_file);
    result_add_meta_item(meta_item_date, result_file);
    result_add_meta_item(meta_item_tool, result_file);
    
    // TODO Add input data
    /* 
     * <input>
<item name="log-file" title="name of the log file, default: result.log" type="MESSAGE" tags="" style="" group="" context="">
/httpd/bioinfo/wum_sessions_v0.7/1104990/jobs/111421/job.log
</item>
<item name="no-disease" title="Excludes: Mutations, miRNA diseases" type="MESSAGE" tags="" style="" group="" context=""/>
<item name="home" title="Variant home path" type="MESSAGE" tags="" style="" group="" context="">/httpd/bioinfo/variant1.0</item>
    */
    result_item_t *input_item_tool = result_item_new("tool", "consequence-type", "tool name", "MESSAGE", "", "", "");
    result_item_t *input_item_outdir = result_item_new("outdir", global_options_data->output_directory, "outdir to save the results", "MESSAGE", "", "", "");
    // TODO log-file
    result_item_t *input_item_vcf_file = result_item_new("vcf-file", global_options_data->vcf_filename, "The VCF variant file", "MESSAGE", "", "", "");
    // TODO excludes
    result_item_t *input_item_species = result_item_new("species", global_options_data->species, "The species of the ids", "MESSAGE", "", "", "");
    // TODO home
    aux_buffer = (char*) calloc (8, sizeof(char));
    sprintf(aux_buffer, "%ld", options_data->num_threads);
    result_item_t *input_item_numthreads = result_item_new("number-threads", aux_buffer, "Number of connections to the web-service", "MESSAGE", "", "", "");
    input_vcf_buffer = (char*) calloc (strlen(global_options_data->vcf_filename), sizeof(char));
    get_filename_from_path(global_options_data->vcf_filename, input_vcf_buffer);
    result_item_t *input_item_vcf_input = result_item_new(input_vcf_buffer, input_vcf_buffer, "VCF input file", "DATA", "", "Input", "");
    
    result_add_input_item(input_item_tool, result_file);
    result_add_input_item(input_item_outdir, result_file);
    result_add_input_item(input_item_vcf_file, result_file);
    result_add_input_item(input_item_species, result_file);
    result_add_input_item(input_item_vcf_input, result_file);
    
    result_item_t *output_item;
    
    // Add output files for summaries
    if (global_options_data->output_filename != NULL && strlen(global_options_data->output_filename) > 0) {
        int filename_len = strlen(global_options_data->output_filename);
        char *passed_filename = (char*) calloc (filename_len + 10, sizeof(char));
        sprintf(passed_filename, "%s.filtered", global_options_data->output_filename);
        output_item = result_item_new(passed_filename, passed_filename, "Filtered Variants", "FILE", "", "Summary", "");
        result_add_output_item(output_item, result_file);
    } else if (options_data->chain != NULL) {
        char *passed_filename = (char*) calloc (strlen(input_vcf_buffer) + 10, sizeof(char));
        sprintf(passed_filename, "%s.filtered", input_vcf_buffer);
        output_item = result_item_new(passed_filename, passed_filename, "Filtered Variants", "FILE", "", "Summary", "");
        result_add_output_item(output_item, result_file);
    } else {
        char *passed_filename = strdup(input_vcf_buffer);
        output_item = result_item_new(passed_filename, passed_filename, "Filtered Variants", "FILE", "", "Summary", "");
        result_add_output_item(output_item, result_file);
    }
    
    output_item = result_item_new("genes_with_variants.txt", "genes_with_variants.txt", "Genes with Variants", 
                                        "FILE", "", "Summary", "");
    result_add_output_item(output_item, result_file);
    output_item = result_item_new("summary", "summary.txt", "Consequence types histogram:", 
                                        "FILE", "HISTOGRAM", "Summary", "");
    result_add_output_item(output_item, result_file);
    
    // Add output files retrieved from the WS
    char *consequence_type, *consequence_type_filename;
    int *count, total_count = 0;
    
    char **keys = (char**) cp_hashtable_get_keys(summary_count);
    int num_keys = cp_hashtable_count(summary_count);
    for (int i = 0; i < num_keys; i++) {
        consequence_type = keys[i];
        count = (int*) cp_hashtable_get(summary_count, consequence_type);
        total_count += *count;
        
        if (strcmp(consequence_type, "all_variants") && strcmp(consequence_type, "summary")) {
            consequence_type_filename = (char*) calloc (strlen(consequence_type) + 5, sizeof(char));
            strncat(consequence_type_filename, consequence_type, strlen(consequence_type));
            strncat(consequence_type_filename, ".txt", 4);
            
            aux_buffer = (char*) calloc (strlen(consequence_type) + 32, sizeof(char));
            sprintf(aux_buffer, "%s (%d)", consequence_type, *count);
            
            output_item = result_item_new(consequence_type_filename, consequence_type_filename, aux_buffer, 
                                          "FILE", "CONSEQUENCE_TYPE_VARIANTS", "Variants by Consequence Type", "");
            result_add_output_item(output_item, result_file);
        }
        
    }
    free(keys);
    
    // Add output data for all_variants
    aux_buffer = (char*) calloc (32, sizeof(char));
    sprintf(aux_buffer, "All (%d)", total_count);
    output_item = result_item_new("all_variants.txt", "all_variants.txt", aux_buffer, 
                                        "FILE", "CONSEQUENCE_TYPE_VARIANTS", "Variants by Consequence Type", "");
    result_add_output_item(output_item, result_file);
    
    // Write and free file
    result_file_write(result_file->filename, result_file);
    result_file_free(result_file);
}
