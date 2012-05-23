#ifndef SPLIT_RUNNER_H
#define SPLIT_RUNNER_H

#include "split.h"

int run_split(global_options_data_t *global_options_data, split_options_data_t *options_data);

static int initialize_output(cp_hashtable **output_files);

static void free_output(cp_hashtable *output_files);

static void free_file_key(char *key);

static void free_file_descriptor(FILE *fd);


#endif
