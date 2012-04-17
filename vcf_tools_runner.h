#ifndef VCF_TOOLS_RUNNER_H
#define VCF_TOOLS_RUNNER_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include <vcf_batch.h>
#include <vcf_filters.h>
#include <vcf_file.h>
#include <vcf_write.h>

#include "main.h"
#include "error.h"
#include "vcf_tools.h"

int execute_vcf_tools(global_options_data_t *global_options_data, vcf_tools_options_data_t *options_data);

#endif
