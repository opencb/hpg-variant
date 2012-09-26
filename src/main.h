#ifndef MAIN_H
#define MAIN_H

/** 
 * @file main.h
 * @brief Entry point of the application
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <commons/log.h>

#include "error.h"
#include "hpg_vcf_tools_utils.h"
#include "filter/filter.h"
// #include "split/split.h"
#include "stats/stats.h"


int vcf_tool_filter(int argc, char *argv[], const char *configuration_file);

int vcf_tool_merge(int argc, char *argv[], const char *configuration_file);

int vcf_tool_sort(int argc, char *argv[], const char *configuration_file);

int vcf_tool_split(int argc, char *argv[], const char *configuration_file);

int vcf_tool_stats(int argc, char *argv[], const char *configuration_file);


#endif
