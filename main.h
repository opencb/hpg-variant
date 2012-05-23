#ifndef MAIN_H
#define MAIN_H

/** 
 * @file main.h
 * @brief Entry point of the application
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <log.h>

#include "error.h"
#include "filter/filter.h"
#include "stats/stats.h"
#include "split/split.h"


int vcf_tool_filter(int argc, char *argv[]);

int vcf_tool_merge(int argc, char *argv[]);

int vcf_tool_sort(int argc, char *argv[]);

int vcf_tool_split(int argc, char *argv[]);

int vcf_tool_stats(int argc, char *argv[]);


#endif
