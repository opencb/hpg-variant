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


int vcf_tool_filter(int argc, char *argv[]);

int vcf_tool_merge(void);

int vcf_tool_sort(void);

int vcf_tool_split(void);

int vcf_tool_stats(int argc, char *argv[]);


#endif
