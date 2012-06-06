#ifndef MAIN_H
#define MAIN_H

/** 
 * @file main.h
 * @brief Entry point of the application
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bioformats/vcf/vcf_util.h>

#include "error.h"
#include "global_options.h"
#include "hpg_variant_utils.h"



/* ***********************
 *     Available tools   *
 * ***********************/

int call(void);

int effect(int argc, char *argv[], const char *configuration_file);

int functional_analysis(void);

int genomic_analysis(int argc, char *argv[], const char *configuration_file);

int pathway_analysis(void);


#endif
