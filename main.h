#ifndef MAIN_H
#define MAIN_H

/** 
 * @file main.h
 * @brief Entry point of the application
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "global_options.h"


int call(void);

int effect(int argc, char *argv[]);

int functional_analysis(void);

int genomic_analysis(void);

int pathway_analysis(void);


#endif
