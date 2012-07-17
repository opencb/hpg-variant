#ifndef VCF_TOOLS_ERROR_H
#define VCF_TOOLS_ERROR_H

// General errors (range 1-100 available)
#define NOT_IMPLEMENTED_TOOL                    1
#define CANT_READ_CONFIG_FILE                   2


#define VCF_FILE_NOT_SPECIFIED                  10

#define HOST_URL_NOT_SPECIFIED                  30
#define VERSION_NOT_SPECIFIED                   31
#define SPECIES_NOT_SPECIFIED                   32

// Filter tool errors
#define EMPTY_LIST_OF_FILTERS                   100

// Merge tool errors
#define MISSING_MODE_NOT_SPECIFIED              200
#define DISCORDANT_CHROMOSOME                   210
#define DISCORDANT_POSITION                     211
#define DISCORDANT_REFERENCE                    212

// Stats tool errors

// Split tool errors
#define NONE_CRITERION_SPECIFIED                400


#endif
