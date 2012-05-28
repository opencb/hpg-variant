#ifndef HPG_VARIANT_ERROR_H
#define HPG_VARIANT_ERROR_H

// General errors (range 1-100 available)
#define NOT_IMPLEMENTED_TOOL                    1
#define CANT_READ_CONFIG_FILE                   2

#define VCF_FILE_NOT_SPECIFIED                  10
#define PED_FILE_NOT_SPECIFIED                  11


// Effect tool errors
#define EFFECT_HOST_URL_NOT_SPECIFIED           100
#define EFFECT_VERSION_NOT_SPECIFIED            101
#define EFFECT_SPECIES_NOT_SPECIFIED            102
#define EFFECT_REGIONS_NOT_SPECIFIED            103


// GWAS tool errors
#define GWAS_TASK_NOT_SPECIFIED                 200

#endif
