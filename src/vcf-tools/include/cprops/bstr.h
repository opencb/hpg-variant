#ifndef _CP_BSTR_H
#define _CP_BSTR_H

#include "collection.h"

typedef struct _cp_bstr
{
	unsigned char *bits;
	int length;
} cp_bstr;

#define BYTECOUNT(b) (((b)->length + 7) >> 3)

CPROPS_DLL
cp_bstr *cp_bstr_create(int length, unsigned char *bits);

/**
 * convenience function for debugging - initialize from a string that looks 
 * like "1010101" (actually any non-null character but '1' is interpreted as '0')
 */
cp_bstr *cstr_to_bstr(char *str);

CPROPS_DLL
void cp_bstr_destroy(cp_bstr *seq);

CPROPS_DLL
cp_bstr *cp_bstr_dup(cp_bstr *seq);

CPROPS_DLL
cp_bstr *cp_bstr_cpy(cp_bstr *dst, cp_bstr *src);

CPROPS_DLL
cp_bstr *cp_bstr_cat(cp_bstr *head, cp_bstr *tail);

/**
 * shift bits left by the specified count, resulting in a shorter bit sequence
 */
CPROPS_DLL
int cp_bstr_shift_left(cp_bstr *seq, int count);

/**
 * compares bit sequences a and b, strcmp semantics. if pos is not null it is
 * set to the index of the first differing bit. if the sequences are identical
 * to the length they are defined, and pos is not null, it is set to the length
 * of the shorter sequence - e.g. 010 and 0101 are identical to the third bit, 
 * hence *pos is set to 3.
 */
CPROPS_DLL
int cp_bstr_cmp(cp_bstr *a, cp_bstr *b, int *pos);

#define cp_bstr_length(seq) (seq)->length

CPROPS_DLL
void cp_bstr_dump(cp_bstr *seq);

CPROPS_DLL
char *cp_bstr_to_string(cp_bstr *seq);

#endif /* _CP_BSTR_H */

