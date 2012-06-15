#ifndef _CP_WTAB_H
#define _CP_WTAB_H

#include "common.h"
#include "collection.h"

__BEGIN_DECLS

#include "config.h"

#include <wchar.h>

/*
 * table entry descriptor. wtab entries map a character to a value, and in 
 * addition allow specifying an attribute for the mapping. this is used by 
 * cp_trie to collapse multiple single child trie edges into a single node.
 */
typedef struct _wtab_node
{
	wchar_t key;
	void *value;
	void *attr;
	struct _wtab_node *next;
} wtab_node;

wtab_node *wtab_node_new(wchar_t key, void *value, void *attr);

/* 
 * the 'owner' parameter is for use by the enclosing data structure. cp_trie
 * uses this to recursively delete a sub-tree.
 */
typedef void *(*wtab_dtr)(void *owner, wtab_node *node);

/*
 * wtab is a hash table implementation specialized for use as a holder for 
 * cp_trie edges. 
 */
typedef struct _wtab
{
	int size;
	int items;
	wtab_node **table;
} wtab;

wtab *wtab_new(int size);

void wtab_delete(wtab *t);
void wtab_delete_custom(wtab *t, 
		                void *owner, 
						wtab_dtr dtr);
wtab_node *wtab_put(wtab *t, wchar_t key, void *value, void *attr);
wtab_node *wtab_get(wtab *t, wchar_t key);
void *wtab_remove(wtab *t, wchar_t key);

int wtab_count(wtab *t);

int wtab_callback(wtab *t, cp_callback_fn fn, void *prm);

__END_DECLS

#endif
