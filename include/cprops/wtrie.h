#ifndef _CP_WRTREE_H
#define _CP_WRTREE_H

/**
 * @addtogroup cp_wtrie
 */
/** @{ */

#include "common.h"

__BEGIN_DECLS

#include "config.h"

#include <string.h>
#include <errno.h>
#include "hashtable.h"
#include "vector.h"
#include "wtab.h"
#include "mempool.h"

#define WNODE_MATCH(n, i) ((n)->others ? wtab_get((n)->others, *i) : NULL)
#define BRANCH_COUNT(node) wtab_count((node)->others)

typedef int (*cp_wtrie_match_fn)(void *leaf); 

CPROPS_DLL struct _cp_wtrie;

/* 
 * cp_wtrie nodes can have any number of subnodes mapped by an wtab - a hash
 * table designed for wide character keys
 */
typedef CPROPS_DLL struct _cp_wtrie_node 
{ 
	wtab *others; 
	void *leaf; 
} cp_wtrie_node; 

CPROPS_DLL
cp_wtrie_node *cp_wtrie_node_new(void *leaf, cp_mempool *pool); 
CPROPS_DLL
void *cp_wtrie_node_delete(struct _cp_wtrie *grp, cp_wtrie_node *node);
CPROPS_DLL
void cp_wtrie_delete_mapping(struct _cp_wtrie *grp, wtab_node *map_node);

CPROPS_DLL
void cp_wtrie_node_unmap(struct _cp_wtrie *grp, cp_wtrie_node **node); 

/**
 * @file wtrie.h
 * cp_wtrie is a wide character trie implementation. Tries allow for prefix matching
 * with O(m) = O(1) time (m being the length of the key). Used to store key - 
 * value mappings, tries have certain advantages over hashtables in that worse 
 * case behavior is still O(1) and no hash function is needed. cp_wtrie is 
 * technically a compact trie in that collapses unused character paths. 
 */
typedef CPROPS_DLL struct _cp_wtrie
{ 
	cp_wtrie_node *root;                /**< root node           */
	int path_count;                    /**< number of enwtries   */

	int mode;                          /**< collection mode     */

	cp_copy_fn copy_leaf;              /**< leaf copy function  */
	cp_destructor_fn delete_leaf;      /**< leaf destructor     */

	cp_lock *lock;                     /**< collection lock     */
	cp_thread txowner;                 /**< transaction owner   */
	int txtype;                        /**< lock type           */

	cp_mempool *mempool; 			   /**< memory pool         */
} cp_wtrie; 

/** 
 * create a new cp_wtrie object with the specified collection mode and
 * leaf management functions
 */
CPROPS_DLL
cp_wtrie *cp_wtrie_create_wtrie(int mode, 
		                     cp_copy_fn copy_leaf, 
							 cp_destructor_fn delete_leaf);

/** create a new cp_wtrie object with the specified collection mode */
CPROPS_DLL
cp_wtrie *cp_wtrie_create(int mode);
/** delete a cp_wtrie object */
CPROPS_DLL
int cp_wtrie_destroy(cp_wtrie *grp); 
/** add a mapping to a wtrie */
CPROPS_DLL
int cp_wtrie_add(cp_wtrie *grp, wchar_t *key, void *leaf); 
/** remove a mapping from a wtrie */
CPROPS_DLL
int cp_wtrie_remove(cp_wtrie *grp, wchar_t *key, void **leaf); 
/** return the mapping for the longest prefix of the given key */
CPROPS_DLL
int cp_wtrie_prefix_match(cp_wtrie *grp, wchar_t *key, void **leaf);
/** return the mapping for the given key if any */
CPROPS_DLL
void *cp_wtrie_exact_match(cp_wtrie *grp, wchar_t *key);
/** return a vector containing exact match and any prefix matches */
CPROPS_DLL
cp_vector *cp_wtrie_fetch_matches(cp_wtrie *grp, wchar_t *key);
/** return a vector containing all entries in subtree under path given by key */
CPROPS_DLL
cp_vector *cp_wtrie_submatch(cp_wtrie *grp, wchar_t *key);

/** return the number of stored items */
CPROPS_DLL
int cp_wtrie_count(cp_wtrie *grp);

CPROPS_DLL
void cp_wtrie_set_root(cp_wtrie *grp, void *leaf); 

CPROPS_DLL
int cp_wtrie_lock(cp_wtrie *grp, int type);
#define cp_wtrie_rdlock(grp) (cp_wtrie_lock(grp, COLLECTION_LOCK_READ))
#define cp_wtrie_wrlock(grp) (cp_wtrie_lock(grp, COLLECTION_LOCK_WRITE))
CPROPS_DLL
int cp_wtrie_unlock(cp_wtrie *grp);

/* get the current collection mode */
CPROPS_DLL
int cp_wtrie_get_mode(cp_wtrie *grp);
/* sets the bits defined by mode on the wtrie mode */
CPROPS_DLL
int cp_wtrie_set_mode(cp_wtrie *grp, int mode);
/* clears the bits defined by mode on the wtrie mode */
CPROPS_DLL
int cp_wtrie_unset_mode(cp_wtrie *grp, int mode);

CPROPS_DLL
void cp_wtrie_dump(cp_wtrie *grp);

/* set wtrie to use given mempool or allocate a new one if pool is NULL */
CPROPS_DLL
int cp_wtrie_use_mempool(cp_wtrie *tree, cp_mempool *pool);

/* set wtrie to use a shared memory pool */
CPROPS_DLL
int cp_wtrie_share_mempool(cp_wtrie *tree, cp_shared_mempool *pool);

__END_DECLS

/** @} */

#endif

