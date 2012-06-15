#ifndef _CP_RRTREE_H
#define _CP_RRTREE_H

/**
 * @addtogroup cp_rtrie
 */
/** @{ */

#include "common.h"

__BEGIN_DECLS

#include "config.h"

#include <string.h>
#include <errno.h>
#include "vector.h"
#include "bstr.h"
#include "mempool.h"


typedef int (*cp_rtrie_match_fn)(void *leaf); 

CPROPS_DLL struct _cp_rtrie;

/* 
 * cp_rtrie nodes can have just about 2 subnodes, which are a sequence of 1 
 * or more bits. 
 */
typedef CPROPS_DLL struct _cp_rtrie_node 
{ 
	cp_bstr *zero;
	cp_bstr *one;
	struct _cp_rtrie_node *node_zero;
	struct _cp_rtrie_node *node_one;
	void *leaf; 
} cp_rtrie_node; 

CPROPS_DLL
cp_rtrie_node *cp_rtrie_node_new(void *leaf, cp_mempool *pool); 
CPROPS_DLL
void *cp_rtrie_node_delete(struct _cp_rtrie *grp, cp_rtrie_node *node);

CPROPS_DLL
void cp_rtrie_node_unmap(struct _cp_rtrie *grp, cp_rtrie_node **node); 

/**
 * @file rtrie.h
 * cp_rtrie is a character rtrie implementation. Tries allow for prefix matching
 * with O(m) = O(1) time (m being the length of the key). Used to store key - 
 * value mappings, rtries have certain advantages over hashtables in that worse 
 * case behavior is still O(1) and no hash function is needed. cp_rtrie is 
 * technically a compact rtrie in that collapses unused character paths. 
 */
typedef CPROPS_DLL struct _cp_rtrie
{ 
	cp_rtrie_node *root;                /**< root node           */
	int path_count;                    /**< number of enrtries   */

	int mode;                          /**< collection mode     */

	cp_copy_fn copy_leaf;              /**< leaf copy function  */
	cp_destructor_fn delete_leaf;      /**< leaf destructor     */

	cp_lock *lock;                     /**< collection lock     */
	cp_thread txowner;                 /**< transaction owner   */
	int txtype;                        /**< lock type           */

	cp_mempool *mempool; 			   /**< memory pool         */
} cp_rtrie; 

/** 
 * create a new cp_rtrie object with the specified collection mode and
 * leaf management functions
 */
CPROPS_DLL
cp_rtrie *cp_rtrie_create_rtrie(int mode, 
		                     cp_copy_fn copy_leaf, 
							 cp_destructor_fn delete_leaf);

/** create a new cp_rtrie object with the specified collection mode */
CPROPS_DLL
cp_rtrie *cp_rtrie_create(int mode);
/** delete a cp_rtrie object */
CPROPS_DLL
int cp_rtrie_destroy(cp_rtrie *grp); 
/** add a mapping to a rtrie */
CPROPS_DLL
int cp_rtrie_add(cp_rtrie *grp, cp_bstr *key, void *leaf);
/** remove a mapping from a rtrie */
CPROPS_DLL
int cp_rtrie_remove(cp_rtrie *grp, cp_bstr *key, void **leaf);
/** return the mapping for the longest prefix of the given key */
CPROPS_DLL
int cp_rtrie_prefix_match(cp_rtrie *grp, cp_bstr *key, void **leaf);
/** return the mapping for the given key if any */
CPROPS_DLL
void *cp_rtrie_exact_match(cp_rtrie *grp, cp_bstr *key);
/** return a vector containing exact match and any prefix matches */
CPROPS_DLL
cp_vector *cp_rtrie_fetch_matches(cp_rtrie *grp, cp_bstr *key);
/** return a vector containing all enrtries in subtree under path given by key */
CPROPS_DLL
cp_vector *cp_rtrie_submatch(cp_rtrie *grp, cp_bstr *key);

/** return the number of stored items */
CPROPS_DLL
int cp_rtrie_count(cp_rtrie *grp);

CPROPS_DLL
void cp_rtrie_set_root(cp_rtrie *grp, void *leaf); 

CPROPS_DLL
int cp_rtrie_lock(cp_rtrie *grp, int type);
#define cp_rtrie_rdlock(grp) (cp_rtrie_lock(grp, COLLECTION_LOCK_READ))
#define cp_rtrie_wrlock(grp) (cp_rtrie_lock(grp, COLLECTION_LOCK_WRITE))
CPROPS_DLL
int cp_rtrie_unlock(cp_rtrie *grp);

/* get the current collection mode */
CPROPS_DLL
int cp_rtrie_get_mode(cp_rtrie *grp);
/* sets the bits defined by mode on the rtrie mode */
CPROPS_DLL
int cp_rtrie_set_mode(cp_rtrie *grp, int mode);
/* clears the bits defined by mode on the rtrie mode */
CPROPS_DLL
int cp_rtrie_unset_mode(cp_rtrie *grp, int mode);

CPROPS_DLL
void cp_rtrie_dump(cp_rtrie *grp);

/* set rtrie to use given mempool or allocate a new one if pool is NULL */
CPROPS_DLL
int cp_rtrie_use_mempool(cp_rtrie *tree, cp_mempool *pool);

/* set rtrie to use a shared memory pool */
CPROPS_DLL
int cp_rtrie_share_mempool(cp_rtrie *tree, cp_shared_mempool *pool);

__END_DECLS

/** @} */

#endif

