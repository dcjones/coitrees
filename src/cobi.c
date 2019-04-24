
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cobi.h"

void cobi_init(node_s* nodes, size_t numnodes);


// structure used when computing VeB ordering.
typedef struct reindex_s
{
    // TODO: maybe we don't need this and can just use dfs order.
    // Like do a traversal of the sorted intervals, and assign indexes as we go.
    uint32_t i;
    uint32_t depth;
    uint32_t dfs;
} reindex_s;


static int compare_depth(const void* i_, const void* j_)
{
    reindex_s* i = (reindex_s*) i_;
    reindex_s* j = (reindex_s*) j_;
    return i->depth < j->depth ? -1 :
           i->depth > j->depth ?  1 : 0;
}


static int compare_dfs(const void* i_, const void* j_)
{
    reindex_s* i = (reindex_s*) i_;
    reindex_s* j = (reindex_s*) j_;
    return i->dfs < j->dfs ? -1 :
           i->dfs > j->dfs ?  1 : 0;
}


size_t init_reindex_subtree(
        reindex_s* idxs, size_t from, size_t to, size_t depth, size_t dfs)
{
    if (from > to) return 0;

    /*fprintf(stderr, "(%d, %d)\n", (int) from, (int) to);*/

    size_t root_idx = from + (to - from)/2;
    idxs[root_idx].depth = depth;
    idxs[root_idx].dfs = dfs++;

    if (root_idx > from) {
        dfs = init_reindex_subtree(
            idxs, from, root_idx-1, depth+1, dfs);
    }

    if (root_idx < to) {
        dfs = init_reindex_subtree(
            idxs, root_idx+1, to, depth+1, dfs);
    }

    return dfs;
}


// reorder nodes in tree so that all nodes with depth <= pivot are on
// one side and those with depth > on the other (like quisort iterations), given
// temporary space `tmp` at least as large as `tree`.
// Returns the index of the last node in the left partition.
// This is basically stable_partition in c++ stl.
size_t depth_partition(reindex_s* tree, reindex_s* tmp, size_t n, uint32_t pivot)
{
    // count elemnts <= pivot
    size_t left_size = 0;
    for (size_t i = 0; i < n; ++i) {
        if (tree[i].depth <= pivot) ++left_size;
    }

    size_t l = 0, r = left_size;
    for (size_t i = 0; i < n; ++i) {
        if (tree[i].depth <= pivot) {
            tmp[l++] = tree[i];
        } else {
            tmp[r++] = tree[i];
        }
    }

    memcpy(tree, tmp, n*sizeof(reindex_s));

    return left_size;
}


// recursively reorder tree to put it in VeB order. Assumes tree has already
// been sorted on dfs number.
void veb_order_recursion(
        reindex_s* tree, reindex_s* tmp, size_t n,
        uint32_t min_depth, uint32_t max_depth)
{
    if (min_depth == max_depth) {
        assert(n == 1);
        return;
    }

    uint32_t pivot_depth = min_depth + (max_depth - min_depth)/2;
    size_t top_size = depth_partition(tree, tmp, n, pivot_depth);

    // recurse on top subtree
    veb_order_recursion(tree, tmp, top_size, min_depth, pivot_depth);

    // find and recurse on bottom subtrees
    uint32_t bottom_subtree_depth = pivot_depth+1;
    size_t i = top_size;
    while (i < n) {
        size_t j = i+1;
        while (j < n && tree[j].depth != bottom_subtree_depth) ++j;
        veb_order_recursion(tree + i, tmp, j-i, bottom_subtree_depth, max_depth);
        i = j;
    }
}


// compute veb order for a tree of size n
void veb_order(size_t n)
{
    // initialized indexes
    reindex_s* idxs_tmp = malloc(sizeof(reindex_s) * n);
    for (int i = 0; i < n; ++i) idxs_tmp[i].i = i;
    init_reindex_subtree(idxs_tmp, 0, n-1, 0, 0);

    // put indexes in dfs order
    reindex_s* idxs = malloc(sizeof(reindex_s) * n);
    for (int i = 0; i < n; ++i) idxs[idxs_tmp[i].dfs] = idxs_tmp[i];

    // compute veb order
    uint32_t max_depth = 0;
    for (int i = 0; i < n; ++i) {
        if (idxs[i].depth > max_depth) max_depth = idxs[i].depth;
    }
    veb_order_recursion(idxs, idxs_tmp, n, 0, max_depth);

    // now use this to reorder the actual tree
    //
    free(idxs);
}

