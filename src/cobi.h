
#pragma once

#include <stddef.h>
#include <stdint.h>

typedef struct node_s
{
    int32_t first;
    int32_t last;

    int32_t subtree_first;
    int32_t subtree_last;

    uint32_t lindex;
    uint32_t rindex;
} node_s;


void cobi_init(node_s* nodes, size_t numnodes);

void veb_order(size_t n);

