//
// Created by mainf on 2025/3/8.
//

#include "../src/toposort.h"

#include <math.h>
#include <../src/symbolic_analysis.h>
#include <../src/base/base_math.h>
#include <../src/base/plot.h>

#include "../src/etree.h"

INDEX_TYPE *eq_order(INDEX_TYPE n, const INDEX_TYPE parent[], INDEX_TYPE max_size) {
    INDEX_TYPE *order = lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE *kid_count = lu_calloc(n + 1, sizeof(INDEX_TYPE));
    INDEX_TYPE *sub_tree_size = lu_calloc(n + 1, sizeof(INDEX_TYPE));
    INDEX_TYPE *stack = lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE top = 0; //head指向stack下一个元素要放的位置
    INDEX_TYPE rear = n; //rear表示栈底, head==rear是true表示队空
    INDEX_TYPE top_order = -1;
    for (INDEX_TYPE i = 0; i < n; sub_tree_size[i++] = 1);
    kid_count++; //为了方便处理根是-1
    for (INDEX_TYPE i = 0; i < n; kid_count[parent[i++]]++);
    // 将所有叶子节点入栈
    for (INDEX_TYPE i = n - 1; i >= 0; i--) {
        if (kid_count[i] == 0) {
            // 叶子节点
            stack[top++] = i;
        }
    }

    while (top != rear) {
        INDEX_TYPE node;
        node = stack[top = (top - 1) % n]; //dequeue
        //printf("%lld ", node);
        order[++top_order] = node;
        // 获取父节点，并减少其子节点数
        INDEX_TYPE parent_node = parent[node];
        kid_count[parent_node]--, sub_tree_size[parent_node] += sub_tree_size[node];
        // 如果父节点的子节点数为0，检查其值并入栈
        while (parent_node > 0 && kid_count[parent_node] == 0) {
            if (sub_tree_size[parent_node] > max_size) {
                stack[--rear] = parent_node;
                sub_tree_size[parent_node] = 1;
                break;
            }
            order[top_order++] = parent_node;
            parent_node = parent[parent_node];
            kid_count[parent_node]--, sub_tree_size[parent_node] += sub_tree_size[node];
        }
    }

    lu_free(stack);
    lu_free(kid_count);
    lu_free(kid_count - 1);
    lu_free(sub_tree_size);
    return order;
}
