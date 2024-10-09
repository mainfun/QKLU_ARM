//
// Created by mainf on 2024/7/17.
//

#include <stdio.h>
#include <stdlib.h>
#include "base/log.h"
#include "forest.h"

// 初始化 Forest 结构体
Forest *init_forest(INDEX_TYPE node_count) {
    Forest *forest = (Forest *) lu_malloc(sizeof(Forest));

    forest->first_child = (INDEX_TYPE *) lu_calloc(node_count, sizeof(INDEX_TYPE));
    forest->next_sibling = (INDEX_TYPE *) lu_calloc(node_count, sizeof(INDEX_TYPE));
    forest->last_sibling = (INDEX_TYPE *) lu_calloc(node_count, sizeof(INDEX_TYPE)); // 分配尾指针的内存
    forest->parent = (INDEX_TYPE *) lu_calloc(node_count, sizeof(INDEX_TYPE));
    forest->child_count = (INDEX_TYPE *) lu_calloc(node_count, sizeof(INDEX_TYPE));
    forest->node_count = node_count;

    // 初始化所有节点的父节点和子节点计数
    for (INDEX_TYPE i = 0; i < node_count; i++) {
        forest->first_child[i] = -1; // -1 表示没有子节点
        forest->next_sibling[i] = -1; // -1 表示没有兄弟节点
        forest->last_sibling[i] = -1; // -1 表示没有兄弟节点
        forest->parent[i] = -1; // -1 表示没有父节点
        forest->child_count[i] = 0;
    }

    return forest;
}



// 销毁 Forest 结构体
void free_forest(Forest *forest) {
    if (forest != NULL) {
        lu_free(forest->first_child);
        lu_free(forest->next_sibling);
        lu_free(forest->last_sibling);
        lu_free(forest->parent);
        lu_free(forest->child_count);
        lu_free(forest);
    }
}

// 添加孩子节点
void add_child(Forest *forest, INDEX_TYPE parent, INDEX_TYPE child) {
    if (parent >= forest->node_count || child >= forest->node_count) {
        LOG_ERROR("Node index out of bounds\n");
    }
    if(forest->parent[child] == parent) return;
    forest->parent[child] = parent;
    forest->child_count[parent]++;

    if (forest->first_child[parent] == -1) {
        forest->first_child[parent] = child;
        forest->last_sibling[parent] = child;
    } else {
        INDEX_TYPE last_sibling = forest->last_sibling[parent];
        forest->next_sibling[last_sibling] = child;
        forest->last_sibling[parent] = child;
    }
}

// 获取父节点
INDEX_TYPE get_parent(Forest *forest, INDEX_TYPE node) {
    if (node >= forest->node_count) {
        LOG_ERROR("Node index out of bounds\n");
    }
    return forest->parent[node];
}

// 获取第一个孩子
INDEX_TYPE get_first_child(Forest *forest, INDEX_TYPE node) {
    if (node >= forest->node_count) {
        LOG_ERROR("Node index out of bounds\n");
    }
    return forest->first_child[node];
}

// 获取下一个兄弟节点
INDEX_TYPE get_next_sibling(Forest *forest, INDEX_TYPE node) {
    if (node >= forest->node_count) {
        LOG_ERROR("Node index out of bounds\n");
    }
    return forest->next_sibling[node];
}

// 获取孩子数量
INDEX_TYPE get_child_count(Forest *forest, INDEX_TYPE node) {
    if (node >= forest->node_count) {
        LOG_ERROR("Node index out of bounds\n");
    }
    return forest->child_count[node];
}

// 遍历节点的子节点
void traverse_children(Forest *forest, INDEX_TYPE node) {
    if (node >= forest->node_count) {
        LOG_ERROR("Node index out of bounds\n");
    }

    INDEX_TYPE child = forest->first_child[node];
    while (child != -1) {
        LOG_DEBUG("Child of node %lld: %lld", node, child);
        child = forest->next_sibling[child];
    }
}

//todo
void get_parallel_leaves(Forest *forest) {
    INDEX_TYPE *original_child_count = (INDEX_TYPE *) lu_malloc(forest->node_count * sizeof(INDEX_TYPE));
    if (original_child_count == NULL) {
        LOG_ERROR("Failed to allocate memory for original_child_count");
    }
    // 备份 child_count 数组
    for (INDEX_TYPE i = 0; i < forest->node_count; i++) {
        original_child_count[i] = forest->child_count[i];
    }

    // 用于存储叶子节点的队列
    INDEX_TYPE *queue = (INDEX_TYPE *) lu_malloc(forest->node_count * sizeof(INDEX_TYPE));
    INDEX_TYPE front = 0, rear = 0;

    // 初始化队列，找到所有初始叶子节点
    for (INDEX_TYPE i = 0; i < forest->node_count; i++) {
        if (forest->child_count[i] == 0 && forest->parent[i] != -1) {
            queue[rear++] = i;
        }
    }

    INDEX_TYPE total_nodes = forest->node_count;
    INDEX_TYPE prefix_sum_index = 0;
    //prefix_sum[prefix_sum_index++] = 0;
    INDEX_TYPE current_level_size = rear - front;
    while (current_level_size > 0) {
        INDEX_TYPE leaf_count = 0;
        LOG_DEBUG("\nLeaf node: ");
        for (INDEX_TYPE i = 0; i < current_level_size; i++) {
            INDEX_TYPE leaf = queue[front++];
            LOG_DEBUG("%lld, ", leaf);
            INDEX_TYPE parent = forest->parent[leaf];
            if (parent != -1) {//如果有父节点
                forest->child_count[parent]--;//度减一
                if (forest->child_count[parent] == 0) {
                    queue[rear++] = parent;
                }
            }
            leaf_count++;
            total_nodes--;
        }
        //prefix_sum[prefix_sum_index] = prefix_sum[prefix_sum_index - 1] + leaf_count;
        prefix_sum_index++;
        current_level_size = rear - front;
    }

    // 恢复 child_count 数组
    lu_free(forest->child_count);
    forest->child_count = original_child_count;
}