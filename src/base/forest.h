//
// Created by mainf on 2024/7/17.
//

#ifndef QKLU_FOREST_H
#define QKLU_FOREST_H
#ifdef __cplusplus
extern "C" {
#endif
#include "matrix.h"

// Forest 结构体定义
typedef struct {
    INDEX_TYPE *first_child;   // 指向第一个孩子的指针数组
    INDEX_TYPE *next_sibling;  // 指向下一个兄弟节点的指针数组
    INDEX_TYPE *last_sibling;  // 指向最后一个兄弟节点的指针数组
    INDEX_TYPE *parent;        // 指向父节点的指针数组
    INDEX_TYPE *child_count;   // 存储每个节点的孩子数量的数组
    INDEX_TYPE node_count;     // 节点数量
} Forest;

/**
 * 初始化 Forest 结构体
 * @param node_count 节点数量
 * @return 返回初始化后的 Forest 指针
 */
Forest *init_forest(INDEX_TYPE node_count);

/**
 * 销毁 Forest 结构体，释放内存
 * @param forest 指向 Forest 结构体的指针
 */
void free_forest(Forest *forest);

/**
 * 添加孩子节点
 * @param forest 指向 Forest 结构体的指针
 * @param parent 父节点的索引
 * @param child 要添加的孩子节点的索引
 */
void add_child(Forest *forest, INDEX_TYPE parent, INDEX_TYPE child);

/**
 * 获取父节点
 * @param forest 指向 Forest 结构体的指针
 * @param node 节点的索引
 * @return 返回父节点的索引，如果节点不存在则返回 -1
 */
INDEX_TYPE get_parent(Forest *forest, INDEX_TYPE node);

/**
 * 获取第一个孩子
 * @param forest 指向 Forest 结构体的指针
 * @param node 节点的索引
 * @return 返回第一个孩子的索引，如果没有孩子则返回 -1
 */
INDEX_TYPE get_first_child(Forest *forest, INDEX_TYPE node);

/**
 * 获取下一个兄弟节点
 * @param forest 指向 Forest 结构体的指针
 * @param node 节点的索引
 * @return 返回下一个兄弟节点的索引，如果没有兄弟节点则返回 -1
 */
INDEX_TYPE get_next_sibling(Forest *forest, INDEX_TYPE node);

/**
 * 获取孩子数量
 * @param forest 指向 Forest 结构体的指针
 * @param node 节点的索引
 * @return 返回孩子数量
 */
INDEX_TYPE get_child_count(Forest *forest, INDEX_TYPE node);

/**
 * 遍历节点的子节点
 * @param forest 指向 Forest 结构体的指针
 * @param node 节点的索引
 */
void traverse_children(Forest *forest, INDEX_TYPE node);

/**
 * 输出计算序列
 * @param forest
 */
void get_parallel_leaves(Forest *forest);

#ifdef __cplusplus
}
#endif
#endif //QKLU_FOREST_H