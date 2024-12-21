//
// Created by mainf on 2024/7/17.
//
#include <stdio.h>
#include "base/forest.h"

// 示例程序
int main() {
    Forest* forest = init_forest(10); // 初始化森林，假设有9个节点

    add_child(forest, 7, 0);
    add_child(forest, 7, 1); // 添加 1 作为 7 的孩子
    add_child(forest, 4, 2);
    add_child(forest, 4, 2);
    add_child(forest, 5, 3);
    add_child(forest, 6, 4);
    add_child(forest, 6, 5);
    add_child(forest, 8, 6);
    add_child(forest, 8, 7);
    add_child(forest, 9, 8);
    printf("Child count of node 4: %lld\n", get_child_count(forest, 9));
    printf("parent of node 9: %lld\n", get_parent(forest, 9));

    //postorder_traversal(forest,8);

    for (INDEX_TYPE i = 0; i < forest->node_count; i++) {
        if (forest->parent[i] == -1) { // 假设 -1 表示根节点
            //postorder_traversal_recursion(forest, i);
            //postorder_traversal(forest, i);
            //partitioning(forest,i,2,2);
            printf("\n");
        }
    }

    free_forest(forest); // 销毁森林

    return 0;
}