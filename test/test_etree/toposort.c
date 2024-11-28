#include <stdio.h>
#include <stdlib.h>

#include "../../src/base/malloc.c"

#define INDEX_TYPE long long

void toposort_etree(INDEX_TYPE n, const INDEX_TYPE parent[],
                    const INDEX_TYPE values[], const INDEX_TYPE threshold) {
    INDEX_TYPE *big_stack = lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE *small_stack = lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE *kid_count = lu_calloc(n + 1, sizeof(INDEX_TYPE));
    INDEX_TYPE top_big = -1, top_small = -1; // 栈顶指针

    kid_count++; //为了方便处理根是-1
    for (INDEX_TYPE i = 0; i < n; kid_count[parent[i++]]++);

    // 将所有叶子节点入栈
    for (INDEX_TYPE i = n - 1; i >= 0; i--) {
        if (kid_count[i] == 0) { // 叶子节点
            if (values[i] > threshold) {
                big_stack[++top_big] = i; // 入大值栈
            } else {
                small_stack[++top_small] = i; // 入小值栈
            }
        }
    }

    // 处理栈中的节点
    while (top_big != -1 || top_small != -1) {
        // 从小值栈取出一个节点，若空则从大值栈取
        INDEX_TYPE node;
        if (top_small != -1) {
            node = small_stack[top_small--]; // 弹出小值栈
        } else {
            node = big_stack[top_big--]; // 弹出大值栈
        }
        printf("%lld ", node);
        // 获取父节点，并减少其子节点数
        INDEX_TYPE parent_node = parent[node];
        kid_count[parent_node]--;

        // 如果父节点的子节点数为0，检查其值并入栈
        if (kid_count[parent_node] == 0) {
            if (values[parent_node] > 100) {
                big_stack[++top_big] = parent_node; // 入大值栈
            } else {
                small_stack[++top_small] = parent_node; // 入小值栈
            }
        }
    }
    lu_free(small_stack);
    lu_free(big_stack);
    lu_free(kid_count - 1);
}


int main() {
    INDEX_TYPE parent[] = {2, 2, -1, 4, 6, 6, -1, 9, 9, -1, 11, 13, 13, 16, 15, 16, -1};
    INDEX_TYPE values[16] = {0};
    values[3] = 200;
    values[9] = 200;
    values[11] = 200;
    toposort_etree(17, parent, values, 100);
}
