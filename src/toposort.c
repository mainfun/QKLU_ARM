//
// Created by mainf on 2024/11/28.
//

#include "toposort.h"

#include <math.h>
#include <symbolic_analysis.h>
#include <base/base_math.h>
#include <base/plot.h>

#include "etree.h"

INDEX_TYPE *toposort_by_etree(INDEX_TYPE n, const INDEX_TYPE parent[],
                              const INDEX_TYPE values[], const INDEX_TYPE threshold,
                              INDEX_TYPE *cut_point) {
    INDEX_TYPE *order = lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE *big_stack = lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE *small_stack = lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE *kid_count = lu_calloc(n + 1, sizeof(INDEX_TYPE));
    INDEX_TYPE top_big = -1, top_small = -1, top_order = -1; // 栈顶指针

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
    int is_first_pop_big_stack=true;
    // 处理栈中的节点
    while (top_big != -1 || top_small != -1) {
        // 从小值栈取出一个节点，若空则从大值栈取
        INDEX_TYPE node;
        if (top_small != -1) {
            node = small_stack[top_small--]; // 弹出小值栈
        } else {
            node = big_stack[top_big--]; // 弹出大值栈
            if (is_first_pop_big_stack) {
                *cut_point=top_order;
                is_first_pop_big_stack=false;
            }
        }
        //printf("%lld ", node);
        order[++top_order] = node;
        // 获取父节点，并减少其子节点数
        INDEX_TYPE parent_node = parent[node];
        kid_count[parent_node]--;

        // 如果父节点的子节点数为0，检查其值并入栈
        if (parent_node > 0 && kid_count[parent_node] == 0) {
            if (values[parent_node] > threshold) {
                big_stack[++top_big] = parent_node; // 入大值栈
            } else {
                small_stack[++top_small] = parent_node; // 入小值栈
            }
        }
    }
    lu_free(small_stack);
    lu_free(big_stack);
    lu_free(kid_count - 1);
    return order;
}

void calc_row_wide(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai, INDEX_TYPE *values,INDEX_TYPE n) {
    for (INDEX_TYPE i = 0; i < n; ++i) {
        INDEX_TYPE idx = Ai[Ap[i]];
        values[i] = i - idx;
    }
}

INDEX_TYPE calc_threshold(INDEX_TYPE n) {
    return (INDEX_TYPE)(sqrt(n));
}

INDEX_TYPE *reorder_toposort(SparseMatrix *A, INDEX_TYPE n, INDEX_TYPE *cut_point) {
    clock_t start_time = clock();
    INDEX_TYPE *parent = lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE *b_ptr, *b_idx;
    INDEX_TYPE bnz = A->nnz * 2;
    //todo 对角线添加
    a_plus_at(n, A->nnz, A->row_pointers, A->col_indices, &bnz, &b_ptr, &b_idx);
    //csr2image_block(b_ptr, b_idx,"arc130_aat.jpg", 0, 130, 0, 130);
    LOG_DEBUG("a_plus_at nnz : %lld",bnz);
    LOG_DEBUG("a nnz : %lld",A->nnz-n);
    sp_symetree(b_ptr, b_ptr + 1, b_idx, n, parent);
    INDEX_TYPE *values = lu_malloc(n * sizeof(INDEX_TYPE));
    calc_row_wide(A->row_pointers, A->col_indices, values, n);
    INDEX_TYPE threshold = calc_threshold(n);
    LOG_DEBUG("threshold: %lld", threshold);
    INDEX_TYPE *order = toposort_by_etree(n, parent, values, threshold, cut_point);
    LOG_DEBUG("cut_point: %lld", *cut_point);
    INDEX_TYPE *iorder = lu_malloc(n * sizeof(INDEX_TYPE));
    for (int i = 0; i < n; ++i) {
        iorder[order[i]] = i;
    }
    lu_free(values);
    lu_free(order);
    lu_free(parent);
    LOG_INFO("reorder_toposort elapsed time: %lf ms", ((double) (clock() - start_time)) / CLOCKS_PER_SEC * 1000.0);
    return iorder;
}
