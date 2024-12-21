//
// Created by mainf on 2024/12/2.
//
#include <iostream>
#include <queue>
#include <etree.h>
#include <base/matrix.h>
#include <symbolic_analysis.h>
#include "toposort.h"

void reverse_queue(std::queue<INDEX_TYPE> &q) {
    std::stack<INDEX_TYPE> s;

    // 将队列元素移动到栈中
    while (!q.empty()) {
        s.push(q.front());
        q.pop();
    }

    // 将栈中的元素再移动回队列中
    while (!s.empty()) {
        q.push(s.top());
        s.pop();
    }
}

///行宽：L的行第一个非零元的列号-对角线的列号
void calc_row_wide(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai, INDEX_TYPE *values,INDEX_TYPE n) {
    for (INDEX_TYPE i = 0; i < n; ++i) {
        INDEX_TYPE idx = Ai[Ap[i]];
        values[i] = i - idx;
    }
}

INDEX_TYPE calc_threshold(INDEX_TYPE n) {
    return (INDEX_TYPE) (sqrt(n));
}

/**
 * @param parent 父亲表示法的森林
 * @param top_order order大小
 * @param order 拓扑序列
 * @param kid_count 节点的度
 * @param values 行宽
 * @param threshold 行宽大小判定的阈值
 * @param small_queue 度为0的行宽小的队列
 * @param big_queue 大行宽队列
 */
void toposort(const INDEX_TYPE parent[],
              const INDEX_TYPE values[],
              const INDEX_TYPE threshold,
              INDEX_TYPE *top_order,
              INDEX_TYPE *order,
              INDEX_TYPE *kid_count,
              std::queue<INDEX_TYPE> &small_queue,
              std::queue<INDEX_TYPE> &big_queue) {
    while (!small_queue.empty()) {
        INDEX_TYPE node = small_queue.front();
        small_queue.pop();
        //std::cout << values[node] << " ";
        order[(*top_order)++] = node;
        // 获取父节点，并减少其子节点数
        INDEX_TYPE parent_node = parent[node];
        kid_count[parent_node]--;
        // 如果父节点的子节点数为0，检查其值并入队
        while (parent_node > 0 && kid_count[parent_node] == 0) {
            if (values[parent_node] > threshold) {
                big_queue.push(parent_node); // 加入大值
                break;
            }
            //std::cout << values[parent_node] << " ";
            order[(*top_order)++] = parent_node;
            parent_node = parent[parent_node];
            kid_count[parent_node]--;
        }
    }
}

/**
 * @param parent 父亲表示法的森林
 * @param top_order order大小
 * @param order 拓扑序列
 * @param kid_count 节点的度
 * @param queue 度为0的节点
 */
void toposort_2(const INDEX_TYPE parent[],
                INDEX_TYPE *top_order,
                INDEX_TYPE *order,
                INDEX_TYPE *kid_count,
                std::queue<INDEX_TYPE> &queue) {
    while (!queue.empty()) {
        INDEX_TYPE node = queue.front();
        queue.pop();
        order[(*top_order)++] = node;
        // 获取父节点，并减少其子节点数
        INDEX_TYPE parent_node = parent[node];
        kid_count[parent_node]--;
        // 如果父节点的子节点数为0
        while (parent_node > 0 && kid_count[parent_node] == 0) {
            order[(*top_order)++] = parent_node;
            parent_node = parent[parent_node];
            kid_count[parent_node]--;
        }
    }
}

INDEX_TYPE *toposort_help(INDEX_TYPE n,
                          const INDEX_TYPE parent[],
                          SparseMatrix *A,
                          INDEX_TYPE *cut_point) {
    INDEX_TYPE *kid_count_0 = (INDEX_TYPE *) lu_calloc(n + 1, sizeof(INDEX_TYPE));
    INDEX_TYPE *kid_count = kid_count_0 + 1;
    INDEX_TYPE *values = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE *order = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    calc_row_wide(A->row_pointers, A->col_indices, values, n);
    INDEX_TYPE threshold = calc_threshold(n);

    std::queue<INDEX_TYPE> small_queue;
    std::queue<INDEX_TYPE> big_queue;
    // 计算每个节点的子节点数
    for (INDEX_TYPE i = 0; i < n; kid_count[parent[i++]]++);
    // 将所有叶子节点入栈或入队
    for (INDEX_TYPE i = 0; i < n; i++) {
        if (kid_count[i] == 0) { // 叶子节点
            if (values[i] > threshold) {
                big_queue.push(i); // 加入大值队列
            } else {
                small_queue.push(i); // 入小值队列
            }
        }
    }
    INDEX_TYPE top_order = 0;
    toposort(parent, values, threshold, &top_order, order, kid_count, small_queue, big_queue);
    LOG_INFO("cut_point:%lld", top_order);
    *cut_point = top_order;
    //新的values
    INDEX_TYPE total_nnz = 0;
    toposort_2(parent, &top_order, order, kid_count, big_queue);
    lu_free(values);
    lu_free(kid_count_0);
    return order;
}


INDEX_TYPE *reorder_toposort(SparseMatrix *A, INDEX_TYPE n, INDEX_TYPE *cut_point) {
    clock_t start_time = clock();
    INDEX_TYPE *parent = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE *b_ptr, *b_idx;
    INDEX_TYPE bnz = A->nnz * 2;
    //todo 对角线添加
    a_plus_at(n, A->nnz, A->row_pointers, A->col_indices, &bnz, &b_ptr, &b_idx);
    LOG_DEBUG("a_plus_at nnz : %lld", bnz);
    LOG_DEBUG("a nnz : %lld", A->nnz-n);
    sp_symetree(b_ptr, b_ptr + 1, b_idx, n, parent);
    INDEX_TYPE threshold = calc_threshold(n);
    LOG_DEBUG("threshold: %lld", threshold);
    INDEX_TYPE *order = toposort_help(n, parent, A, cut_point);
    INDEX_TYPE *iorder = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    for (int i = 0; i < n; ++i) {
        iorder[order[i]] = i;
    }
    lu_free(order);
    lu_free(parent);
    LOG_INFO("reorder_toposort elapsed time: %lf ms", ((double) (clock() - start_time)) / CLOCKS_PER_SEC * 1000.0);
    return iorder;
}
