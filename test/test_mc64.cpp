//
// Created by mainf on 2024/9/27.
//
#include <pivoting.h>
#include <symbolic_analysis.h>
#include <gtest/gtest.h>
#include "base/forest.h"

TEST(TEST_ETREE, TEST_ETREE_CREATE) {
    SparseMatrix matrix;
    INDEX_TYPE row_ptr[] = {0, 4, 7, 10, 14, 18, 22, 26, 32, 37}; // 每行非零元素的起始位置
    INDEX_TYPE col_indices[] = {
        0, 6, 7, 8,
        1, 3, 5,
        2, 4, 7,
        1, 3, 7, 8,
        2, 4, 5, 7,
        1, 4, 5, 8,
        0, 6, 7, 8,
        0, 2, 3, 4, 6, 7,
        0, 3, 5, 6, 8
    }; // 非零元素的列索引
    ELE_TYPE values[] = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
        26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37
    }; // 非零元素的值
    matrix.row_pointers = row_ptr;
    matrix.col_indices = col_indices;
    matrix.csr_values = values;
    matrix.num_row = 9;
    matrix.num_col = 9;
    matrix.nnz = 37;

    INDEX_TYPE *perm = NULL;
    INDEX_TYPE *iperm = NULL;
    ELE_TYPE *row_scale = NULL;
    ELE_TYPE *col_scale = NULL;
    // static_pivoting(&matrix, &perm, &iperm, &row_scale, &col_scale);
    // tmp_matrix = init_csr_matrix(n, n, nnz);
    // apply_pivoting(matrix, tmp_matrix, perm, row_scale, col_scale);
    // info->mc64_perm = perm;
    // info->mc64_iperm = iperm;
    // info->Dr = row_scale;
    // info->Dc = col_scale;
}
