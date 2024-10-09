//
// Created by mainf on 2024/9/27.
//
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

    PreprocessInfo info;
    INDEX_TYPE diag[] = {0, 4, 7, 11, 15, 20, 23, 31, 37};
    info.diag_index_csr = diag;
    create_row_etree_sym(&matrix, &info);
    EXPECT_EQ(info.row_etree->node_count, 9);
    INDEX_TYPE expect_parent[] = {6, 3, 4, 5, 5, 7, 7, 8, -1};
    for (INDEX_TYPE i = 0; i < info.row_etree->node_count; ++i) {
        EXPECT_EQ(info.row_etree->parent[i], expect_parent[i]);
    }
}
