//
// Created by mainf on 2024/9/16.
//
#include <gtest/gtest.h>
#include "base/matrix.h"
#include "base/file.h"

TEST(SparseMatrixTest, LoadMatrix) {
    SparseMatrix *m = load_matrix_csr("../res/test02.mtx", false);
    ASSERT_NE(m, nullptr);
    // 预期的数组内容
    INDEX_TYPE expected_row_ptr[] = {0, 4, 7, 10, 14, 18, 22, 26, 32, 37};  // 每行非零元素的起始位置
    INDEX_TYPE expected_col_index[] = {0, 6, 7, 8,
                                   1, 3, 5,
                                   2, 4, 7,
                                   1, 3, 7, 8,
                                   2, 4, 5, 7,
                                   1, 4, 5, 8,
                                   0, 6, 7, 8,
                                   0, 2, 3, 4, 6, 7,
                                   0, 3, 5, 6, 8};  // 非零元素的列索引
    ELE_TYPE expected_values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                                  15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                                  26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37};  // 非零元素的值
    EXPECT_EQ(m->num_row, 9);
    EXPECT_EQ(m->num_col, 9);
    EXPECT_EQ(m->num_row_block, 0);
    EXPECT_EQ(m->num_col_block, 0);
    EXPECT_EQ(m->nnz, 37);
    EXPECT_EQ(m->csc_values, nullptr);
    EXPECT_EQ(m->col_pointers, nullptr);
    EXPECT_EQ(m->row_indices, nullptr);
    EXPECT_EQ(m->sub_matrices, nullptr);
    // 验证 value 数组
    for (INDEX_TYPE i = 0; i < m->nnz; ++i) {
        EXPECT_DOUBLE_EQ(m->csr_values[i], expected_values[i]);
    }

    // 验证 row_ptr 数组
    for (size_t i = 0; i < m->num_row; ++i) {
        EXPECT_EQ(m->row_pointers[i], expected_row_ptr[i]);
    }

    // 验证 col_index 数组
    for (size_t i = 0; i < m->nnz; ++i) {
        EXPECT_EQ(m->col_indices[i], expected_col_index[i]);
    }

    ASSERT_NE(m, nullptr);  // 确保 m 不是 nullptr
}