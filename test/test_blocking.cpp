#include <gtest/gtest.h>
#include <cstdio>
#include "base/file.h"

TEST(BlockingTest, Test02GetSubmatricesNNZ) {
    SparseMatrix matrix;
    INDEX_TYPE row_ptr[] = {0, 4, 7, 10, 14, 18, 22, 26, 32, 37};  // 每行非零元素的起始位置
    INDEX_TYPE col_indices[] = {0, 6, 7, 8,
                            1, 3, 5,
                            2, 4, 7,
                            1, 3, 7, 8,
                            2, 4, 5, 7,
                            1, 4, 5, 8,
                            0, 6, 7, 8,
                            0, 2, 3, 4, 6, 7,
                            0, 3, 5, 6, 8};  // 非零元素的列索引
    ELE_TYPE values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                         15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                         26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37};  // 非零元素的值
    matrix.row_pointers=row_ptr;
    matrix.col_indices=col_indices;
    matrix.csr_values=values;
    matrix.num_row=9;
    matrix.num_col=9;
    matrix.nnz=37;

    INDEX_TYPE num_col_block = 3;
    INDEX_TYPE num_row_block = 3;

    INDEX_TYPE *block_nnz_arr = get_sub_matrices_nnz(&matrix, 3, 3);
    INDEX_TYPE expect_values[] = {3, 3, 4, 3, 5, 4, 4, 4, 7};
    for (INDEX_TYPE i = 0; i < num_row_block * num_col_block; ++i) {
        EXPECT_EQ(expect_values[i], block_nnz_arr[i]);
    }
}