#include <layer_matrix.h>
#include <base/matrix.h>
//
// Created by mainf on 2024/10/15.
//

// void lower_solver_block(const L2Matrix l2, ELE_TYPE *x, const ELE_TYPE *b, const INDEX_TYPE n) {
//     for (INDEX_TYPE i = 0; i < n; ++i) {
//         x[i] = b[i];
//     }
//     for (INDEX_TYPE j = 0; j < n; ++j) {
//         if (x[j] == 0) continue;
//         for (INDEX_TYPE i = Lp[j]; i < Lp[j + 1]; i++) {
//             x[Li[i]] -= Lx[i] * x[j];
//         }
//     }
// }

int main() {
    INDEX_TYPE row_ptr[] = {0, 3, 5, 8, 9, 11}; // 每行非零元素的起始位置
    INDEX_TYPE col_indices[] = {0, 2, 4, 1, 3, 0, 2, 3, 1, 0, 2}; // 非零元素的列索引
    ELE_TYPE values[] = {
        1, 0, 2, 0, 3,
        0, 4, 0, 5, 0,
        6, 0, 7, 8, 0,
        0, 9, 0, 0, 0,
        10, 0, 11, 0, 0
    };
    ELE_TYPE csr_values[] = {
        1,    2,    3,
           4,    5,
        6,    7, 8,
           9,
        10,   11,
    };
    BlockMatrix m1 = {
        .dense_values = values,
        .row_pointers = row_ptr,
        .col_indices = col_indices,
        .nnz = 11
    };
    return EXIT_SUCCESS;
}
