//
// Created by mainf on 2024/9/15.
//

#include "matrix.h"

SparseMatrix *init_sparse_matrix(INDEX_TYPE num_row, INDEX_TYPE num_col, INDEX_TYPE nnz) {
    SparseMatrix *matrix = (SparseMatrix *) lu_malloc(sizeof(SparseMatrix));
    matrix->num_row = num_row;
    matrix->num_col = num_col;
    matrix->nnz = nnz;
    matrix->row_pointers = NULL;
    matrix->col_indices = NULL;
    matrix->csr_values = NULL;

    matrix->col_pointers = NULL;
    matrix->row_indices = NULL;
    matrix->csc_values = NULL;
    matrix->sub_matrices = NULL;

    matrix->num_row_block = 0;
    matrix->num_col_block = 0;
    return matrix;
}

SparseMatrix *init_csr_matrix(INDEX_TYPE num_row, INDEX_TYPE num_col, INDEX_TYPE nnz) {
    SparseMatrix *matrix = init_sparse_matrix(num_row, num_col, nnz);

    matrix->row_pointers = (INDEX_TYPE *) lu_malloc((num_row + 1) * sizeof(INDEX_TYPE));
    matrix->col_indices = (INDEX_TYPE *) lu_malloc(nnz * sizeof(INDEX_TYPE));
    matrix->csr_values = (ELE_TYPE *) lu_malloc(nnz * sizeof(ELE_TYPE));

    return matrix;
}

SparseMatrix *init_csc_matrix(INDEX_TYPE num_row, INDEX_TYPE num_col, INDEX_TYPE nnz) {
    SparseMatrix *matrix = init_sparse_matrix(num_row, num_col, nnz);
    matrix->col_pointers = (INDEX_TYPE *) lu_malloc((num_col + 1) * sizeof(INDEX_TYPE));
    matrix->row_indices = (INDEX_TYPE *) lu_malloc(nnz * sizeof(INDEX_TYPE));
    matrix->csc_values = (ELE_TYPE *) lu_malloc(nnz * sizeof(ELE_TYPE));
    return matrix;
}

void free_sparse_matrix(SparseMatrix *matrix) {
    if (matrix == NULL) return;
    lu_free(matrix->csr_values);
    lu_free(matrix->csc_values);
    if (matrix->row_pointers != matrix->col_pointers && matrix->row_indices != matrix->col_indices) {
        //CSR
        lu_free(matrix->row_pointers);
        lu_free(matrix->col_indices);
        //CSC
        lu_free(matrix->col_pointers);
        lu_free(matrix->row_indices);
    } else {
        lu_free(matrix->row_pointers);
        lu_free(matrix->col_indices);
    }

    if (matrix->sub_matrices != NULL) {
        for (INDEX_TYPE i = 0; i < matrix->num_col_block * matrix->num_row_block; i++) {
            //printf("free sub_matrices:%lld\n",i);
            free_sparse_matrix(matrix->sub_matrices[i]);
        }
        lu_free(matrix->sub_matrices);
    }
    lu_free(matrix);
}

/**
 * @param a >0
 * @param b >0
 */
INDEX_TYPE ceil_div(INDEX_TYPE a, INDEX_TYPE b) {
    return (a + b - 1) / b;
}

INDEX_TYPE *get_sub_matrices_nnz(SparseMatrix *m, INDEX_TYPE block_width, INDEX_TYPE block_height) {
    INDEX_TYPE *Ap = m->row_pointers;
    INDEX_TYPE *Ai = m->col_indices;

    INDEX_TYPE num_row_block = ceil_div(m->num_row, block_height);
    INDEX_TYPE num_col_block = ceil_div(m->num_col, block_width);
    LOG_DEBUG("num_row_block: %lld, num_col_block: %lld", num_row_block, num_col_block);

    INDEX_TYPE *block_nnz_arr = (INDEX_TYPE *) calloc(num_col_block * num_row_block, sizeof(INDEX_TYPE));
    for (INDEX_TYPE i = 0; i < m->num_row; ++i) {
        //i / block_height 是块的行号，Ai[j] / block_width是块的列号
        INDEX_TYPE block_row_ptr = (i / block_height) * num_col_block;
        for (INDEX_TYPE j = Ap[i]; j < Ap[i + 1]; ++j) {
            block_nnz_arr[block_row_ptr + Ai[j] / block_width]++;
        }
    }
    int count = 0;
    for (INDEX_TYPE i = 0; i < num_row_block * num_col_block; ++i) {
        if (block_nnz_arr[i] > 0) {
            count++;
        }
    }
    LOG_INFO("count=%d ", count);
    //    LOG_DEBUG("block nnz: \n");
    //    for (LU_INT i = 0; i < num_row_block; ++i) {
    //        for (LU_INT j = 0; j < num_col_block; ++j) {
    //            printf("%4lld ", block_nnz_arr[i * num_col_block + j]);
    //        }
    //        printf("\n");
    //    }
    return block_nnz_arr;
}


void blocking_csr_matrix(SparseMatrix *matrix, INDEX_TYPE block_width, INDEX_TYPE block_height) {
    INDEX_TYPE *block_nnz_arr = get_sub_matrices_nnz(matrix, block_width, block_height);

    INDEX_TYPE num_row_block = ceil_div(matrix->num_row, block_height);
    INDEX_TYPE num_row_remainder = matrix->num_row % block_height;

    INDEX_TYPE num_col_block = ceil_div(matrix->num_col, block_width);
    INDEX_TYPE num_col_remainder = matrix->num_col % block_width;

    /**---------------分配空间---------------**/
    INDEX_TYPE n = num_col_block * num_row_block;
    SparseMatrix **sub_matrices = (SparseMatrix **) lu_malloc(n * sizeof(SparseMatrix *));
    num_row_remainder = num_row_remainder ? num_row_remainder : block_height;
    num_col_remainder = num_col_remainder ? num_col_remainder : block_width;
    INDEX_TYPE index = 0;
    for (INDEX_TYPE i = 0; i < num_row_block; i++) {
        for (INDEX_TYPE j = 0; j < num_col_block; ++j) {
            if (block_nnz_arr[index] > 0) {
                INDEX_TYPE csr_num_row = i == num_row_block - 1 ? num_row_remainder : block_height;
                INDEX_TYPE csr_num_col = j == num_col_block - 1 ? num_col_remainder : block_width;
                sub_matrices[index] = init_csr_matrix(csr_num_row, csr_num_col, block_nnz_arr[index]);
                sub_matrices[index]->nnz = 0;
            } else {
                sub_matrices[index] = NULL;
            }
            index++;
        }
    }
    /**--------------end分配空间---------------**/
    // 遍历原始矩阵并填充子块
    for (INDEX_TYPE i = 0; i < matrix->num_row; i++) {
        INDEX_TYPE row_block_idx = i / block_height;
        for (INDEX_TYPE row_nz_idx = matrix->row_pointers[i]; row_nz_idx < matrix->row_pointers[i + 1]; row_nz_idx++) {
            INDEX_TYPE col_idx = matrix->col_indices[row_nz_idx];
            INDEX_TYPE col_block_idx = col_idx / block_width;

            // 获取该元素对应的子块矩阵
            SparseMatrix *sub_matrix = sub_matrices[row_block_idx * num_col_block + col_block_idx];

            // 计算该子块中的行和列相对位置
            INDEX_TYPE local_row = i % block_height;
            INDEX_TYPE local_col = col_idx % block_width;

            // 添加元素到子块中（这里需要根据具体实现调整）
            sub_matrix->row_pointers[local_row + 1]++;
            sub_matrix->col_indices[sub_matrix->nnz] = local_col;
            sub_matrix->csr_values[sub_matrix->nnz] = matrix->csr_values[row_nz_idx];
            sub_matrix->nnz++;
        }
    }

    // 更新每个子块的行指针数组
    for (INDEX_TYPE i = 0; i < num_row_block; i++) {
        for (INDEX_TYPE j = 0; j < num_col_block; j++) {
            SparseMatrix *sub_matrix = sub_matrices[i * num_col_block + j];
            if (sub_matrix != NULL) {
                for (INDEX_TYPE r = 1; r <= sub_matrix->num_row; r++) {
                    sub_matrix->row_pointers[r] += sub_matrix->row_pointers[r - 1];
                }
            }
        }
    }
    lu_free(block_nnz_arr);

    matrix->num_row_block = num_row_block;
    matrix->num_col_block = num_col_block;
    matrix->block_width = block_width;
    matrix->block_height = block_height;
    matrix->sub_matrices = sub_matrices;
}

void print_matrix_csr(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai, INDEX_TYPE n) {
    if (n > 20) {
        for (INDEX_TYPE r = 0; r < 100; ++r) {
            for (INDEX_TYPE _j = Ap[r]; _j < Ap[r + 1]; ++_j) {
                INDEX_TYPE c = Ai[_j];
                printf("%lld, %lld\n", r, c);
            }
        }
        return;
    }
    INDEX_TYPE i, j, p, found;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            found = 0;
            for (p = Ap[i]; p < Ap[i + 1]; p++) {
                INDEX_TYPE col_index = Ai[p];
                if (col_index == j) {
                    printf(" 1 ");
                    found = 1;
                    break;
                }
            }
            if (!found) {
                printf(" 0 ");
            }
        }
        printf("\n");
    }
}

ELE_TYPE *csc2dense(const SparseMatrix *A) {
    const INDEX_TYPE *Ap = A->col_pointers;
    const INDEX_TYPE *Ai = A->row_indices;
    const ELE_TYPE *Ax = A->csc_values;
    INDEX_TYPE n = A->num_row;
    ELE_TYPE *D;
    D = (ELE_TYPE *) calloc(n * n, sizeof(ELE_TYPE));
    for (INDEX_TYPE c = 0; c < n; ++c) {
        for (INDEX_TYPE j = Ap[c]; j < Ap[c + 1]; j++) {
            D[Ai[j] * n + c] = Ax == NULL ? 1 : Ax[j];
        }
    }
    return D;
}

ELE_TYPE *csr2dense(const SparseMatrix *A) {
    const INDEX_TYPE *Ap = A->row_pointers;
    const INDEX_TYPE *Ai = A->col_indices;
    const ELE_TYPE *Ax = A->csr_values;
    INDEX_TYPE n = A->num_row;

    ELE_TYPE *D;
    D = (ELE_TYPE *) calloc(n * n, sizeof(ELE_TYPE));
    for (INDEX_TYPE r = 0; r < n; ++r) {
        for (INDEX_TYPE j = Ap[r]; j < Ap[r + 1]; j++) {
            D[r * n + Ai[j]] = Ax == NULL ? 1 : Ax[j];
        }
    }
    return D;
}

ELE_TYPE *csr2dense_v2(const INDEX_TYPE *Ap_start, const INDEX_TYPE *Ap_end, const INDEX_TYPE *Ai, const ELE_TYPE *Ax,INDEX_TYPE n) {
    ELE_TYPE *D = (ELE_TYPE *) calloc(n * n, sizeof(ELE_TYPE));
    for (INDEX_TYPE r = 0; r < n; ++r) {
        for (INDEX_TYPE j = Ap_start[r]; j < Ap_end[r]; j++) {
            D[r * n + Ai[j]] = Ax[j];
        }
    }
    return D;
}

ELE_TYPE *csc2dense_v2(const INDEX_TYPE *Ap_start, const INDEX_TYPE *Ap_end, const INDEX_TYPE *Ai, const ELE_TYPE *Ax,INDEX_TYPE n) {
    ELE_TYPE *D = (ELE_TYPE *) calloc(n * n, sizeof(ELE_TYPE));
    for (INDEX_TYPE c = 0; c < n; ++c) {
        for (INDEX_TYPE j = Ap_start[c]; j < Ap_end[c]; j++) {
            D[Ai[j] * n + c] = Ax[j];
        }
    }
    return D;
}

void print_dense_matrix(ELE_TYPE *A, INDEX_TYPE n) {
    for (INDEX_TYPE i = 0; i < n; ++i) {
        for (INDEX_TYPE i = 0; i < n; ++i) {
            printf("---------------");
        }
        printf("\n|");
        for (INDEX_TYPE j = 0; j < n; ++j) {
            if (A[i * n + j] == 0) printf("              |");
            else printf("%13lf |", A[i * n + j]);
        }
        printf("\n");
    }
    for (INDEX_TYPE i = 0; i < n; ++i) {
        printf("---------------");
    }
    printf("\n");
}

void print_matrix_csc(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai, INDEX_TYPE n) {
    ELE_TYPE *A = (ELE_TYPE *) lu_calloc(n * n, sizeof(ELE_TYPE));
    for (INDEX_TYPE c = 0; c < n; ++c) {
        for (INDEX_TYPE j = Ap[c]; j < Ap[c + 1]; j++) {
            A[Ai[j] * n + c] = 1;
        }
    }
    for (INDEX_TYPE i = 0; i < n; ++i) {
        for (INDEX_TYPE j = 0; j < n; ++j) {
            printf("%3d", A[i * n + j] == 0 ? 0 : 1);
        }
        printf("\n");
    }
    lu_free(A);
}

void SpMV_csr(const SparseMatrix *A, const ELE_TYPE *x, ELE_TYPE *y) {
    for (INDEX_TYPE i = 0; i < A->num_row; i++) {
        y[i] = 0.0;
        for (INDEX_TYPE j = A->row_pointers[i]; j < A->row_pointers[i + 1]; j++) {
            y[i] += A->csr_values[j] * x[A->col_indices[j]];
        }
    }
}
