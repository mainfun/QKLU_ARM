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

ELE_TYPE *csr2dense_v2(const INDEX_TYPE *Ap_start, const INDEX_TYPE *Ap_end, const INDEX_TYPE *Ai, const ELE_TYPE *Ax,
                       INDEX_TYPE n) {
    ELE_TYPE *D = (ELE_TYPE *) calloc(n * n, sizeof(ELE_TYPE));
    for (INDEX_TYPE r = 0; r < n; ++r) {
        for (INDEX_TYPE j = Ap_start[r]; j < Ap_end[r]; j++) {
            D[r * n + Ai[j]] = Ax[j];
        }
    }
    return D;
}

ELE_TYPE *csc2dense_v2(const INDEX_TYPE *Ap_start, const INDEX_TYPE *Ap_end, const INDEX_TYPE *Ai, const ELE_TYPE *Ax,
                       INDEX_TYPE n) {
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
