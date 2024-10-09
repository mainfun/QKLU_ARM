#include <math.h>

#include "numerical.h"
#include "base/matrix.h"
#include "preprocess.h"

#define THRESHOLD 1e-8

/**
 * LU分解的朴素实现(48ms~65ms) CSC
 */
void factor(const SparseMatrix *A, const PreprocessInfo *A_pattern) {
    INDEX_TYPE n = A->num_row;
    INDEX_TYPE elimination_count = 0;
    //稀疏转稠密
    ELE_TYPE *D;
    D = (ELE_TYPE *) lu_calloc(n * n, sizeof(ELE_TYPE));
    //求对角元素索引
    //包括填充元的矩阵的对角位置
    INDEX_TYPE *Ap = A->row_pointers;
    INDEX_TYPE *Ai = A->col_indices;
    ELE_TYPE *Ax = A->csr_values;
    //填充
    for (INDEX_TYPE r = 0; r < n; ++r) {
        ELE_TYPE *p = D + r * n; //取行指针
        for (INDEX_TYPE j = Ap[r]; j < Ap[r + 1]; j++) {
            //            __builtin_prefetch(p + Ai[j + 8], 1, 3);
            p[Ai[j]] = Ax[j];
        }
    }
    LOG_INFO("开始消元:");
    clock_t factor_time = clock();
    //向下高斯消元
    for (INDEX_TYPE i = 0; i < n; ++i) {
        //枚举列
        ELE_TYPE pivot = D[i * n + i]; // diag value
        if (fabs(pivot) < THRESHOLD) pivot = THRESHOLD;
        D[i * n + i] = pivot;
        //#pragma omp parallel for num_threads(6)
        for (INDEX_TYPE p = A_pattern->L->col_pointers[i]; p < A_pattern->L->col_pointers[i + 1]; p++) {
            INDEX_TYPE j = A_pattern->L->row_indices[p]; //行号
            //printf("第%lld行消第%lld行\n", i + 1, j + 1);
            ELE_TYPE scale = D[j * n + i] / pivot;
            //L的列
            D[j * n + i] = scale;
            ELE_TYPE *pivot_row_ptr = D + i * n;
            ELE_TYPE *eli_row_ptr = D + j * n;
            //优化访存
            for (INDEX_TYPE k = A_pattern->U->row_pointers[i]; k < A_pattern->U->row_pointers[i + 1]; k++) {
                INDEX_TYPE c = A_pattern->U->col_indices[k];
                eli_row_ptr[c] -= scale * pivot_row_ptr[c];
                elimination_count++;
            }
        }
        //printf("elimination_count:%lld\n", elimination_count);
    }
    LOG_INFO("LU factor elapsed time: %lf ms", ((double) (clock() - factor_time)) / CLOCKS_PER_SEC * 1000.0); \
    //写回
    SparseMatrix *L = A_pattern->L;
    SparseMatrix *U = A_pattern->U;
    // L->csc_values = lu_malloc(L->nnz * sizeof(ELE_TYPE));
    // U->csr_values = lu_malloc(U->nnz * sizeof(ELE_TYPE));
    // U->csc_values = lu_malloc(U->nnz * sizeof(ELE_TYPE));
    LOG_DEBUG("lnz=%lld,unz=%lld\n", L->nnz, U->nnz);
    INDEX_TYPE l_count_csc = 0;
    INDEX_TYPE l_count_csr = 0;
    INDEX_TYPE u_count_csr = 0;
    INDEX_TYPE u_count_csc = 0;
    L->num_row = L->num_col = U->num_row = U->num_col = n;
    for (INDEX_TYPE i = 0; i < n; i++) {
        //U csr
        for (INDEX_TYPE j = A_pattern->U->row_pointers[i]; j < A_pattern->U->row_pointers[i + 1]; j++) {
            INDEX_TYPE index = A_pattern->U->col_indices[j];
            U->csr_values[u_count_csr++] = D[i * n + index];
            //printf("D[%lld][%lld] ", i, index);
        }
        //U csc
        // for (LU_INT j = A_pattern->U->col_pointers[i]; j < A_pattern->U->col_pointers[i + 1]; j++) {
        //     LU_INT index = A_pattern->U->row_indices[j];
        //     U->csc_values[u_count_csc++] = D[index * n + i];
        // }
        //L csc
        for (INDEX_TYPE j = A_pattern->L->col_pointers[i]; j < A_pattern->L->col_pointers[i + 1]; j++) {
            INDEX_TYPE index = A_pattern->L->row_indices[j];
            L->csc_values[l_count_csc++] = D[index * n + i];
        }
        //L csr
        for (INDEX_TYPE j = A_pattern->L->row_pointers[i]; j < A_pattern->L->row_pointers[i + 1]; j++) {
            INDEX_TYPE index = A_pattern->L->col_indices[j];
            L->csc_values[l_count_csr++] = D[i * n + index];
        }
    }
    LOG_DEBUG("消元次数为   ::::%lld\n", elimination_count);
    if (A->num_row < 10) {
        printf("U:\n");
        print_matrix_csr(U->row_pointers, U->col_indices, n);
        printf("D:\n");
        printf("L:\n");
        print_matrix_csr(L->row_pointers, L->col_indices, n);
        printf("D:\n");
    }
    //my_lu free
    lu_free(D);
}
