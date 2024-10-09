//
// Created by mainf on 2024/9/17.
//
#include "base/matrix.h"
#include "preprocess.h"

void GETRF(SparseMatrix **m, INDEX_TYPE i) {
    // LU 分解的实现
    LOG_DEBUG("GETRF called");
}

void GESSM(SparseMatrix **m, INDEX_TYPE j, INDEX_TYPE k) {
    // 上三角解的实现
    LOG_DEBUG("GESSM called");
}

void TSTRF(SparseMatrix **m, INDEX_TYPE i, INDEX_TYPE k) {
    // 下三角解的实现
    LOG_DEBUG("TSTRF called");
}

void SSSSM(SparseMatrix **m, INDEX_TYPE i, INDEX_TYPE j, INDEX_TYPE k) {
    // 更新子矩阵的实现
    LOG_DEBUG("SSSSM called");
}

void factorization_parallel_v0(SparseMatrix *A) {
    INDEX_TYPE num_row_block = A->num_row_block;
    INDEX_TYPE num_col_block = A->num_col_block;
    INDEX_TYPE L_index[num_row_block];
    INDEX_TYPE U_index[num_col_block];
    for (INDEX_TYPE i = 0; i < A->num_col; ++i) {//枚举列
        GETRF(A->sub_matrices, i); // LU 分解
        //符号计算 L
        INDEX_TYPE L_index_length = 0;
        for (int64_t j = i; j < num_row_block; ++j) {
            (A->sub_matrices[j * num_col_block + i] != NULL) ? (L_index[L_index_length++] = j) : 0;
        }
        //符号计算 U
        INDEX_TYPE U_index_length = 0;
        for (int64_t j = i; j < num_col_block; ++j) {
            (A->sub_matrices[j * num_row_block + i] != NULL) ? (U_index[U_index_length++] = j) : 0;
        }
        #pragma omp parallel for
        for (INDEX_TYPE p = 0; p < U_index_length; p++) {
            GESSM(A->sub_matrices, U_index[p], i); // 上三角解
        }
        #pragma omp parallel for
        for (INDEX_TYPE k = 0; k < L_index_length; k++) {
            TSTRF(A->sub_matrices, i, L_index[k]); // 下三角解
            for (INDEX_TYPE p = 0; p < U_index_length; p++) {
                SSSSM(A->sub_matrices, i, L_index[k], U_index[p]); // 更新子矩阵
            }
        }
    }
}

void factorization_parallel_v1(SparseMatrix *A) {
    INDEX_TYPE num_row_block = A->num_row_block;
    INDEX_TYPE num_col_block = A->num_col_block;
    INDEX_TYPE L_index[num_row_block];
    INDEX_TYPE U_index[num_col_block];

    for (INDEX_TYPE i = 0; i < A->num_col; ++i) { // 枚举列
        #pragma omp task depend(inout: A->sub_matrices[i * num_col_block + i])
        {
            GETRF(A->sub_matrices, i); // LU 分解
        }

        // 符号计算 L
        INDEX_TYPE L_index_length = 0;
        for (int64_t j = i; j < num_row_block; ++j) {
            if (A->sub_matrices[j * num_col_block + i] != NULL) {
                L_index[L_index_length++] = j;
            }
        }

        // 符号计算 U
        INDEX_TYPE U_index_length = 0;
        for (int64_t j = i; j < num_col_block; ++j) {
            if (A->sub_matrices[j * num_row_block + i] != NULL) {
                U_index[U_index_length++] = j;
            }
        }

        // 任务：处理上三角部分 U
        for (INDEX_TYPE p = 0; p < U_index_length; p++) {
            #pragma omp task depend(inout: A->sub_matrices[i * num_col_block + U_index[p]]) depend(in: A->sub_matrices[i * num_col_block + i])
            {
                GESSM(A->sub_matrices, U_index[p], i); // 上三角解
            }
        }

        // 任务：处理下三角部分 L
        for (INDEX_TYPE k = 0; k < L_index_length; k++) {
            #pragma omp task depend(inout: A->sub_matrices[L_index[k] * num_col_block + i]) depend(in: A->sub_matrices[i * num_col_block + i])
            {
                TSTRF(A->sub_matrices, i, L_index[k]); // 下三角解
            }

            // 任务：更新子矩阵
            for (INDEX_TYPE p = 0; p < U_index_length; p++) {
                #pragma omp task depend(inout: A->sub_matrices[L_index[k] * num_col_block + U_index[p]]) \
                                  depend(in: A->sub_matrices[L_index[k] * num_col_block + i], A->sub_matrices[i * num_col_block + U_index[p]])
                {
                    SSSSM(A->sub_matrices, i, L_index[k], U_index[p]); // 更新子矩阵
                }
            }
        }
    }

    // 等待所有任务完成
    #pragma omp taskwait
}