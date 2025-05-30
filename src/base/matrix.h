//
// Created by mainf on 2024/9/15.
//
#ifndef QKLU_MATRIX_H
#define QKLU_MATRIX_H
#ifdef __cplusplus
extern "C" {
#endif
#include <stdbool.h>
#include <stdlib.h>
#include "malloc.h"
#include "log.h"

#define ELE_TYPE double
#define INDEX_TYPE long long

typedef struct SPARSE_MATRIX {
    INDEX_TYPE num_row; // 矩阵的行数量
    INDEX_TYPE num_col; // 矩阵的列数量
    INDEX_TYPE nnz; // 非零元素的总数
    /**CSR**/
    INDEX_TYPE *row_pointers; // 每行的起始索引
    INDEX_TYPE *col_indices; // 非零元素对应的列索引
    ELE_TYPE *csr_values; // 矩阵非零元素的值
    /**CSC**/
    INDEX_TYPE *col_pointers;
    INDEX_TYPE *row_indices;
    ELE_TYPE *csc_values; // 矩阵非零元素的值
    /**子矩阵块**/
    struct SPARSE_MATRIX **sub_matrices;
    INDEX_TYPE num_row_block;
    INDEX_TYPE num_col_block;
    INDEX_TYPE block_width;
    INDEX_TYPE block_height;
} SparseMatrix;

/**
 * 初始化
 * @param num_row 行个数=列大小=矩阵高
 * @param num_col 列个数=行大小=矩阵宽
 * @param nnz 非零元个数
 * @return SparseMatrix
 */
SparseMatrix *init_sparse_matrix(INDEX_TYPE num_row, INDEX_TYPE num_col, INDEX_TYPE nnz);

/**
 * 初始化
 * @param num_row 行个数=列大小=矩阵高
 * @param num_col 列个数=行大小=矩阵宽
 * @param nnz 非零元个数
 * @return SparseMatrix
 */
SparseMatrix *init_csr_matrix(INDEX_TYPE num_row, INDEX_TYPE num_col, INDEX_TYPE nnz);

/**
 * 初始化
 * @param num_row 行个数=列大小=矩阵高
 * @param num_col 列个数=行大小=矩阵宽
 * @param nnz 非零元个数
 * @return SparseMatrix
 */
SparseMatrix *init_csc_matrix(INDEX_TYPE num_row, INDEX_TYPE num_col, INDEX_TYPE nnz);

/**
 * 递归的释放本矩阵和它的全部分块矩阵的空间
 * @param matrix 稀疏矩阵
 */
void free_sparse_matrix(SparseMatrix *matrix);

void print_matrix_csr(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai, INDEX_TYPE n);

void print_matrix_csc(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai, INDEX_TYPE n);

ELE_TYPE *csc2dense_v2(const INDEX_TYPE *Ap_start, const INDEX_TYPE *Ap_end,
                       const INDEX_TYPE *Ai, const ELE_TYPE *Ax, INDEX_TYPE n);

ELE_TYPE *csr2dense_v2(const INDEX_TYPE *Ap_start, const INDEX_TYPE *Ap_end,
                       const INDEX_TYPE *Ai, const ELE_TYPE *Ax, INDEX_TYPE n);

void print_dense_matrix(ELE_TYPE *A, INDEX_TYPE n);

ELE_TYPE *csc2dense(const SparseMatrix *A);

ELE_TYPE *csr2dense(const SparseMatrix *A);

/**
 * A是CSR格式，Ax=y
 * @param A SparseMatrix
 * @param x ELE_TYPE
 * @param y ELE_TYPE
 */
void SpMV_csr(const SparseMatrix *A, const ELE_TYPE *x, ELE_TYPE *y);

/**
* @brief calc A'+A-E.
* @param[in]  n A n
* @param[in]  nz A nnz
* @param[in]  col_ptr A column pointer
* @param[in]  row_idx A row indices
* @param[out] bnz A'+A nnz
* @param[out] b_col_ptr A'+A column pointer
* @param[out] b_row_ind A'+A row indices
**/
void a_plus_at(INDEX_TYPE n, INDEX_TYPE nz,
               const INDEX_TYPE *col_ptr, const INDEX_TYPE *row_idx,
               INDEX_TYPE *bnz, INDEX_TYPE **b_col_ptr, INDEX_TYPE **b_row_ind);

///A'+A
void a_plus_at_v2(INDEX_TYPE n, INDEX_TYPE nz,
                  const INDEX_TYPE *col_ptr, const INDEX_TYPE *row_idx,
                  INDEX_TYPE *bnz, INDEX_TYPE **b_col_ptr, INDEX_TYPE **b_row_ind);
#ifdef __cplusplus
}
#endif
#endif //QKLU_MATRIX_H
