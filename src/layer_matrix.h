//
// Created by mainf on 2024/10/9.
//

#ifndef LAYER_MATRIX_H
#define LAYER_MATRIX_H

#include <base/matrix.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    SPARSE, DENSE
} SparseMatrixFormat;

typedef struct {
    int *row_idx;
    int *col_idx;
    int *row_pointers;
    int *col_indices;
    ELE_TYPE *values;
    int *col_pointers;
    int *row_indices;
    int *offset;
    int d_size; //值数组的大小
    ELE_TYPE *csc_values;

    SparseMatrixFormat format;
    int u_nnz; // 非零元素的总数
    int l_nnz; // 非零元素的总数
    int side; //行（列）大小
    int num_row; //行个数（列大小）
    int num_col;
} BlockMatrix;


typedef struct {
    /**CSR**/
    int *row_pointers;
    int *col_indices;
    int *col_pointers;
    int *row_indices;
    BlockMatrix **block_matrices;
    int *diag_index;
    int *diag_index_csc;
    int num_row_block; //一列有多少行块
    int num_col_block;
    int block_count; //非0块个数
    INDEX_TYPE side; //整个矩阵行（列）大小
    INDEX_TYPE *block_side_sum; //block_side_sum[i]是前i个分块的边长的累加和
} L2Matrix;

typedef struct {
    int num_row_block;
    int num_col_block;
    int block_width;
    int block_height;
    int block_count; //非0块个数
    int *row_pointers;
    int *col_indices;
    L2Matrix **sub_matrices;
} L3Matrix;

BlockMatrix *get_diag_block(const L2Matrix *l2, int i);

BlockMatrix *get_block(const L2Matrix *l2, int i, int j);

void print_csr(BlockMatrix *mat, int n, int d_size);

void print_csc(BlockMatrix *mat, int n, int d_size);

BlockMatrix *init_BlockMatrix();

L2Matrix *init_L2Matrix();

/**
* @param Lp 符号分析结构
* @param Li 符号分析结构
* @param Up 符号分析结构
* @param Ui 符号分析结构
* @param l2 结果
* @param BLOCK_SIDE 块矩阵边长
* @param n 整个矩阵边长
*
*/
void csr2L2Matrix(const INDEX_TYPE *Lp, const INDEX_TYPE *Li,
                  const INDEX_TYPE *Up, const INDEX_TYPE *Ui,
                  const INDEX_TYPE *Ap, const INDEX_TYPE *Ai, const ELE_TYPE *Ax,
                  L2Matrix *l2, int BLOCK_SIDE, INDEX_TYPE n);

void calc_block_side(const INDEX_TYPE *Lp, const INDEX_TYPE *Li,
                     L2Matrix *l2, int BLOCK_SIDE, INDEX_TYPE n);

#ifdef __cplusplus
}
#endif
#endif //LAYER_MATRIX_H
