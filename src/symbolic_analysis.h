//
// Created by mainf on 2024/5/6.
//
#ifndef QKLU_SYMBOLIC_ANALYSIS_H
#define QKLU_SYMBOLIC_ANALYSIS_H
#ifdef __cplusplus
extern "C" {
#endif
#include <stdlib.h>
#include "base/malloc.h"
#include "base/matrix.h"
#include "base/forest.h"
#include "preprocess.h"

/**
 * 对称矩阵的符号分析
 * @param matrix 矩阵
 * @param info 预处理信息
 */
void symbolic_calc_sym(SparseMatrix *matrix, PreprocessInfo *info);

/**
 * 对称矩阵生成行消元树
 * @param csr_mat CSR格式矩阵
 * @param preprocess_info 预处理信息
 */
void create_row_etree_sym(SparseMatrix *csr_mat, PreprocessInfo *preprocess_info);

/**
 * 获取index数组中对角线元素的下标
 * @param A CSR和CSC的稀疏矩阵
 * @param preprocess_info 预处理信息
 */
void get_diag_index(SparseMatrix *A, PreprocessInfo *preprocess_info);

INDEX_TYPE *get_diag_index_v2(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai, INDEX_TYPE n);

/**
 * 返回L的列结构
 * @param csr_mat CSR格式矩阵
 * @param preprocess_info 预处理信息
 */
void get_l_col_patterns(const SparseMatrix *csr_mat, PreprocessInfo *preprocess_info);

/**
* L的行结构
*/
void get_l_row_patterns(const SparseMatrix *csr_mat, PreprocessInfo *info);

/**
* 将CSR格式矩阵转换为CSC格式
*/
void csr2csc_pattern(SparseMatrix *m);

void csc2csr_pattern(SparseMatrix *m);

void csc2csr(SparseMatrix *m);

void csr2csc(SparseMatrix *m);

/**
* 非对称符号分析
*/
void fill_in_2_no_sort_pruneL(INDEX_TYPE n, INDEX_TYPE nnz,
                              INDEX_TYPE *ai, INDEX_TYPE *ap,
                              INDEX_TYPE **L_rowpointer, INDEX_TYPE **L_columnindex,
                              INDEX_TYPE **U_rowpointer, INDEX_TYPE **U_columnindex);

#ifdef __cplusplus
}
#endif
#endif //QKLU_SYMBOLIC_ANALYSIS_H
