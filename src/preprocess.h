//
// Created by mainf on 2024/8/16.
//
#ifndef QKLU_PREPROCESS_H
#define QKLU_PREPROCESS_H

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

#include "base/matrix.h"
#include "base/forest.h"

typedef struct PREPROCESS_INFO {
    struct SPARSE_MATRIX *pattern; //L+U的结构
    double max_nz; //预估的L的非零元上限
    Forest *row_etree; //行消元树
    struct SPARSE_MATRIX *L; //L的结构
    struct SPARSE_MATRIX *U; //U的结构
    bool pattern_sym; //结构对称
    bool numeric_sym; //数组对称
    double pattern_sym_rate; // 结构对称元素个数/nz*100%
    double numeric_sym_rate; // 数值对称元素个数/nz*100%
    INDEX_TYPE pattern_asym_number; //结构不对称个数
    INDEX_TYPE numeric_asym_number; //数值不对称个数
    INDEX_TYPE *diag_index_csr; //对角线索引
    INDEX_TYPE *diag_index_csc; //对角线索引
    INDEX_TYPE *reorder_perm; //重排序的置换矩阵
    INDEX_TYPE *reorder_iperm; //重排序的置换矩阵的转置矩阵
    INDEX_TYPE *mc64_perm; //选主元的置换矩阵
    INDEX_TYPE *mc64_iperm;
    double *Dc; //选主元的缩放矩阵,col_scale
    double *Dr; //row_scale
} PreprocessInfo;

PreprocessInfo *init_preprocess_info();

void free_preprocess_info(PreprocessInfo *info);

/**
 * CPU单线程计算的预处理
 * @param matrix 原矩阵
 * @param info 预处理信息
 * @param is_static_pivoting 是否静态选主元
 * @param is_reorder 是否重排序
 * @param is_symbolic_calc 是否符号分析
*/
SparseMatrix *preprocess(SparseMatrix *matrix, PreprocessInfo *info, bool is_static_pivoting, bool is_reorder,
                bool is_symbolic_calc);

/**
 * 获取矩阵对称信息
 * @param matrix 原矩阵
 * @param info 预处理信息
 */
void calculate_asymmetry(const SparseMatrix *matrix, PreprocessInfo *info);

#ifdef __cplusplus
}
#endif //__cplusplus
#endif //QKLU_PREPROCESS_H
