//
// Created by mainf on 2024/4/15.
//

#ifndef QKLU_REORDER_H
#define QKLU_REORDER_H

#include "base/matrix.h"
#include "base/malloc.h"
#include "symbolic_analysis.h"
#include "preprocess.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * AMD重排序
 * @param matrix
 * @param preprocess_info 符号计算信息
 */
INDEX_TYPE *reorder_csr_amd(const SparseMatrix *matrix, PreprocessInfo *preprocess_info);

 /**
  * 求PAP^T
  * @param A
  * @param R PAP^T
  * @param iperm P^T
  */
 void apply_permutation(const SparseMatrix *A, SparseMatrix *R, const INDEX_TYPE iperm[]);

#ifdef __cplusplus
}
#endif

#endif //QKLU_REORDER_H
