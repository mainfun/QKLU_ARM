//
// Created by mainf on 2024/8/18.
//

#ifndef MY_LU_CHECK_H
#define MY_LU_CHECK_H
#ifdef __cplusplus
extern "C" {
#endif

#include "base/matrix.h"
#include "base/vector_util.h"

/**
 * 解向量误差检验
 * @param A 系数矩阵
 * @param x 解向量
 * @param b 值向量
 * @return || Ax - b || / || b ||
 */
ELE_TYPE check_solving(const SparseMatrix *A, ELE_TYPE x[], ELE_TYPE b[]);

ELE_TYPE check_lower_solving(const SparseMatrix *L, const ELE_TYPE x[], ELE_TYPE b[]);

ELE_TYPE check_upper_solving(const SparseMatrix *U, const ELE_TYPE x[], ELE_TYPE b[]);

/**
 * LU分解误差检验
 * @param A SparseMatrix
 * @param L lower matrix
 * @param U upper matrix
 */
void check_lu(const SparseMatrix *A, const SparseMatrix *L, const SparseMatrix *U);

#ifdef __cplusplus
}
#endif
#endif //MY_LU_CHECK_H
