//
// Created by mainf on 2024/8/17.
//

#ifndef MY_LU_SOLVING_H
#define MY_LU_SOLVING_H
#ifdef __cplusplus
extern "C" {
#endif

#include <layer_matrix.h>

#include "base/matrix.h"
#include "preprocess.h"
#include "base/file.h"
#include "numerical.h"

void solve(SparseMatrix *A, ELE_TYPE *x, const ELE_TYPE *b);

void solve_no_pivoting(SparseMatrix *A, ELE_TYPE *x, const ELE_TYPE *b);

void lower_solver(const SparseMatrix *L, ELE_TYPE *x, const ELE_TYPE *b);

void upper_solver(const SparseMatrix *U, ELE_TYPE *x, const ELE_TYPE *b);

void upper_solver_csc(const SparseMatrix *U, ELE_TYPE *x, const ELE_TYPE *b);

void forward_substitution(const SparseMatrix *L, const INDEX_TYPE *reorder_iperm,
                          const INDEX_TYPE *reorder_perm, const ELE_TYPE *Dr,
                          const INDEX_TYPE *mc64_perm, const INDEX_TYPE *mc64_iperm,
                          const ELE_TYPE *b, ELE_TYPE *y, INDEX_TYPE n);

void backward_substitution(const SparseMatrix *U, const INDEX_TYPE *reorder_perm,
                           const INDEX_TYPE *reorder_iperm, const ELE_TYPE *Dc,
                           ELE_TYPE *x, const ELE_TYPE *y, INDEX_TYPE n);

void backward_substitution_block(const L2Matrix *U, const INDEX_TYPE *reorder_perm,
                                 const INDEX_TYPE *reorder_iperm, const ELE_TYPE *Dc,
                                 ELE_TYPE *x, const ELE_TYPE *y, INDEX_TYPE n,
                                 const INDEX_TYPE *toposort_iperm);

void forward_substitution_block(const L2Matrix *L, const INDEX_TYPE *reorder_iperm,
                                const INDEX_TYPE *reorder_perm, const ELE_TYPE *Dr,
                                const INDEX_TYPE *mc64_perm, const INDEX_TYPE *mc64_iperm,
                                const ELE_TYPE *b, ELE_TYPE *y, INDEX_TYPE n,
                                const INDEX_TYPE *toposort_iperm);
#ifdef __cplusplus
}
#endif
#endif //MY_LU_SOLVING_H
