//
// Created by mainf on 2024/5/4.
//

#ifndef QKLU_SORT_H
#define QKLU_SORT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "base/matrix.h"

void bubble_sort(INDEX_TYPE Ai[], ELE_TYPE Ax[], INDEX_TYPE n);

void sort_sparse_matrix(INDEX_TYPE Ai[], ELE_TYPE Ax[], INDEX_TYPE n);

#ifdef __cplusplus
}
#endif

#endif //QKLU_SORT_H
