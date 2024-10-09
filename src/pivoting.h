//
// Created by mainf on 2024/5/27.
//

#ifndef QKLU_PIVOTING_H
#define QKLU_PIVOTING_H

#ifdef __cplusplus
extern "C" {
#endif

#include "base/matrix.h"
#include "math.h"
#include "float.h"
#include "sys/types.h"

void static_pivoting(SparseMatrix *S, INDEX_TYPE **perm, INDEX_TYPE **iperm,
                     ELE_TYPE **row_scale, ELE_TYPE **col_scale);

#ifdef __cplusplus
}
#endif

#endif //QKLU_PIVOTING_H
