//
// Created by mainf on 2024/5/5.
//

#ifndef QKLU_VECTOR_UTIL_H
#define QKLU_VECTOR_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/types.h>
#include "matrix.h"

void random_int_vector(INDEX_TYPE vector[], INDEX_TYPE size);
void random_vector(ELE_TYPE vector[], INDEX_TYPE size);
double l2_norm(const ELE_TYPE *vector, INDEX_TYPE size);
double l1_norm(const ELE_TYPE *vector, INDEX_TYPE size);
double l_inf_norm(const ELE_TYPE *vector, INDEX_TYPE size);
void vector_sub(const ELE_TYPE *x, const ELE_TYPE *y, ELE_TYPE *r, INDEX_TYPE n);

#ifdef __cplusplus
}
#endif
#endif //QKLU_VECTOR_UTIL_H