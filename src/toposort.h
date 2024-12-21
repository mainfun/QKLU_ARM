//
// Created by mainf on 2024/11/28.
//

#ifndef TOPOSORT_H
#define TOPOSORT_H

#ifdef __cplusplus
extern "C" {
#endif
#include<base/malloc.h>
#include <base/matrix.h>

INDEX_TYPE *reorder_toposort(SparseMatrix *A, INDEX_TYPE n, INDEX_TYPE *cut_point);

#ifdef __cplusplus
}
#endif

#endif //TOPOSORT_H
