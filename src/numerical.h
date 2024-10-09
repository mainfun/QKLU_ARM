//
// Created by mainf on 2024/9/17.
//

#ifndef QKLU_NUMERICAL_H
#define QKLU_NUMERICAL_H

#include "base/matrix.h"
#include "preprocess.h"


/**
* LU分解单核版本
* L csc格式，U csr格式
*/
void factor(const SparseMatrix *A,ELE_TYPE *Lx, ELE_TYPE *Ux,
            const INDEX_TYPE *Lp_start, const INDEX_TYPE *Lp_end, const INDEX_TYPE *Li,
            const INDEX_TYPE *Up_start, const INDEX_TYPE *Up_end, const INDEX_TYPE *Ui);

#endif //QKLU_NUMERICAL_H
