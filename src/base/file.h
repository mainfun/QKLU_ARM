//
// Created by mainf on 2024/4/16.
//

#ifndef QKLU_FILE_H
#define QKLU_FILE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <jpeglib.h>
#include <math.h>
#include <stdlib.h>
#include "matrix.h"
#include "file.h"

/**
头行的格式为：
%MatrixMarket object format field symmetry
头行必须是文件的第一行，并且必须以 %MatrixMarket 或 %%MatrixMarket 开始。其余4个字段的取值为：

object：通常是 matrix，另一个合法值是 vector，可以认为是矩阵的特殊形式；
format：coordinate 或者 array，前者用于稀疏矩阵，后者用于稠密矩阵或者向量；
field：合法值是 real, double, complex, integer or pattern，一般用 double 或者 complex，单精度使用 real，整数矩阵使用 integer；如果是 pattern，必须使用coordinate格式，并且只列出非零元所在的行列；
symmetry：是否对称，合法值包括 general（一般矩阵）、symmetric（对阵矩阵）、skew-symmetric（转置取负）和hermitian（复数对称）。出了general，其余三种格式要求只列出下三角数据。
*/
typedef struct {
    INDEX_TYPE rows;
    INDEX_TYPE cols;
    INDEX_TYPE not_zeros;
    char object[20];
    char format[20];
    char field[20];
    char symmetry[20];
} MtxInfo;

FILE *file_open(const char *file_name);

SparseMatrix *load_matrix_csr(const char *file_name, bool base_0);

SparseMatrix *load_matrix_csc_order(const char *file_name, bool base_0);

void csr2mtx(const char *file_name, const SparseMatrix *A);

void csc2mtx(const char *file_name, const SparseMatrix *A);

void L2mtx(const char *file_name, const SparseMatrix *A);

void int_vector2mtx(const char *file_name, const INDEX_TYPE *vector, INDEX_TYPE n);

void fp_vector2mtx(const char *file_name, const ELE_TYPE *vector, INDEX_TYPE n);

#ifdef __cplusplus
}
#endif

#endif //QKLU_FILE_H
