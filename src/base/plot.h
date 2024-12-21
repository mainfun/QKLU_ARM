//
// Created by mainf on 2024/9/15.
//

#ifndef QKLU_PLOT_H
#define QKLU_PLOT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <jpeglib.h>
#include <math.h>
#include <stdlib.h>
#include "matrix.h"

/**
 * @param matrix
 * @param filename
 * @param img_width image max width
 * @param img_height image max height
 */
void csr2image(const SparseMatrix *matrix, const char *filename, int img_width, int img_height);

void csr2_RGB_image(const SparseMatrix *matrix, const char *filename, int img_width, int img_height);

/**
 * 画出矩阵的局部分块的结构
 */
void csr2image_block(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai,
                     const char *filename,
                     INDEX_TYPE row_start, INDEX_TYPE row_end,
                     INDEX_TYPE col_start, INDEX_TYPE col_end);

void csr2image_block_cut(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai,
                         const char *filename,
                         INDEX_TYPE row_start, INDEX_TYPE row_end,
                         INDEX_TYPE col_start, INDEX_TYPE col_end,
                         const INDEX_TYPE *cut_points, INDEX_TYPE cut_points_count);

///带坐标系
void csr2_RGB_image_coor(const SparseMatrix *matrix, const char *filename,
                         int img_width, int img_height);

#ifdef __cplusplus
}
#endif

#endif //QKLU_PLOT_H
