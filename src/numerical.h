//
// Created by mainf on 2024/9/17.
//

#ifndef QKLU_NUMERICAL_H
#define QKLU_NUMERICAL_H
#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include <layer_matrix.h>

#define BLOCK_SIDE 80
#define A(i,j) a[offset1[i]+(j)]
#define B(i,j) b[offset2[i]+(j)]
#define C(i,j) c[offset3[i]+(j)]
#define L(i,j) lx[offset_l[i]+(j)]
#define DA(i,j) a[BLOCK_SIDE*(i)+(j)]
#define DB(i,j) b[BLOCK_SIDE*(i)+(j)]
#define DC(i,j) c[BLOCK_SIDE*(i)+(j)]
#define DL(i,j) lx[(i)*BLOCK_SIDE+(j)]

/**
 * @details 分块LU分解，CPU单线程运行
 * @param l2 layer 2 matrix
*/
void block_factor(L2Matrix *l2);

void block_parallel_factor(L2Matrix *l2, int cut);

#ifdef __cplusplus
}
#endif //__cplusplus
#endif //QKLU_NUMERICAL_H
