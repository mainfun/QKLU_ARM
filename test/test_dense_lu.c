//
// Created by mainf on 2024/10/11.
//

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <base/matrix.h>

#define ELE_TYPE double
#define BLOCK_SIDE 9
#define A(i,j) a[(i)*BLOCK_SIDE+(j)]
#define B(i,j) b[(i)*BLOCK_SIDE+(j)]
#define C(i,j) c[(i)*BLOCK_SIDE+(j)]

void LUD(ELE_TYPE *a) {
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        ELE_TYPE pivot = A(i, i);
        for (int j = i + 1; j < BLOCK_SIDE; ++j) {
            ELE_TYPE l = A(j, i) /= pivot;
            for (int k = i + 1; k < BLOCK_SIDE; ++k) {
                A(j, k) -= l * A(i, k);
            }
        }
    }
}

void generateDiagonalDominantMatrix(ELE_TYPE *a) {
    // 初始化随机数生成器
    srand(time(NULL));

    for (int i = 0; i < BLOCK_SIDE; i++) {
        int sum = 0; // 记录非对角线元素的和
        for (int j = 0; j < BLOCK_SIDE; j++) {
            if (i == j) {
                // 生成随机的对角线元素，确保其大于其他元素的和
                int diagonal_value = rand() % 10 + 1; // 生成 1 到 10 的随机数
                A(i, j) = diagonal_value + sum + (rand() % 5 + 1); // 确保对角线元素占优
            } else {
                // 生成随机的非对角线元素
                A(i, j) = rand() % 5; // 生成 0 到 4 的随机数
                sum += abs(A(i, j)); // 更新非对角线元素的和
            }
        }
    }
}

void printMatrix(ELE_TYPE *a) {
    for (int i = 0; i < BLOCK_SIDE; i++) {
        for (int j = 0; j < BLOCK_SIDE; j++) {
            printf("%14.10lf ", A(i, j));
        }
        printf("\n");
    }
}

void verify_lu(const ELE_TYPE *a, const ELE_TYPE *original_matrix) {
    ELE_TYPE *b = (ELE_TYPE *) malloc(BLOCK_SIDE * BLOCK_SIDE * sizeof(ELE_TYPE));
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int k = i; k < BLOCK_SIDE; ++k) {
            B(i, k) = A(i, k);
        }
    }
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        //算B[i]
        for (int j = 0; j < i; ++j) {
            for (int k = j; k < BLOCK_SIDE; ++k) {
                //B[i][k]+=L[i][j]*U[j][k]
                B(i, k) += A(i, j) * A(j, k);
            }
        }
    }
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int k = 0; k < BLOCK_SIDE; ++k) {
            B(i, k) -= original_matrix[i * BLOCK_SIDE + k];
        }
    }
    printf("B:\n");
    printMatrix(b);
    free(b);
}

int main() {
    // 定义矩阵
    ELE_TYPE *matrix = (ELE_TYPE *) malloc(BLOCK_SIDE * BLOCK_SIDE * sizeof(ELE_TYPE));
    ELE_TYPE *original_matrix = (ELE_TYPE *) malloc(BLOCK_SIDE * BLOCK_SIDE * sizeof(ELE_TYPE));
    // 生成随机对角占优矩阵
    generateDiagonalDominantMatrix(matrix);
    memcpy(original_matrix, matrix, BLOCK_SIDE * BLOCK_SIDE * sizeof(ELE_TYPE));

    // 打印矩阵
    printf("Generated Diagonal Dominant Matrix:\n");
    printMatrix(matrix);
    LUD(matrix);
    printf("L+U-E:\n");
    printMatrix(matrix);

    verify_lu(matrix, original_matrix);
    free(matrix);
    return 0;
}
