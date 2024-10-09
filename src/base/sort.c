//
// Created by mainf on 2024/5/4.
//

#include "sort.h"

/**
 * 通过冒泡排序调整稀疏矩阵的index数组，使其变有序
 * @param Ai index
 * @param Ax value
 * @param n 矩阵纬度
 */
void bubbleSort(INDEX_TYPE Ai[], ELE_TYPE Ax[], INDEX_TYPE n) {
    INDEX_TYPE i, j, temp;
    ELE_TYPE value_temp;
    INDEX_TYPE swapped;
    for (i = 0; i < n - 1; i++) {
        swapped = 0; // 没有发生交换的标志
        // 内层循环进行相邻元素比较
        for (j = 0; j < n - i - 1; j++) {
            if (Ai[j] > Ai[j + 1]) {
                // 交换
                temp = Ai[j];
                Ai[j] = Ai[j + 1];
                Ai[j + 1] = temp;
                // 交换
                if (Ax != NULL) {
                    value_temp = Ax[j];
                    Ax[j] = Ax[j + 1];
                    Ax[j + 1] = value_temp;
                }
                swapped = 1; // 发生了交换
            }
        }
        // 如果在某次遍历中没有发生交换，则数组已经有序
        if (swapped == 0)
            break;
    }
}