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
void bubble_sort(INDEX_TYPE Ai[], ELE_TYPE Ax[], INDEX_TYPE n) {
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




/**
 * 辅助函数：交换两个元素
 */
void swap(INDEX_TYPE *a, INDEX_TYPE *b) {
    INDEX_TYPE temp = *a;
    *a = *b;
    *b = temp;
}

/**
 * 快速排序的分区函数
 */
INDEX_TYPE partition(INDEX_TYPE Ai[], ELE_TYPE Ax[], INDEX_TYPE low, INDEX_TYPE high) {
    INDEX_TYPE pivot = Ai[high]; // 选择最后一个元素作为基准
    INDEX_TYPE i = low - 1; // 小于基准的元素的索引

    for (INDEX_TYPE j = low; j < high; j++) {
        if (Ai[j] < pivot) {
            i++;
            swap(&Ai[i], &Ai[j]); // 交换 Ai
            if (Ax != NULL) {
                swap((INDEX_TYPE*)&Ax[i], (INDEX_TYPE*)&Ax[j]); // 交换 Ax
            }
        }
    }
    swap(&Ai[i + 1], &Ai[high]); // 将基准元素放到正确位置
    if (Ax != NULL) {
        swap((INDEX_TYPE*)&Ax[i + 1], (INDEX_TYPE*)&Ax[high]); // 交换基准的 Ax
    }
    return i + 1; // 返回基准的索引
}

/**
 * 非递归快速排序
 */
void quick_sort_non_recursive(INDEX_TYPE Ai[], ELE_TYPE Ax[], INDEX_TYPE n) {
    INDEX_TYPE *stack = (INDEX_TYPE *)lu_malloc(n * sizeof(INDEX_TYPE)); // 创建栈
    INDEX_TYPE top = -1; // 栈顶指针

    // 初始化栈
    stack[++top] = 0; // 左边界
    stack[++top] = n - 1; // 右边界

    while (top >= 0) {
        // 弹出右边界和左边界
        INDEX_TYPE high = stack[top--];
        INDEX_TYPE low = stack[top--];

        // 进行分区
        INDEX_TYPE pi = partition(Ai, Ax, low, high);

        // 如果左边子数组有元素，则推入栈中
        if (pi - 1 > low) {
            stack[++top] = low;
            stack[++top] = pi - 1;
        }

        // 如果右边子数组有元素，则推入栈中
        if (pi + 1 < high) {
            stack[++top] = pi + 1;
            stack[++top] = high;
        }
    }

    lu_free(stack); // 释放栈内存
}

/**
 * 快速排序的入口函数
 */
void sort_sparse_matrix(INDEX_TYPE Ai[], ELE_TYPE Ax[], INDEX_TYPE n) {
    quick_sort_non_recursive(Ai, Ax, n);
}