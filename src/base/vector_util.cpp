//
// Created by mainf on 2024/5/5.
//

#include "vector_util.h"
#include <random>
#include <cmath>

void random_int_vector(INDEX_TYPE vector[], INDEX_TYPE size) {
    std::random_device rd;  // 用于获得随机种子
    std::mt19937 gen(rd()); // 以随机种子初始化 Mersenne Twister 引擎
    std::uniform_int_distribution<> dis(0, 99); // 生成0到99之间的随机数
    for (INDEX_TYPE i = 0; i < size; ++i) {
        vector[i] = dis(gen);
    }
}

void random_vector(ELE_TYPE vector[], INDEX_TYPE size) {
    std::random_device rd;  // 用于获得随机种子
    std::mt19937 gen(rd()); // 以随机种子初始化 Mersenne Twister 引擎
    std::uniform_real_distribution<> dis(10, 100); // 生成0到99之间的随机数
    for (INDEX_TYPE i = 0; i < size; ++i) {
        vector[i] = dis(gen);
    }
}

/// 计算向量的 L2 范数
double l2_norm(const ELE_TYPE *vector, INDEX_TYPE size) {
    ELE_TYPE sum = 0.0;
    for (INDEX_TYPE i = 0; i < size; i++) {
        sum += vector[i] * vector[i];
    }
    return sqrt(sum);
}

/// 计算向量的 L1 范数
double l1_norm(const ELE_TYPE *vector, INDEX_TYPE size) {
    double sum = 0.0;
    for (INDEX_TYPE i = 0; i < size; i++) {
        sum += fabs(vector[i]);
    }
    return sum;
}

// 计算向量的 L∞ 范数
double l_inf_norm(const ELE_TYPE *vector, INDEX_TYPE size) {
    ELE_TYPE max_val = 0.0;
    for (INDEX_TYPE i = 0; i < size; i++) {
        max_val = std::max(fabs(vector[i]), max_val);
    }
    return max_val;
}

/// r=x-y
void vector_sub(const ELE_TYPE *x, const ELE_TYPE *y, ELE_TYPE *r, INDEX_TYPE n) {
    for (INDEX_TYPE i = 0; i < n; ++i) {
        r[i] = x[i] - y[i];
    }
}