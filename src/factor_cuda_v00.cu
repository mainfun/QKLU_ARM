//
// Created by mainf on 2024/9/19.
//
extern "C" {
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define N 1024        // 矩阵大小
#define NB 32         // 块大小
#define NUM_BLOCKS (N / NB)

// 内核函数：执行块LU分解（GETRF）
__global__ void GETRF(float *A, int n, int nb) {
    __shared__ float block[NB][NB];
    const unsigned int tx = threadIdx.x;
    const unsigned int ty = threadIdx.y;

    // 加载当前块到共享内存
    block[ty][tx] = A[blockIdx.y * nb * n + blockIdx.x * nb + ty * n + tx];
    __syncthreads();

    // 进行LU分解（不选主元）
    for (int k = 0; k < nb; ++k) {
        if (tx >= k && ty == k) {
            block[tx][k] /= block[k][k]; // 计算L部分
        }
        __syncthreads();
        if (tx > k && ty > k) {
            block[ty][tx] -= block[ty][k] * block[k][tx]; // 更新U部分
        }
        __syncthreads();
    }

    // 将结果写回全局内存
    A[blockIdx.y * nb * n + blockIdx.x * nb + ty * n + tx] = block[ty][tx];
}

// 内核函数：更新行块（GESSM）
__global__ void GESSM(float *A, float *L, int n, int nb, int k) {
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int col_offset = (k + 1) * nb;

    for (int i = k + 1; i < n / nb; ++i) {
        // 更新行块
        A[(i * nb + ty) * n + col_offset + tx] -= L[ty * n + col_offset + tx] * A[k * nb * n + col_offset + tx];
    }
}

// 内核函数：更新列块（TSTRF）
__global__ void TSTRF(float *A, float *U, int n, int nb, int k) {
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int row_offset = (k + 1) * nb;

    for (int i = k + 1; i < n / nb; ++i) {
        // 更新列块
        A[row_offset * n + i * nb + tx] -= U[ty * n + row_offset + tx] * A[row_offset * n + k * nb + tx];
    }
}

// 内核函数：更新剩余矩阵（SSSM）
__global__ void SSSM(float *A, float *L, float *U, int n, int nb, int k) {
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int row_offset = (k + 1) * nb;
    int col_offset = (k + 1) * nb;

    for (int i = k + 1; i < n / nb; ++i) {
        for (int j = k + 1; j < n / nb; ++j) {
            // Schur补矩阵更新
            A[(i * nb + ty) * n + (j * nb + tx)] -= L[ty * n + k * nb + tx] * U[k * nb * n + j * nb + tx];
        }
    }
}
}

// 主函数
int main() {
    // 初始化矩阵A
    float *h_A = (float *) malloc(N * N * sizeof(float));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                h_A[i * N + j] = 2.0f;
            } else {
                h_A[i * N + j] = 1.0f;
            }
        }
    }

    // 分配GPU内存
    float *d_A;
    cudaMalloc((void **) &d_A, N * N * sizeof(float));
    cudaMemcpy(d_A, h_A, N * N * sizeof(float), cudaMemcpyHostToDevice);

    // 配置CUDA网格和块
    dim3 threadsPerBlock(NB, NB);
    dim3 numBlocks(NUM_BLOCKS, NUM_BLOCKS);

    // 分块LU分解流程
    for (int k = 0; k < N / NB; ++k) {
        // 1. GETRF: 分解主块
        GETRF<<<dim3(1, 1), threadsPerBlock>>>(d_A, N, NB);
        cudaDeviceSynchronize();

        // 2. GESSM: 更新行块
        GESSM<<<dim3(NUM_BLOCKS - (k + 1), 1), threadsPerBlock>>>(d_A, d_A, N, NB, k);
        cudaDeviceSynchronize();

        // 3. TSTRF: 更新列块
        TSTRF<<<dim3(1, NUM_BLOCKS - (k + 1)), threadsPerBlock>>>(d_A, d_A, N, NB, k);
        cudaDeviceSynchronize();

        // 4. SSSM: 更新剩余块
        SSSM<<<dim3(NUM_BLOCKS - (k + 1), NUM_BLOCKS - (k + 1)), threadsPerBlock>>>(d_A, d_A, d_A, N, NB, k);
        cudaDeviceSynchronize();
    }

    // 将结果复制回主机
    cudaMemcpy(h_A, d_A, N * N * sizeof(float), cudaMemcpyDeviceToHost);

    // 输出部分结果
    for (int i = 0; i < min(N, 10); i++) {
        for (int j = 0; j < min(N, 10); j++) {
            printf("%f ", h_A[i * N + j]);
        }
        printf("\n");
    }

    // 释放内存
    free(h_A);
    cudaFree(d_A);

    return 0;
}
