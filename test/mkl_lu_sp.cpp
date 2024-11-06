//
// Created by mainf on 2024/10/15.
//
#include <mkl.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <chrono>

#include "base/file.h"

int main(int argc, char *argv[]) {
    // if (argc < 2) {
    //     std::cerr << "Usage: " << argv[0] << " <matrix_file.mtx>" << std::endl;
    //     return 1;
    // }

    // const char *filename = argv[1];
    MKL_INT n, nnz;
    MKL_INT *ia = nullptr, *ja = nullptr;
    double *a = nullptr;
    printf("MKL_INT=%d\n",sizeof(MKL_INT));
    // 读取矩阵文件
    //read_mtx_file(filename, n, nnz, ia, ja, a);
    SparseMatrix *m = load_matrix_csr("/Users/mainf/其他/mtx/onetone1.mtx", false);
    n=m->num_row;
    nnz=m->nnz;
    ia=m->row_pointers;
    ja=m->col_indices;
    a=m->csr_values;
    printf("read_mtx_file end\n");
    // PARDISO 变量
    void *pt[64] = {0};
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    double ddum; // Double dummy
    MKL_INT idum; // Integer dummy
    MKL_INT nrhs = 1; // 右手边向量数量
    double *b = new double[n];
    double *x = new double[n];

    // 初始化PARDISO
    std::memset(iparm, 0, 64 * sizeof(MKL_INT));
    iparm[0] = 1;   // 使用默认设置
    iparm[1] = 2;   // 使用多线程
    iparm[7] = 2;   // 输出到屏幕
    iparm[9] = 13;

    maxfct = 1;     // 求解一个矩阵
    mnum = 1;       // 使用第一个矩阵
    msglvl = 0;     // 不输出
    error = 0;

    MKL_INT mtype = 11; // 对称正定矩阵 (如果非对称，则改为 13)

    // 填充b向量
    for (MKL_INT i = 0; i < n; ++i) {
        b[i] = 1.0;
        x[i] = 0.0;
    }

    // 计时开始
    auto start = std::chrono::high_resolution_clock::now();

    // PARDISO: 分析阶段
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        std::cerr << "Pardiso analysis error: " << error << std::endl;
        return 1;
    }else printf("分析阶段end\n");

    // PARDISO: 因式分解阶段
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        std::cerr << "Pardiso factorization error: " << error << std::endl;
        return 1;
    }

    // 计时结束
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "LU decomposition time: " << elapsed.count() << " seconds." << std::endl;

    // 释放PARDISO内存
    phase = -1; // 释放内存
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);

    // 释放分配的内存
    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;
    delete[] x;

    return 0;
}
