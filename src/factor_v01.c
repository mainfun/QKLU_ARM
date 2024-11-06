#include <math.h>

#include "numerical.h"
#include "base/matrix.h"
#include "preprocess.h"

#define THRESHOLD 1e-8
/**
 * 循环展开
 */
void factor(const SparseMatrix *A, ELE_TYPE *Lx, ELE_TYPE *Ux,
            const INDEX_TYPE *Lp_start, const INDEX_TYPE *Lp_end, const INDEX_TYPE *Li,
            const INDEX_TYPE *Up_start, const INDEX_TYPE *Up_end, const INDEX_TYPE *Ui) {
    INDEX_TYPE n = A->num_row;
    INDEX_TYPE elimination_count = 0;
    //稀疏转稠密
    LOG_INFO("factor_v01 开始消元:");
    ELE_TYPE *D = csr2dense(A);
    clock_t factor_time = clock();
    //------------------------------------向下高斯消元------------------------------------
    ELE_TYPE scales[n];
    for (INDEX_TYPE i = 0; i < n; ++i) {//枚举列
        ELE_TYPE pivot = D[i * n + i]; // diag value
        if (pivot < THRESHOLD && pivot > -THRESHOLD) pivot = THRESHOLD;
        D[i * n + i]=pivot;
        for (INDEX_TYPE p = Lp_start[i]; p < Lp_end[i]; p++) {
            INDEX_TYPE j = Li[p];//行号
            //printf("第%lld行消第%lld行\n", i + 1, j + 1);
            //L的列
            D[j * n + i] /= pivot;
            scales[j] = D[j * n + i];
        }
        INDEX_TYPE k;
        ELE_TYPE *pivot_row_ptr = D + i * n;
        for (k = Up_start[i]+1; k < Up_end[i] - 6; k += 6) {
            //一次消多列
            INDEX_TYPE c0 = Ui[k];
            INDEX_TYPE c1 = Ui[k + 1];
            INDEX_TYPE c2 = Ui[k + 2];
            INDEX_TYPE c3 = Ui[k + 3];
            INDEX_TYPE c4 = Ui[k + 4];
            INDEX_TYPE c5 = Ui[k + 5];
            ELE_TYPE p0 = pivot_row_ptr[c0];
            ELE_TYPE p1 = pivot_row_ptr[c1];
            ELE_TYPE p2 = pivot_row_ptr[c2];
            ELE_TYPE p3 = pivot_row_ptr[c3];
            ELE_TYPE p4 = pivot_row_ptr[c4];
            ELE_TYPE p5 = pivot_row_ptr[c5];
            for (INDEX_TYPE p = Lp_start[i]; p < Lp_end[i]; p++) {
                INDEX_TYPE j = Li[p];//行号
                ELE_TYPE scale = scales[j];
                ELE_TYPE *eli_row_ptr = D + j * n;
                eli_row_ptr[c0] -= scale * p0;
                eli_row_ptr[c1] -= scale * p1;
                eli_row_ptr[c2] -= scale * p2;
                eli_row_ptr[c3] -= scale * p3;
                eli_row_ptr[c4] -= scale * p4;
                eli_row_ptr[c5] -= scale * p5;
            }
        }
        // //余量
        // for (; k < Up_end[i]; ++k) {
        //     for (INDEX_TYPE p = Lp_start[i] + 1; p < Lp_end[i]; p++) {
        //         INDEX_TYPE j = Li[p];//行号
        //         //L的列
        //         ELE_TYPE scale = D[j * n + i];
        //         ELE_TYPE *eli_row_ptr = D + j * n;
        //         INDEX_TYPE c0 = Ui[k];
        //         eli_row_ptr[c0] -= scale * pivot_row_ptr[c0];
        //     }
        // }
    }
    //---------------------------------向下高斯消元--end------------------------------------
    LOG_INFO("LU factor elapsed time: %lf ms", ((double) (clock() - factor_time)) / CLOCKS_PER_SEC * 1000.0);
    // if (A->nnz < 100) { LOG_INFO("D:"), print_dense_matrix(D, n); }
    LOG_DEBUG("消元次数为   ::::%lld\n", elimination_count);
    //写回
    INDEX_TYPE l_count_csc = 0;
    INDEX_TYPE u_count_csr = 0;
    for (INDEX_TYPE i = 0; i < n; i++) {
        //U csr
        for (INDEX_TYPE j = Up_start[i]; j < Up_end[i]; j++) {
            INDEX_TYPE index = Ui[j];
            Ux[u_count_csr++] = D[i * n + index];
            //printf("D[%lld][%lld]=%lf ", i, index,D[i * n + index]);
        }
        //L csc
        for (INDEX_TYPE j = Lp_start[i]; j < Lp_end[i]; j++) {
            INDEX_TYPE index = Li[j];
            Lx[l_count_csc++] = D[index * n + i];
        }
    }
    // if (A->nnz < 100) {
    //     LOG_INFO("U:");
    //     print_dense_matrix(csr2dense_v2(Up_start,Up_end,Ui,Ux,n), n);
    //     LOG_INFO("L:\n");
    //     print_dense_matrix(csc2dense_v2(Lp_start,Lp_end,Li,Lx,n), n);
    // }
    lu_free(D);
}
