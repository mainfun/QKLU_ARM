#include <math.h>

#include "numerical.h"
#include "base/matrix.h"
#include "preprocess.h"
#include "base/base_math.h"

#define THRESHOLD 1e-8
/**
 * LU分解的朴素实现(48ms~65ms)
 */
void GETRF(const SparseMatrix *A, ELE_TYPE *Lx, ELE_TYPE *Ux,
           const INDEX_TYPE *Lp_start, const INDEX_TYPE *Lp_end, const INDEX_TYPE *Li,
           const INDEX_TYPE *Up_start, const INDEX_TYPE *Up_end, const INDEX_TYPE *Ui) {
    INDEX_TYPE n = A->num_row;
    INDEX_TYPE elimination_count = 0;
    //稀疏转稠密
    ELE_TYPE *D = csr2dense(A);
    LOG_INFO("开始消元:");
    clock_t factor_time = clock();
    //向下高斯消元
    for (INDEX_TYPE i = 0; i < n; ++i) {
        //枚举列
        ELE_TYPE pivot = D[i * n + i]; // diag value
        if (fabs(pivot) < THRESHOLD) pivot = THRESHOLD;
        D[i * n + i] = pivot;
        //#pragma omp parallel for num_threads(6)
        for (INDEX_TYPE p = Lp_start[i]; p < Lp_end[i]; p++) {
            INDEX_TYPE j = Li[p]; //行号
            //printf("第%lld行消第%lld行\n", i + 1, j + 1);
            ELE_TYPE scale = D[j * n + i] / pivot;
            //L的列
            D[j * n + i] = scale;
            ELE_TYPE *pivot_row_ptr = D + i * n;
            ELE_TYPE *eli_row_ptr = D + j * n;
            //优化访存
            for (INDEX_TYPE k = Up_start[i] + 1; k < Up_end[i]; k++) {
                INDEX_TYPE c = Ui[k];
                eli_row_ptr[c] -= scale * pivot_row_ptr[c];
                elimination_count++;
            }
        }
        //printf("elimination_count:%lld\n", elimination_count);
    }
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

/**
 * LU分解的分块实现
 */
void factor(const SparseMatrix *A, ELE_TYPE *Lx, ELE_TYPE *Ux,
            const INDEX_TYPE *Lp_start, const INDEX_TYPE *Lp_end, const INDEX_TYPE *Li,
            const INDEX_TYPE *Up_start, const INDEX_TYPE *Up_end, const INDEX_TYPE *Ui) {
    INDEX_TYPE n = A->num_row;
    INDEX_TYPE elimination_count = 0;
    LOG_INFO("开始消元:（分块法）");
    clock_t factor_time = clock();
    const int block_height = 100;
    const int block_width = 100;
    INDEX_TYPE num_row_block = CEIL_DIV(A->num_row, block_height); //一行多少块
    INDEX_TYPE num_row_remainder = A->num_row % block_height; //不足一行的余量

    INDEX_TYPE num_col_block = CEIL_DIV(A->num_col, block_width);
    INDEX_TYPE num_col_remainder = A->num_col % block_width;

    INDEX_TYPE *Lp_block_start = (INDEX_TYPE *) lu_malloc(sizeof(INDEX_TYPE) * block_width);
    INDEX_TYPE *Lp_block_end = (INDEX_TYPE *) lu_malloc(sizeof(INDEX_TYPE) * block_width);
    INDEX_TYPE *Up_block_start = (INDEX_TYPE *) lu_malloc(sizeof(INDEX_TYPE) * block_height);
    INDEX_TYPE *Up_block_end = (INDEX_TYPE *) lu_malloc(sizeof(INDEX_TYPE) * block_height);

    memcpy(Lp_block_start, Lp_start, block_width);
    memcpy(Up_block_start, Up_start, block_height);
    for (INDEX_TYPE k = 0; k < A->num_row; k += block_height) {
        //U csr 分块
        for (INDEX_TYPE i = 0; i < block_height; i++) {
            INDEX_TYPE j;
            for (j = Up_start[i]; j < Up_end[i]; j++) {
                INDEX_TYPE col_idx = Ui[j];
                if (col_idx >= k+block_width) {
                    Up_block_end[i] = j;
                    break;
                }
            }
            Up_block_end[i] = j;
        }
        //L csc 分块
        for (INDEX_TYPE i = 0; i < block_width; i++) {
            INDEX_TYPE j;
            for (j = Lp_start[i]; j < Lp_end[i]; j++) {
                INDEX_TYPE row_idx = Ui[j];
                if (row_idx >= k+block_height) {
                    Lp_block_end[i] = j;
                    break;
                }
            }
            Lp_block_end[i] = j;
        }
        //GETRF();
        memcpy(Lp_block_start, Lp_block_end, block_width);
        memcpy(Up_block_start, Up_block_end, block_height);
        //GESSM

        //TSTRF
        //SSSSM
    }
    lu_free(Up_block_start);
    lu_free(Up_block_end);
    lu_free(Lp_block_start);
    lu_free(Lp_block_end);
    LOG_INFO("LU factor elapsed time: %lf ms", ((double) (clock() - factor_time)) / CLOCKS_PER_SEC * 1000.0);
}
