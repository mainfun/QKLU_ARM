//
// Created by mainf on 2024/8/18.
//
#include "check.h"
#include "math.h"

#include "preprocess.h"

ELE_TYPE check_solving(const SparseMatrix *A, ELE_TYPE x[], ELE_TYPE b[]) {
    INDEX_TYPE n = A->num_row;
    ELE_TYPE *r = (ELE_TYPE *) lu_malloc(n * sizeof(ELE_TYPE));
    ELE_TYPE *y = (ELE_TYPE *) lu_malloc(n * sizeof(ELE_TYPE));
    SpMV_csr(A, x, y); //y=Ax
    vector_sub(y, b, r, n); //r=Ax-b
    ELE_TYPE l2 = l2_norm(r, n); //|| Ax - b ||
    ELE_TYPE result = l2 / l2_norm(b, n);
    LOG_INFO("|| Ax - b || / || b || = %le", result);
    lu_free(y);
    lu_free(r);
    return result;
}

double check_lower_solving(const SparseMatrix *L, const ELE_TYPE x[], ELE_TYPE b[]) {
    INDEX_TYPE n = L->num_row;
    ELE_TYPE *y = (ELE_TYPE *) lu_malloc(n * sizeof(ELE_TYPE));
    ELE_TYPE *r = (ELE_TYPE *) lu_malloc(n * sizeof(ELE_TYPE));
    /* y=(L+E)x */
    for (INDEX_TYPE i = 0; i < n; ++i) {
        y[i] = 0.0;
    }
    for (INDEX_TYPE col = 0; col < n; ++col) {
        // 首先加上单位矩阵的贡献，即x[col]
        y[col] += x[col];
        // 遍历该列中的每一个非零元素
        for (INDEX_TYPE idx = L->col_pointers[col]; idx < L->col_pointers[col + 1]; ++idx) {
            INDEX_TYPE row = L->row_indices[idx];
            ELE_TYPE value = L->csc_values[idx];
            y[row] += value * x[col];
        }
    }
    vector_sub(y, b, r, n);
    ELE_TYPE l2 = l2_norm(r, n);
    LOG_INFO("check_lower_solver: l2_norm: %lg", l2);
    return l2;
}

double check_upper_solving(const SparseMatrix *U, const ELE_TYPE x[], ELE_TYPE b[]) {
    INDEX_TYPE n = U->num_row;
    ELE_TYPE *y = (ELE_TYPE *) lu_malloc(n * sizeof(ELE_TYPE));
    ELE_TYPE *r = (ELE_TYPE *) lu_malloc(n * sizeof(ELE_TYPE));
    /* y=Ux */
    SpMV_csr(U, x, y);
    vector_sub(y, b, r, n);
    ELE_TYPE l2 = l2_norm(r, n);
    LOG_INFO("check_upper_solver: l2_norm: %lg", l2);
    return l2;
}

/**
 * @param RMSE 均方根误差
 * @param max_diff 绝对值最大误差
 */
void calc_RMSE_and_max_diff(const SparseMatrix *A, const SparseMatrix *L,
                            const SparseMatrix *U, ELE_TYPE *RMSE, ELE_TYPE *max_diff) {
    LOG_INFO("check lu start.");
    const INDEX_TYPE *Ap = A->row_pointers;
    const INDEX_TYPE *Ai = A->col_indices;
    const ELE_TYPE *Ax = A->csr_values;

    const INDEX_TYPE *Lp = L->row_pointers;
    const INDEX_TYPE *Li = L->col_indices;
    const ELE_TYPE *Lx = L->csr_values;

    const INDEX_TYPE *Up = U->row_pointers;
    const INDEX_TYPE *Ui = U->col_indices;
    const ELE_TYPE *Ux = U->csr_values;

    ELE_TYPE *R = (ELE_TYPE *) lu_malloc(A->num_col * sizeof(ELE_TYPE));
    ELE_TYPE sum_squared_diff = 0.0;
    ELE_TYPE max_absolute_diff = 0.0;
    for (INDEX_TYPE i = 0; i < L->num_row; ++i) {
        //初始化
        memset(R, 0, sizeof(ELE_TYPE) * A->num_col);
        for (INDEX_TYPE j = Ap[i]; j < Ap[i + 1]; j++) {
            R[Ai[j]] = -Ax[j];
        }
        for (INDEX_TYPE j = Up[i]; j < Up[i + 1]; j++) {
            R[Ui[j]] = Ux[j];
        }
        //计算
        for (INDEX_TYPE j = Lp[i]; j < Lp[i + 1]; j++) {
            ELE_TYPE Lv = Lx[j];
            INDEX_TYPE l_col = Li[j];
            for (INDEX_TYPE k = Up[l_col]; k < Up[l_col + 1]; k++) {
                ELE_TYPE Uv = Ux[k];
                INDEX_TYPE col = Ui[k];
                R[col] += Lv * Uv;
            }
        }
        //计算误差
        for (INDEX_TYPE j = Ap[i]; j < Ap[i + 1]; j++) {
            ELE_TYPE v = R[Ai[j]];
            //LOG_DEBUG("%lf ",v);
            sum_squared_diff += v * v;
            double abs_diff = fabs(v);
            if (abs_diff > max_absolute_diff) {
                max_absolute_diff = abs_diff;
            }
        }
    }
    lu_free(R);
    *RMSE = sqrt(sum_squared_diff / (ELE_TYPE) A->nnz);
    *max_diff = max_absolute_diff;
}

void check_lu(const SparseMatrix *A, const SparseMatrix *L, const SparseMatrix *U) {
    ELE_TYPE RMSE = 0;
    ELE_TYPE max = 0;
    calc_RMSE_and_max_diff(A, L, U, &RMSE, &max);
    LOG_DEBUG("RMSE=%le, max=%le", RMSE, max);
}
