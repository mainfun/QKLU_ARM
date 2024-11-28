//
// Created by mainf on 2024/8/16.
//

#include "preprocess.h"

#include <toposort.h>

#include "pivoting.h"
#include "reorder.h"
#include "base/file.h"

PreprocessInfo *init_preprocess_info() {
    PreprocessInfo *info = (PreprocessInfo *) lu_malloc(sizeof(PreprocessInfo));
    info->pattern = NULL;
    info->diag_index_csc = NULL;
    info->diag_index_csr = NULL;
    info->max_nz = 0;
    info->row_etree = NULL;
    info->L = NULL;
    info->U = NULL;
    info->numeric_asym_number = 0;
    info->numeric_sym_rate = 0;
    info->numeric_sym = false;
    info->pattern_asym_number = 0;
    info->pattern_sym_rate = 0;
    info->pattern_sym = false;
    info->reorder_iperm = NULL;
    info->reorder_perm = NULL;
    info->Dc = NULL;
    info->Dr = NULL;
    info->mc64_iperm = NULL;
    info->mc64_perm = NULL;
    info->Dc = NULL;
    info->Dr = NULL;
    return info;
}

void free_preprocess_info(PreprocessInfo *info) {
    lu_free(info->diag_index_csc);
    lu_free(info->diag_index_csr);

    lu_free(info->row_etree);
    free_sparse_matrix(info->L);
    free_sparse_matrix(info->U);

    lu_free(info->mc64_iperm);
    lu_free(info->mc64_perm);
    lu_free(info->reorder_iperm);
    lu_free(info->reorder_perm);
    lu_free(info->Dc);
    lu_free(info->Dr);
    lu_free(info);
}

void apply_pivoting(SparseMatrix *A, SparseMatrix *R,
                    const INDEX_TYPE perm[],
                    const ELE_TYPE row_scale[], const ELE_TYPE col_scale[]) {
    INDEX_TYPE *Ap = A->row_pointers;
    INDEX_TYPE *Ai = A->col_indices;
    ELE_TYPE *Ax = A->csr_values;
    INDEX_TYPE n = A->num_row;
    //求PAP^T
    //计算Rp
    for (INDEX_TYPE i = 0; i < n; i++) {
        INDEX_TYPE row_num = Ap[i + 1] - Ap[i];
        R->row_pointers[perm[i] + 1] = row_num;
    }
    R->row_pointers[0] = 0;
    for (INDEX_TYPE i = 0; i < n; i++) {
        R->row_pointers[i + 1] += R->row_pointers[i];
    }
    for (INDEX_TYPE i = 0; i < n; i++) {
        INDEX_TYPE row = perm[i];
        ELE_TYPE rs = row_scale[i];
        INDEX_TYPE tmp_index = R->row_pointers[row];
        for (INDEX_TYPE j = Ap[i]; j < Ap[i + 1]; j++) {
            INDEX_TYPE col = Ai[j];
            R->col_indices[tmp_index] = col;
            R->csr_values[tmp_index++] = Ax[j] * rs * col_scale[col];
        }
    }
}

//获取矩阵对称信息
void calculate_asymmetry(const SparseMatrix *matrix, PreprocessInfo *info) {
    INDEX_TYPE pattern_asym_count = 0;
    INDEX_TYPE numeric_asym_count = 0;

    for (INDEX_TYPE i = 0; i < matrix->num_row; ++i) {
        INDEX_TYPE j = info->diag_index_csr[i];
        while (j < matrix->row_pointers[i + 1]) {
            INDEX_TYPE csr_col = matrix->col_indices[j];
            ELE_TYPE csr_value = matrix->csr_values[j];
            INDEX_TYPE k = info->diag_index_csc[i];
            while (k < matrix->col_indices[i + 1]) {
                INDEX_TYPE csc_row = matrix->row_indices[k];
                ELE_TYPE csc_value = matrix->csc_values[k];
                if (csr_col == csc_row) {
                    pattern_asym_count++;
                    if (csr_value == csc_value) {
                        numeric_asym_count++;
                    }
                } else if (csr_col > csc_row) {
                    k++;
                } else {
                    j++;
                }
            }
        }
    }

    info->pattern_asym_number = pattern_asym_count;
    info->numeric_asym_number = numeric_asym_count;
    info->numeric_sym_rate = (double) numeric_asym_count / (double) (matrix->nnz - matrix->num_row);
    info->pattern_sym_rate = (double) pattern_asym_count / (double) (matrix->nnz - matrix->num_row);
    info->numeric_sym = !info->numeric_asym_number;
    info->pattern_sym = !info->pattern_asym_number;
    LOG_DEBUG("\nnumerical_asym = %lld, pattern_asym = %lld\n", numeric_asym_count, pattern_asym_count);
    LOG_DEBUG("\nnumeric asymmetry占比%.2lf%%\npattern asymmetry占比%.2lf%%\n",
              info->numeric_sym_rate, info->pattern_sym_rate);
}

SparseMatrix *preprocess(SparseMatrix *matrix, PreprocessInfo *info,
                         bool is_static_pivoting, bool is_reorder, bool is_symbolic_calc) {
    clock_t preprocess_time = clock();
    INDEX_TYPE n = matrix->num_row;
    INDEX_TYPE nnz = matrix->nnz;
    // //CSR->CSC
    // if (matrix->col_pointers == NULL) matrix->col_pointers = matrix->row_pointers;
    // if (matrix->row_indices == NULL) matrix->row_indices = matrix->col_indices;
    // //CSC->CSR
    // if (matrix->row_pointers == NULL) matrix->row_pointers = matrix->col_pointers;
    // if (matrix->col_indices == NULL) matrix->col_indices = matrix->row_indices;
    // get_diag_index(matrix, info);
    /**
     * 选主元 matrix->tmp_matrix
     */
    if (matrix->nnz < 100) { LOG_INFO("原矩阵:"), print_dense_matrix(csr2dense(matrix), matrix->num_row); }
    SparseMatrix *tmp_matrix;
    if (is_static_pivoting) {
        INDEX_TYPE *perm = NULL;
        INDEX_TYPE *iperm = NULL;
        ELE_TYPE *row_scale = NULL;
        ELE_TYPE *col_scale = NULL;
        static_pivoting(matrix, &perm, &iperm, &row_scale, &col_scale);
        tmp_matrix = init_csr_matrix(n, n, nnz);
        apply_pivoting(matrix, tmp_matrix, perm, row_scale, col_scale);
        info->mc64_perm = perm;
        info->mc64_iperm = iperm;
        info->Dr = row_scale;
        info->Dc = col_scale;
        // if (matrix->nnz < 100) {
        //     LOG_INFO("pivoting:"), print_matrix_csr(tmp_matrix->row_pointers, tmp_matrix->col_indices,
        //                                             tmp_matrix->num_row);
        //     printf("\nperm:");
        //     for (INDEX_TYPE i = 0; i < n; i++) {
        //         printf("%ld ", perm[i]);
        //     }
        //     printf("\nDr:");
        //     for (INDEX_TYPE i = 0; i < n; i++) {
        //         printf("%lf ", row_scale[i]);
        //     }
        //     printf("\nDc:");
        //     for (INDEX_TYPE i = 0; i < n; i++) {
        //         printf("%lf ", col_scale[i]);
        //     }
        //     printf("\n");
        // }
        if (matrix->nnz < 100) {
            LOG_INFO("pivoting values:"), print_dense_matrix(csr2dense(tmp_matrix), tmp_matrix->num_col);
        }
    } else {
        tmp_matrix = matrix;
    }
    //csr2mtx("mc64.mtx", tmp_matrix);
    //int_vector2mtx("mc64_perm.mtx", perm, n);
    /**
     * 重排序 tmp_matrix->A
     */
    SparseMatrix *A = init_csr_matrix(n, n, nnz);
    if (is_reorder) {
        SparseMatrix *temp2 = init_csr_matrix(n, n, A->nnz);

        reorder_csr_amd(tmp_matrix, info);
        apply_permutation(tmp_matrix, temp2, info->reorder_iperm);

        INDEX_TYPE *toposort_iperm = reorder_toposort(temp2, n);
        apply_permutation(temp2, A, toposort_iperm);
        info->toposort_iperm = toposort_iperm;

        if (is_static_pivoting) {
            free_sparse_matrix(tmp_matrix);
        }
        if (matrix->nnz < 100) { LOG_INFO("reorder:"), print_matrix_csr(A->row_pointers, A->col_indices, A->num_row); }
        if (matrix->nnz < 100) { LOG_INFO("reorder values:"), print_dense_matrix(csr2dense(A), A->num_row); }
    } else {
        A = tmp_matrix;
    }
    //    csr2mtx("reorder01.mtx", A);
    //    fp_vector2mtx("Dr.mtx", info->Dr, n);
    //    fp_vector2mtx("Dc.mtx", info->Dc, n);
    //    int_vector2mtx("amd_perm01.mtx", info->reorder_perm, n);
    //    int_vector2mtx("amd_iperm01.mtx", info->reorder_iperm, n);
    //    int_vector2mtx("mc64_perm01.mtx", perm, n);
    /**
     * 符号计算
     */
    if (is_symbolic_calc) {
        LOG_INFO("符号计算开始......");
        clock_t symbolic_calc_time = clock();
        //symbolic_calc_sym(A, info);
        csr2csc_pattern(A);
        info->L = init_sparse_matrix(A->num_row, A->num_row, A->nnz);
        info->U = init_sparse_matrix(A->num_row, A->num_row, A->nnz);
        //反着来的
        fill_in_2_no_sort_pruneL(A->num_row, A->nnz, A->col_indices, A->row_pointers,
                                 &info->U->row_pointers, &info->U->col_indices,
                                 &info->L->row_pointers, &info->L->col_indices);
        // if (matrix->nnz < 100) { LOG_INFO("L:"), print_matrix_csr(info->L->row_pointers, info->L->col_indices, A->num_row); }
        // if (matrix->nnz < 100) { LOG_INFO("U:"), print_matrix_csr(info->U->row_pointers, info->U->col_indices, A->num_row); }
        //pattern=L+U
        info->U->nnz = info->U->row_pointers[A->num_col] + n;
        info->L->nnz = info->L->row_pointers[A->num_col] - n;
        const INDEX_TYPE lu_nnz = info->U->nnz + info->L->nnz;
        LOG_DEBUG("l+u-e_nnz=%ld", lu_nnz);
        info->pattern = init_csr_matrix(A->num_row, A->num_col, lu_nnz);
        // info->pattern->col_pointers = (INDEX_TYPE *) lu_malloc((A->num_row + 1) * sizeof(INDEX_TYPE));
        // info->pattern->row_indices = (INDEX_TYPE *) lu_malloc(A->nnz * sizeof(INDEX_TYPE));
        INDEX_TYPE index_size = 0;
        info->pattern->row_pointers[0] = index_size;
        for (INDEX_TYPE i = 0; i < A->num_row; i++) {
            for (INDEX_TYPE j = info->L->row_pointers[i]; j < info->L->row_pointers[i + 1]; j++) {
                info->pattern->col_indices[index_size++] = info->L->col_indices[j];
            }
            for (INDEX_TYPE j = info->U->row_pointers[i]; j < info->U->row_pointers[i + 1]; j++) {
                info->pattern->col_indices[index_size++] = info->U->col_indices[j];
            }
            info->pattern->row_pointers[i + 1] = index_size;
        }
        info->U->nnz = info->U->row_pointers[A->num_col];
        csr2csc_pattern(info->U);
        info->L->nnz = info->L->row_pointers[A->num_col];
        csr2csc_pattern(info->L);
        csr2csc_pattern(info->pattern);
        get_diag_index(info->pattern, info);
        LOG_INFO("符号计算 elapsed time: %lf ms", ((double) (clock() - symbolic_calc_time)) / CLOCKS_PER_SEC * 1000.0);
    }
    if (matrix->nnz < 100) {
        LOG_INFO("L+U:"), print_matrix_csc(info->pattern->col_pointers, info->pattern->row_indices, A->num_row);
    }

    // matrix->csr_values = A->csr_values;
    // matrix->csc_values = A->csc_values;
    // matrix->row_pointers = A->row_pointers;
    // matrix->col_pointers = A->col_pointers;
    // matrix->row_indices = A->row_indices;
    // matrix->col_indices = A->col_indices;
    //lu_free(A);
    LOG_INFO("预处理 elapsed time: %lf ms", ((double) (clock() - preprocess_time)) / CLOCKS_PER_SEC * 1000.0);
    /**
     *计算矩阵结构信息
     */
    //calculate_asymmetry(matrix, info);
    return A;
}
