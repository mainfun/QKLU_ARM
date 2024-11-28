//
// Created by mainf on 2024/10/10.
//
#include <layer_matrix.h>
#include <preprocess.h>
#include <base/file.h>
#include <cblas.h>
#include <check.h>
#include <solving.h>
#include <symbolic_analysis.h>
#include <numerical.h>
#include <omp.h>
#include <reorder.h>
#include <base/plot.h>

#include "toposort.h"

int main() {
    SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/tmt_unsym.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/onetone1.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/k3plates.jpg.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("../res/test02.mtx", false);
    const INDEX_TYPE n = original_matrix->num_row;
    PreprocessInfo *info = init_preprocess_info();
    SparseMatrix *A = preprocess(original_matrix, info, true, true, true);
    //csr2image(info->pattern,"tmt_unsym_topo.jpg",2000,2000);
    L2Matrix *l2 = init_L2Matrix();
    calc_block_side(info->L->col_pointers, info->L->row_indices, l2, BLOCK_SIDE, n);
    //csr2_RGB_image(info->pattern,"tmt_unsym.jpg",18356,18356);
    csr2L2Matrix(info->L->row_pointers, info->L->col_indices,
                 info->U->row_pointers, info->U->col_indices,
                 A->row_pointers, A->col_indices,
                 A->csr_values, l2, BLOCK_SIDE, n);
    //---------------------------l2 get csc---------------------------
    l2->col_pointers = (int *) lu_calloc(l2->num_col_block + 1, sizeof(int));
    l2->row_indices = (int *) lu_malloc(l2->block_count * sizeof(int));
    csr2csc_pattern_v2(l2->row_pointers, l2->col_indices,
                       l2->col_pointers, l2->row_indices, l2->block_count, l2->num_row_block);
    l2->diag_index_csc = get_diag_index_i(l2->col_pointers, l2->row_indices, l2->num_row_block);
    //-----------------------end l2 get csc end-------------------------
    block_factor(l2);
    // block_parallel_factor(l2);
    //初始化解向量b
    ELE_TYPE *b = (ELE_TYPE *) malloc(n * sizeof(ELE_TYPE));
    //    random_vector(b, n);
    for (INDEX_TYPE i = 0; i < n; ++i) {
        ELE_TYPE sum = 0;
        for (INDEX_TYPE j = original_matrix->row_pointers[i]; j < original_matrix->row_pointers[i + 1]; ++j) {
            sum += original_matrix->csr_values[j];
        }
        b[i] = sum;
    }
    ELE_TYPE *x = (ELE_TYPE *) malloc(n * sizeof(ELE_TYPE));
    ELE_TYPE *y = (ELE_TYPE *) malloc(n * sizeof(ELE_TYPE));
    double solving_time = omp_get_wtime();
    forward_substitution_block(l2, info->reorder_iperm, info->reorder_perm,
                               info->Dr, info->mc64_perm, info->mc64_iperm, b, y, n,
                               info->toposort_iperm);
    backward_substitution_block(l2, info->reorder_perm, info->reorder_iperm, info->Dc, x, y, n,
                                info->toposort_iperm);
    check_solving(original_matrix, x, b);
    LOG_INFO("solving_time is %lf ms", (omp_get_wtime()-solving_time)*1000);

    //todo free L2
    return EXIT_SUCCESS;
}
