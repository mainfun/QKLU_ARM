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

int main() {
    SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/tmt_unsym.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/onetone1.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("../res/test02.mtx", false);
    const INDEX_TYPE n = original_matrix->num_row;
    PreprocessInfo *info = init_preprocess_info();
    SparseMatrix *A = preprocess(original_matrix, info, true, true, true);
    L2Matrix *l2 = init_L2Matrix();
    //L2Matrix block_matrix_array2D[][];
    // csr2L2Matrix(A_pattern->row_pointers, A_pattern->col_indices, A, l2_matrix, BLOCK_SIDE);
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
    clock_t solving_time = clock();
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
    forward_substitution_block(l2, info->reorder_iperm, info->reorder_perm,
                               info->Dr, info->mc64_perm, info->mc64_iperm, b, y, n);
    backward_substitution_block(l2, info->reorder_perm, info->reorder_iperm, info->Dc, x, y, n);
    check_solving(original_matrix, x, b);


    //todo free L2
    return EXIT_SUCCESS;
}
