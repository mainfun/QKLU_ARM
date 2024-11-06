#include <solving.h>
#include <stdio.h>
#include <symbolic_analysis.h>
#include <base/plot.h>

#include "check.h"
#include "layer_matrix.h"
#include "numerical.h"
#include "preprocess.h"
#include "base/identify_chip.h"
#include "base/file.h"

int main(const int argc, char *argv[]) {
    if (argc < 2) { LOG_ERROR("Usage: %s <matrix_file>", argv[0]); }
    // SparseMatrix *m = load_matrix_csr("/Users/mainf/其他/mtx/g7jac180.mtx", false);
    // SparseMatrix *m = load_matrix_csr("./res/test02.mtx", false);
    // SparseMatrix *m = load_matrix_csc_order("./res/test01.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/tmt_unsym.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr(argv[1], false);
    SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/onetone1.mtx", false);
    clock_t start_main = clock();
    INDEX_TYPE n = original_matrix->num_row;
    PreprocessInfo *info = init_preprocess_info();
    SparseMatrix *A = preprocess(original_matrix, info, true, true, true);
    SparseMatrix *A_pattern = info->pattern;
    //csr2image(A_pattern,"onetone1.jpg",2000,2000);
    //分块
    L2Matrix *l2_matrix = init_L2Matrix();
    INDEX_TYPE block_size = 100;
    //csr2L2Matrix(A_pattern->row_pointers, A_pattern->col_indices, A, l2_matrix, block_size);
    //end 分块
    ELE_TYPE *Lx = lu_malloc(info->L->nnz * sizeof(ELE_TYPE));
    ELE_TYPE *Ux = lu_malloc((info->U->nnz + original_matrix->num_row) * sizeof(ELE_TYPE));
    factor(A, Lx, Ux,
           info->L->col_pointers, info->L->col_pointers + 1, info->L->row_indices,
           info->U->row_pointers, info->U->row_pointers + 1, info->U->col_indices);

    info->L->csc_values = Lx;
    info->U->csr_values = Ux;
    lu_free(info->L->row_pointers);
    lu_free(info->L->col_indices);
    clock_t solving_time = clock();
    csc2csr(info->L);
    csr2csc(info->U);
    // LOG_INFO("L:");
    // print_dense_matrix(csr2dense(info->L),n);
    // LOG_INFO("U:");
    // print_dense_matrix(csr2dense(info->U),n);
    //check_lu(A, info->L, info->U);
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
    forward_substitution(info->L, info->reorder_iperm, info->reorder_perm,
                         info->Dr, info->mc64_perm, info->mc64_iperm, b, y, n);
    backward_substitution(info->U, info->reorder_perm, info->reorder_iperm, info->Dc, x, y, n);
    // LOG_INFO("X:");
    // for (int i = 0; i < n; ++i) {
    //     printf("%lf ",x[i]);
    // }
    check_solving(original_matrix, x, b);
    lu_free(x);
    lu_free(y);
    free_sparse_matrix(A);
    free_sparse_matrix(original_matrix);
    // free_preprocess_info(info);
    LOG_INFO("solving elapsed time: %lf ms", ((double) (clock() - solving_time)) / CLOCKS_PER_SEC * 1000.0);
    LOG_INFO("main() elapsed time: %lf ms", ((double) (clock() - start_main)) / CLOCKS_PER_SEC * 1000.0);
    return EXIT_SUCCESS;
}
