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
#include <sys/time.h>

#include "symbol_calc.h"
#include "part3_solver.h"

#include "toposort.h"

#define TIME_FUNCTION(func_call) ({ \
struct timeval start, end; \
gettimeofday(&start, NULL); \
func_call; \
gettimeofday(&end, NULL); \
double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6; \
printf("Function %s executed in %f ms\n", #func_call, elapsed_time*1000); \
elapsed_time; \
})

SparseMatrix *blocking(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai,
                       INDEX_TYPE row_start, INDEX_TYPE row_end,
                       INDEX_TYPE col_start, INDEX_TYPE col_end) {
    double time = omp_get_wtime();
    INDEX_TYPE n = row_end - row_start;
    INDEX_TYPE nnz = 0;
    for (INDEX_TYPE i = row_start; i < row_end; i++) {
        for (INDEX_TYPE p = Ap[i]; p < Ap[i + 1]; p++) {
            INDEX_TYPE col_index = Ai[p];
            if (col_start <= col_index && col_index < col_end) {
                nnz++;
            }
        }
    }
    SparseMatrix *r = init_csr_matrix(n, n, nnz);
    INDEX_TYPE top = 0;
    r->row_pointers[0] = 0;
    for (INDEX_TYPE i = row_start; i < row_end; i++) {
        for (INDEX_TYPE p = Ap[i]; p < Ap[i + 1]; p++) {
            INDEX_TYPE col_index = Ai[p];
            if (col_start <= col_index && col_index < col_end) {
                r->col_indices[top++] = col_index - col_start;
            }
        }
        r->row_pointers[i - row_start + 1] = top;
    }
    LOG_INFO("blocking time is %lf ms", (omp_get_wtime()-time)*1000);
    return r;
}

int main() {
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/tmt_unsym.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/onetone1.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/k3plates.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/g7jac180.mtx", false);
    SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/ASIC_680k.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/apache2.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("../res/test02.mtx", false);
    const INDEX_TYPE n = original_matrix->num_row;
    double main_time = omp_get_wtime();
    PreprocessInfo *info = init_preprocess_info();
    SparseMatrix *A = preprocess(original_matrix, info, true, true, true);
    SparseMatrix *r = blocking(info->pattern->row_pointers, info->pattern->col_indices,
                               info->cut_point, n, info->cut_point, n);
    // SparseMatrix *r = blocking(info->pattern->row_pointers, info->pattern->col_indices,
    //          info->cut_point, n, info->cut_point, n);
    // INDEX_TYPE _=0;
    // INDEX_TYPE *toposort_iperm = reorder_toposort(r,n-info->cut_point,&_);
    // SparseMatrix *rr = init_csr_matrix(r->num_row, r->num_row, r->nnz);
    // apply_permutation(r, rr, toposort_iperm);
    // csr2image(rr,"tmt_unsym_block_v0.jpg",2000,2000);
    // exit(0);
    // INDEX_TYPE *rp_l = lu_malloc((n + 1) * sizeof(INDEX_TYPE));
    // INDEX_TYPE *rp_u = lu_malloc((n + 1) * sizeof(INDEX_TYPE));
    // INDEX_TYPE *ri_l = lu_malloc((info->L->nnz) * sizeof(INDEX_TYPE));
    // INDEX_TYPE *ri_u = lu_malloc((info->U->nnz) * sizeof(INDEX_TYPE));
    // INDEX_TYPE l_nnz, u_nnz;
    //static_symbol_calc(A, &l_nnz, &u_nnz, rp_l, ri_l, rp_u, ri_u);
    // csr2_RGB_image(info->pattern,"tmt_unsym_topo_v3.jpg",1720,1720);
    // csr2_RGB_image(info->pattern,"g7jac180_topo_v3.jpg",1067,1067);
    // csr2_RGB_image(info->pattern,"onetone1_topo_v5.jpg",750,750);
    // csr2image_block(info->pattern->row_pointers, info->pattern->col_indices,
    // "g7jac180_rb_que.jpg", 38912, 53370, 38912, 53370);
    // csr2image_block(A->row_pointers, A->col_indices,
    //                 "tmt_unsym_rb_que_A.jpg", 854567, 917825, 854567, 917825);
    // csr2image_block(info->pattern->row_pointers, info->pattern->col_indices,
    //                 "onetone1_rb_stack_sym.jpg", 32381, 36057, 32381, 36057);
    // csr2image(info->pattern,"onetone1_topo_v4.jpg",720,720);
    // csr2image_block(info->pattern->row_pointers, info->pattern->col_indices,
    //                 "g7jac180_10%2.jpg", 0, 5337, 53370*0.9, 53370);
    //csr2image(info->pattern,"tmt_unsym_topo500.jpg",2000,2000);
    L2Matrix *l2 = init_L2Matrix();
    calc_block_side(info->L->row_pointers, info->L->col_indices, l2, BLOCK_SIDE, n);
    //csr2_RGB_image_coor(info->pattern,"g7jac180——coor.jpg",1835,1835);
    // csr2image_block_cut(A->row_pointers, A->col_indices,
    //                     "onetone1_cut——1000.jpg",
    //                     0, 1000,
    //                     0, 1000,
    //                     l2->block_side_sum, 30);
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
    // block_parallel_factor(l2, (long) (n - info->cut_point) / 50);

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
    TIME_FUNCTION(umfpack_solver(r,b,x));
    double solving_time = omp_get_wtime();
    forward_substitution_block(l2, info->reorder_iperm, info->reorder_perm,
                               info->Dr, info->mc64_perm, info->mc64_iperm, b, y, n,
                               info->toposort_iperm);
    backward_substitution_block(l2, info->reorder_perm, info->reorder_iperm, info->Dc, x, y, n,
                                info->toposort_iperm);
    check_solving(original_matrix, x, b);
    LOG_INFO("solving_time is %lf ms", (omp_get_wtime()-solving_time)*1000);
    LOG_INFO("main time is %lf ms", (omp_get_wtime()-main_time)*1000);
    //todo free L2
    return EXIT_SUCCESS;
}
