//
// Created by mainf on 2024/11/30.
//

#include <symbol_calc.h>
#include <base/matrix.h>

int main() {
    INDEX_TYPE n = 6;
    INDEX_TYPE ptr[] = {0, 3, 5, 9, 12, 14, 16};
    INDEX_TYPE idx[] = {0, 2, 5, 0, 1, 0, 2, 3, 5, 0, 3, 4, 3, 4, 1, 5};
    SparseMatrix *A = init_sparse_matrix(6, 6, 16);
    A->row_pointers = ptr;
    A->col_indices = idx;
    csr2csc_pattern(A);
    INDEX_TYPE l_nnz = 30, u_nnz = 30;
    INDEX_TYPE *rp_l = lu_malloc((n + 1) * sizeof(INDEX_TYPE));
    INDEX_TYPE *rp_u = lu_malloc((n + 1) * sizeof(INDEX_TYPE));
    INDEX_TYPE *ri_l = lu_malloc(l_nnz * sizeof(INDEX_TYPE));
    INDEX_TYPE *ri_u = lu_malloc(u_nnz * sizeof(INDEX_TYPE));
    static_symbol_calc(A, &l_nnz, &u_nnz, rp_l, ri_l, rp_u, ri_u);
    for (INDEX_TYPE i = 0; i < l_nnz; ++i) {
        printf("%lld ", ri_l[i]);
    }
    printf("\n");
    for (INDEX_TYPE i = 0; i < u_nnz; ++i) {
        printf("%lld ", ri_u[i]);
    }
}
