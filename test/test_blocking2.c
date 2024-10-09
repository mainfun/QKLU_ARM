#include <stdio.h>
#include "base/identify_chip.h"
#include "base/file.h"

void print_csr(SparseMatrix *mat) {
    printf("%lld x %lld\n", mat->num_row, mat->num_col);
    printf("Row ptr: ");
    for (LU_INT i = 0; i <= mat->num_row; i++) {
        printf("%lld ", mat->row_pointers[i]);
    }
    printf("\nCol indices: ");
    for (LU_INT i = 0; i < mat->nnz; i++) {
        printf("%lld ", mat->col_indices[i]);
    }
    printf("\nValues: ");
    for (LU_INT i = 0; i < mat->nnz; i++) {
        printf("%lf ", mat->csr_values[i]);
    }
    printf("\n");
}

void print_block(SparseMatrix *mat) {
    for (int i = 0; i < mat->num_row_block; i++) {
        for (int j = 0; j < mat->num_col_block; j++) {
            printf("Submatrix (%d, %d):\n", i, j);
            print_csr(mat->sub_matrices[i * mat->num_col_block + j]);
        }
    }
}

int main() {
//    log_cpu_info();
////    SparseMatrix *m = load_matrix_csr("/Users/mainf/CLionProjects/my_lu/res/k3plates.mtx", false);
////    blocking_csr_matrix(m, 100, 100);
    SparseMatrix *m = load_matrix_csr("../res/test02.mtx", false);
    blocking_csr_matrix(m, 3, 3);
    print_block(m);
    free_sparse_matrix(m);
    return 0;
}