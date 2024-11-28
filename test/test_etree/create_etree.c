//
// Created by mainf on 2024/11/27.
//
#include <layer_matrix.h>
#include <preprocess.h>
#include <reorder.h>
#include <base/file.h>
#include <symbolic_analysis.h>
#include <base/plot.h>

#include "toposort.h"


int main() {
    SparseMatrix *om = load_matrix_csr("/Users/mainf/其他/mtx/fs_541_1.mtx", false);
    const INDEX_TYPE n = om->num_row;
    PreprocessInfo *info = init_preprocess_info();
    SparseMatrix *A = preprocess(om, info, true, true, true);
    //csr2image(A,"k3plates.jpg",11107,11107);
    INDEX_TYPE *iorder = reorder_toposort(A, n);
    SparseMatrix *order_matrix = init_csr_matrix(n, n, om->nnz);
    apply_permutation(A, order_matrix, iorder);


    csr2image(order_matrix, "fs541_toposort.jpg", 541, 541);
}
