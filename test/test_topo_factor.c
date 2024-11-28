//
// Created by mainf on 2024/11/28.
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
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/tmt_unsym.mtx", false);
    SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/onetone1.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("/Users/mainf/其他/mtx/k3plates.jpg.mtx", false);
    // SparseMatrix *original_matrix = load_matrix_csr("../res/test02.mtx", false);
    const INDEX_TYPE n = original_matrix->num_row;
    PreprocessInfo *info = init_preprocess_info();
    SparseMatrix *A = preprocess(original_matrix, info, true, true, true);





}