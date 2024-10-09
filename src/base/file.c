//
// Created by mainf on 2024/4/17.
//
#include "file.h"
#include "sort.h"
#include "log.h"

FILE *file_open(const char *file_name) {
    FILE *f;
    if ((f = fopen(file_name, "r")) == NULL) {
        LOG_ERROR("Could not open the file: %s", file_name);
    }
    return f;
}

void read_mtx_header(FILE *file, MtxInfo *info) {
    char buffer[1024];
    if (fgets(buffer, sizeof(buffer), file) == NULL) {
        LOG_ERROR("Error reading mtx file header");
    }
    // Parse the header to extract matrix information
    if (sscanf(buffer, "%%%%MatrixMarket %s %s %s %s", info->object, info->format, info->field, info->symmetry) != 4) {
        LOG_ERROR("mtx file header does not contain enough information");
    }

    // Read the header and skip comment lines
    while (fgets(buffer, sizeof(buffer), file)) {
        if (buffer[0] != '%') { // Check if the line is not a comment
            break;
        }
    }

    // Now buffer holds the line with matrix dimensions
    if (sscanf(buffer, "%lld %lld %lld", &(info->rows), &(info->cols), &(info->not_zeros)) != 3) {
        LOG_ERROR("Error reading matrix dimensions");
    }
}

/**
 * 如果列没有排好序，进行排序
 */
void adjust_csr(SparseMatrix *csr_matrix) {
    for (INDEX_TYPE r = 0; r < csr_matrix->num_row; ++r) {
        INDEX_TYPE row_start = csr_matrix->row_pointers[r];
        INDEX_TYPE row_end = csr_matrix->row_pointers[r + 1];
        INDEX_TYPE n = row_end - row_start;
        //todo 排序优化
        bubbleSort(csr_matrix->col_indices + row_start, csr_matrix->csr_values + row_start, n);
    }
}

SparseMatrix *add_transpose_for_sym(INDEX_TYPE *Ap, INDEX_TYPE *Ai, ELE_TYPE *Ax, INDEX_TYPE n, INDEX_TYPE nz);

bool is_sorted_ascending(const SparseMatrix *m);

SparseMatrix *load_matrix_csr(const char *file_name, bool base_0) {
    FILE *f = file_open(file_name);
    INDEX_TYPE n_row, n_col, nz;

    MtxInfo mtxInfo;
    read_mtx_header(f, &mtxInfo);
    if (strcmp(mtxInfo.format, "coordinate") != 0) {
        LOG_WARN("matrix in mtx file not is sparse");
    }

    #ifndef ELE_TYPE_COMPLEX
    if (strcmp(mtxInfo.field, "complex") == 0) {
        LOG_ERROR("complex");
    }
    #endif
    nz = mtxInfo.not_zeros;
    n_row = mtxInfo.rows;
    n_col = mtxInfo.cols;

    double sparse_rate = (double) nz / (double) n_row / (double) n_col * 100;
    LOG_DEBUG("matrix information: row=%lld, col=%lld, nz=%lld, sparse rate=%lf%%", n_row, n_col, nz, sparse_rate);

    SparseMatrix *csr_matrix = init_csr_matrix(n_row, n_col, nz);
    //csr_matrix->sparse_rate = sparse_rate;
    memset(csr_matrix->row_pointers, 0, (n_row + 1) * sizeof(INDEX_TYPE));

    INDEX_TYPE *file_row_index_arr = (INDEX_TYPE *) lu_malloc(nz * sizeof(INDEX_TYPE));
    INDEX_TYPE *file_col_index_arr = (INDEX_TYPE *) lu_malloc(nz * sizeof(INDEX_TYPE));
    ELE_TYPE *file_value_arr = (ELE_TYPE *) lu_malloc(nz * sizeof(ELE_TYPE));

    for (INDEX_TYPE i = 0; i < nz; i++) {
        INDEX_TYPE row, col;
        ELE_TYPE value;
        //todo 复数
        fscanf(f, "%lld %lld %lg\n", &row, &col, &value);
        // 0-based
        if (base_0) {
            row++;
            col++;
        }
        csr_matrix->row_pointers[row]++;
        // adjust from 1-based to 0-based
        row--;
        col--;
        file_row_index_arr[i] = row;
        file_col_index_arr[i] = col;
        file_value_arr[i] = value;
    }

    for (INDEX_TYPE i = 1; i < n_row + 1; ++i) {
        csr_matrix->row_pointers[i] += csr_matrix->row_pointers[i - 1];
    }

    INDEX_TYPE *column_index_count = lu_calloc(n_row + 1, sizeof(INDEX_TYPE));
    for (INDEX_TYPE i = 0; i < nz; ++i) {
        INDEX_TYPE row = file_row_index_arr[i];
        INDEX_TYPE col = file_col_index_arr[i];
        INDEX_TYPE offset = csr_matrix->row_pointers[row] + column_index_count[row]++;
        csr_matrix->col_indices[offset] = col;
        csr_matrix->csr_values[offset] = file_value_arr[i];
    }
    lu_free(column_index_count);
    lu_free(file_row_index_arr);
    lu_free(file_col_index_arr);
    lu_free(file_value_arr);
    fclose(f);
    if (strcmp(mtxInfo.symmetry, "symmetric") == 0) {
        csr_matrix = add_transpose_for_sym(csr_matrix->row_pointers,
                                           csr_matrix->col_indices,
                                           csr_matrix->csr_values,
                                           n_row, nz);
    }
    if (!is_sorted_ascending(csr_matrix)) adjust_csr(csr_matrix);
    return csr_matrix;
}

// 检查列号是否按升序排列
bool is_sorted_ascending(const SparseMatrix *m) {
    for (INDEX_TYPE i = 0; i < m->num_row; ++i) {
        for (INDEX_TYPE j = m->row_pointers[i]; j < m->row_pointers[i + 1] - 1; ++j) {
            if (m->col_indices[j] > m->col_indices[j + 1]) {
                return false;
            }
        }
    }
    return true;
}


SparseMatrix *load_matrix_csc_order(const char *file_name, bool base_0) {
    FILE *f = file_open(file_name);
    INDEX_TYPE n_row, n_col, nz;

    MtxInfo mtxInfo;
    read_mtx_header(f, &mtxInfo);
    if (strcmp(mtxInfo.format, "coordinate") != 0) {
        LOG_WARN("matrix in mtx file not is sparse");
    }

    #ifndef ELE_TYPE_COMPLEX
    if (strcmp(mtxInfo.field, "complex") == 0) {
        LOG_ERROR("complex");
    }
    #endif
    nz = mtxInfo.not_zeros;
    n_row = mtxInfo.rows;
    n_col = mtxInfo.cols;

    double sparse_rate = (double) nz / (double) n_row / (double) n_col * 100;
    LOG_DEBUG("matrix information: row=%lld, col=%lld, nz=%lld, sparse rate=%lf%%", n_row, n_col, nz, sparse_rate);

    SparseMatrix *csc_matrix = init_csc_matrix(n_row, n_col, nz);
    //csr_matrix->sparse_rate = sparse_rate;
    memset(csc_matrix->col_pointers, 0, (n_row + 1) * sizeof(INDEX_TYPE));

    for (INDEX_TYPE i = 0; i < nz; i++) {
        INDEX_TYPE row, col;
        ELE_TYPE value;
        //todo 复数
        fscanf(f, "%lld %lld %lg\n", &row, &col, &value);
        // 0-based
        if (base_0) {
            row++;
            col++;
        }
        csc_matrix->col_pointers[col]++;
        // adjust from 1-based to 0-based
        row--;
        col--;
        csc_matrix->row_indices[i] = row;
        csc_matrix->csc_values[i] = value;
    }

    for (INDEX_TYPE i = 1; i < n_row + 1; ++i) {
        csc_matrix->col_pointers[i] += csc_matrix->col_pointers[i - 1];
    }
    fclose(f);
    if (strcmp(mtxInfo.symmetry, "symmetric") == 0) {
        csc_matrix = add_transpose_for_sym(csc_matrix->row_pointers,
                                           csc_matrix->col_indices,
                                           csc_matrix->csr_values,
                                           n_row, nz);
    }
    return csc_matrix;
}

// C = A + A^T，其中A是三角矩阵
SparseMatrix *add_transpose_for_sym(INDEX_TYPE *Ap, INDEX_TYPE *Ai, ELE_TYPE *Ax, INDEX_TYPE n, INDEX_TYPE nz) {
    SparseMatrix *C= init_csr_matrix(n, n, 2 * nz);  // 最多可能有 2 * A.nnz 个非零元素

    int *row_nnz = (int *) calloc(n, sizeof(int));  // 每行的非零元素计数

    // 遍历A的每一行
    for (INDEX_TYPE i = 0; i < n; i++) {
        for (INDEX_TYPE j = Ap[i]; j < Ap[i + 1]; j++) {
            INDEX_TYPE col = Ai[j];
            ELE_TYPE val = Ax[j];

            // 添加A(i, j)
            C->col_indices[C->row_pointers[i] + row_nnz[i]] = col;
            C->csr_values[C->row_pointers[i] + row_nnz[i]] = val;
            row_nnz[i]++;

            // 添加A^T(i, j) -> A(j, i)
            if (i != col) {  // 防止对角线元素重复计算
                C->col_indices[C->row_pointers[col] + row_nnz[col]] = i;
                C->csr_values[C->row_pointers[col] + row_nnz[col]] = val;
                row_nnz[col]++;
            }
        }
    }

    // 更新C的row_ptr
    for (INDEX_TYPE i = 1; i <= n; i++) {
        C->row_pointers[i] = C->row_pointers[i - 1] + row_nnz[i - 1];
    }
    free(Ap);
    free(Ai);
    free(Ax);
    // 清理临时数组
    free(row_nnz);
    return C;
}

void csr2mtx(const char *file_name, const SparseMatrix *A) {
    INDEX_TYPE *Ap = A->row_pointers;
    INDEX_TYPE *Ai = A->col_indices;
    ELE_TYPE *Ax = A->csr_values;
    FILE *f;
    if ((f = fopen(file_name, "w")) == NULL) {
        LOG_ERROR("Could not write the file: %s", file_name);
    }
    fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%lld %lld %lld\n", A->num_row, A->num_col, A->nnz);
    for (INDEX_TYPE r = 0; r < A->num_row; ++r) {
        for (INDEX_TYPE j = Ap[r]; j < Ap[r + 1]; ++j) {
            INDEX_TYPE c = Ai[j];
            ELE_TYPE v = Ax[j];
            fprintf(f, "%lld %lld %.100lg\n", r, c, v);
        }
    }
    fclose(f);
}

void csc2mtx(const char *file_name, const SparseMatrix *A) {
    INDEX_TYPE *Ap = A->col_pointers;
    INDEX_TYPE *Ai = A->row_indices;
    ELE_TYPE *Ax = A->csc_values;
    FILE *f;
    if ((f = fopen(file_name, "w")) == NULL) {
        LOG_ERROR("Could not write the file: %s", file_name);
    }
    fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%lld %lld %lld\n", A->num_row, A->num_col, A->nnz);
    for (INDEX_TYPE i = 0; i < A->num_row; ++i) {
        for (INDEX_TYPE j = Ap[i]; j < Ap[i + 1]; ++j) {
            INDEX_TYPE r = Ai[j];
            ELE_TYPE v = Ax[j];
            fprintf(f, "%lld %lld %.100lg\n", r, i, v);
        }
    }
    fclose(f);
}

void L2mtx(const char *file_name, const SparseMatrix *A) {
    INDEX_TYPE *Ap = A->col_pointers;
    INDEX_TYPE *Ai = A->row_indices;
    ELE_TYPE *Ax = A->csc_values;
    FILE *f;
    if ((f = fopen(file_name, "w")) == NULL) {
        LOG_ERROR("Could not write the file: %s", file_name);
    }
    fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%lld %lld %lld\n", A->num_row, A->num_col, A->nnz);
    for (INDEX_TYPE c = 0; c < A->num_row; ++c) {
        fprintf(f, "%lld %lld %lf\n", c, c, 1.0);
        for (INDEX_TYPE j = Ap[c]; j < Ap[c + 1]; ++j) {
            INDEX_TYPE r = Ai[j];
            ELE_TYPE v = Ax[j];
            fprintf(f, "%lld %lld %.100lg\n", r, c, v);
        }
    }
    fclose(f);
}

void int_vector2mtx(const char *file_name, const INDEX_TYPE *vector, INDEX_TYPE n) {
    FILE *f;
    if ((f = fopen(file_name, "w")) == NULL) {
        LOG_ERROR("Could not write the file: %s", file_name);
    }
    fprintf(f, "%lld\n", n);
    for (INDEX_TYPE i = 0; i < n; ++i) {
        fprintf(f, "%lld\n", vector[i]);
    }
    fclose(f);
}

void fp_vector2mtx(const char *file_name, const ELE_TYPE *vector, INDEX_TYPE n) {
    FILE *f;
    if ((f = fopen(file_name, "w")) == NULL) {
        LOG_ERROR("Could not write the file: %s", file_name);
    }
    fprintf(f, "%lld\n", n);
    for (INDEX_TYPE i = 0; i < n; ++i) {
        fprintf(f, "%.17lg\n", vector[i]);
    }
    fclose(f);
}