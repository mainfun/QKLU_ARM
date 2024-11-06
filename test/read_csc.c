#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
//
// Created by mainf on 2024/10/15.
//
#define INDEX_TYPE int
#define ELE_TYPE double

INDEX_TYPE *Ap;
INDEX_TYPE *Ai;
ELE_TYPE *Ax;

typedef struct {
    INDEX_TYPE rows;
    INDEX_TYPE cols;
    INDEX_TYPE not_zeros;
    char object[20];
    char format[20];
    char field[20];
    char symmetry[20];
} MtxInfo;


FILE *file_open(const char *file_name) {
    FILE *f;
    if ((f = fopen(file_name, "r")) == NULL) {
        printf("Could not open the file: %s", file_name);
    }
    return f;
}

void read_mtx_header(FILE *file, MtxInfo *info) {
    char buffer[1024];
    if (fgets(buffer, sizeof(buffer), file) == NULL) {
        printf("Error reading mtx file header");
    }
    // Parse the header to extract matrix information
    if (sscanf(buffer, "%%%%MatrixMarket %s %s %s %s", info->object, info->format, info->field, info->symmetry) != 4) {
        printf("mtx file header does not contain enough information");
    }

    // Read the header and skip comment lines
    while (fgets(buffer, sizeof(buffer), file)) {
        if (buffer[0] != '%') {
            // Check if the line is not a comment
            break;
        }
    }

    // Now buffer holds the line with matrix dimensions
    if (sscanf(buffer, "%lld %lld %lld", &(info->rows), &(info->cols), &(info->not_zeros)) != 3) {
        printf("Error reading matrix dimensions");
    }
}

void load_matrix_csc(const char *file_name) {
    FILE *f = file_open(file_name);
    INDEX_TYPE n_row, n_col, nz;

    MtxInfo mtxInfo;
    read_mtx_header(f, &mtxInfo);
    if (strcmp(mtxInfo.format, "coordinate") != 0) {
        printf("matrix in mtx file not is sparse");
    }

    nz = mtxInfo.not_zeros;
    n_row = mtxInfo.rows;
    n_col = mtxInfo.cols;

    double sparse_rate = (double) nz / (double) n_row / (double) n_col * 100;
    printf("matrix information: row=%lld, col=%lld, nz=%lld, sparse rate=%lf%%\n", n_row, n_col, nz, sparse_rate);
    Ap = (INDEX_TYPE *) malloc(sizeof(INDEX_TYPE) * (n_col + 1));
    Ai = (INDEX_TYPE *) malloc(sizeof(INDEX_TYPE) * nz);
    Ax = (ELE_TYPE *) malloc(sizeof(ELE_TYPE) * nz);

    memset(Ap, 0, (n_row + 1) * sizeof(INDEX_TYPE));

    for (INDEX_TYPE i = 0; i < nz; i++) {
        INDEX_TYPE row, col;
        ELE_TYPE value;
        fscanf(f, "%d %d %lg\n", &row, &col, &value);
        Ap[col]++;
        Ai[i] = row;
        Ax[i] = value;
    }
    Ap[0] = 1;
    for (INDEX_TYPE i = 1; i < n_row + 1; ++i) {
        Ap[i] += Ap[i - 1];
    }
    fclose(f);
}


int main() {
    INDEX_TYPE *ia = Ap, *ja = Ai;
    ELE_TYPE *a = Ax;
    load_matrix_csc("/Users/mainf/其他/mtx/onetone1.mtx");
    // load_matrix_csc("../res/test02.mtx");
    for (INDEX_TYPE i = 0; i < 5; ++i) {
        for (INDEX_TYPE j = Ap[i]; j < Ap[i + 1]; j++) {
            printf("%d    %d   %lf\n", Ai[j-1], i+1, Ax[j-1]);
        }
    }
    for (INDEX_TYPE i = 0; i < 5; ++i) {
        printf("%d, ", Ap[i]);
    }
    printf("\n");
    for (INDEX_TYPE i = 0; i < 10; ++i) {
        printf("%d, ", Ai[i]);
    }
    printf("\n");
    for (INDEX_TYPE i = 0; i < 5; ++i) {
        printf("%lf, ", Ax[i]);
    }
}
