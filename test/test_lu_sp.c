//
// Created by mainf on 2024/10/14.
//
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <layer_matrix.h>
#include <math.h>

#define ELE_TYPE double
#define BLOCK_SIDE 5
#define A(r,c) a[(r)*BLOCK_SIDE+(c)]
#define B(r,c) b[(r)*BLOCK_SIDE+(c)]
#define THRESHOLD 1e-8

void LUS_v2(const int *Lp_start, const int *Lp_end, const int *Li,
         const int *Up_start, const int *Up_end, const int *Ui,
         ELE_TYPE *a, int n) {
    for (int i = 0; i < n; ++i) {
        for (int p = Lp_start[i]; p < Lp_end[i]; p++) {
            int j = Li[p];
            //printf("第%lld行消第%lld行\n", j, i);
            ELE_TYPE l = A(i, j) /= A(j, j);
            ELE_TYPE *pivot_row_ptr = a + j * n;
            ELE_TYPE *eli_row_ptr = a + i * n;
            for (int k = Up_start[j] + 1; k < Up_end[j]; k++) {
                int c = Ui[k];
                //printf("eli_row_ptr[%d]-= %lf * pivot_row_ptr[%d]\n", c, l, c);
                eli_row_ptr[c] -= l * pivot_row_ptr[c];
                //a(i,c)-=l(i,j)*u(j,c)
            }
        }
    }
}

void LUS(const int *Lp_start, const int *Lp_end, const int *Li,
         const int *Up_start, const int *Up_end, const int *Ui,
         ELE_TYPE *a, int n) {
    //csr2csc
    for (int i = 0; i < n; ++i) {
        for (int p = Lp_start[i]; p < Lp_end[i]; p++) {
            int j = Li[p];
            printf("第%d行消第%d行\n", i, j);
            ELE_TYPE l = A(j, i) /= A(i, i);
            ELE_TYPE *pivot_row_ptr = a + i * n;
            ELE_TYPE *eli_row_ptr = a + j * n;
            for (int k = Up_start[i] + 1; k < Up_end[i]; k++) {
                int c = Ui[k];
                printf("eli_row_ptr[%d]-= %lf * pivot_row_ptr[%d]\n", c, l, c);
                eli_row_ptr[c] -= l * pivot_row_ptr[c];
            }
        }
    }
}

void LUD(ELE_TYPE *a) {
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        ELE_TYPE pivot = A(i, i);
        for (int j = i + 1; j < BLOCK_SIDE; ++j) {
            ELE_TYPE l = A(j, i) /= pivot;
            for (int k = i + 1; k < BLOCK_SIDE; ++k) {
                A(j, k) -= l * A(i, k);
            }
        }
    }
}

void verify_lu(const ELE_TYPE *a, const ELE_TYPE *original_matrix) {
    ELE_TYPE *b = (ELE_TYPE *) malloc(BLOCK_SIDE * BLOCK_SIDE * sizeof(ELE_TYPE));
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int k = i; k < BLOCK_SIDE; ++k) {
            B(i, k) = A(i, k);
        }
    }
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        //算B[i]
        for (int j = 0; j < i; ++j) {
            for (int k = j; k < BLOCK_SIDE; ++k) {
                //B[i][k]+=L[i][j]*U[j][k]
                B(i, k) += A(i, j) * A(j, k);
            }
        }
    }
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int k = 0; k < BLOCK_SIDE; ++k) {
            B(i, k) -= original_matrix[i * BLOCK_SIDE + k];
        }
    }
    printf("B:\n");
    print_dense_matrix(b,BLOCK_SIDE);
    free(b);
}


int main() {
    int l_row_ptr[] = {0, 0, 0, 1, 3, 4};
    int l_col_indices[] = {0, 0, 2, 1};

    int u_row_ptr[] = {0, 3, 4, 7, 9, 10};
    int u_col_indices[] = {0, 3, 4, 1, 2, 3, 4, 3, 4, 4};

    ELE_TYPE values[] = {
        10, 0, 0, 6, 7,
        0, 20, 0, 0, 0,
        8, 0, 30, 0, 12,
        9, 0, 11, 40, 0,
        0, 10, 0, 0, 50
    };
    ELE_TYPE original[] = {
        10, 0, 0, 6, 7,
        0, 20, 0, 0, 0,
        8, 0, 30, 0, 12,
        9, 0, 11, 40, 0,
        0, 10, 0, 0, 50
    };

    LUS(l_row_ptr, l_row_ptr + 1, l_col_indices, u_row_ptr, u_row_ptr + 1, u_col_indices, values,BLOCK_SIDE);
    // LUD(values);
    print_dense_matrix(values,BLOCK_SIDE);
    verify_lu(values, original);
    return 0;
}
