//
// Created by mainf on 2024/10/13.
//

#include <layer_matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>

#define ELE_TYPE double
#define INDEX_TYPE long long
#define BLOCK_SIDE 5
#define A(i,j) a[(i)*BLOCK_SIDE+(j)]
#define B(i,j) b[(i)*BLOCK_SIDE+(j)]
#define C(i,j) c[(i)*BLOCK_SIDE+(j)]
#define DA(i,j) a[(i)*BLOCK_SIDE+(j)]
#define DB(i,j) b[(i)*BLOCK_SIDE+(j)]
#define DC(i,j) c[(i)*BLOCK_SIDE+(j)]

///用a解b
void DTSTRF_DD(const ELE_TYPE *a, ELE_TYPE *b) {
    for (int i = 0; i < BLOCK_SIDE; i++) {
        for (int j = 0; j < BLOCK_SIDE; j++) {
            //第j行消第i行
            ELE_TYPE l = B(i, j) /= A(j, j);
            for (int k = j + 1; k < BLOCK_SIDE; k++) {
                B(i, k) -= l * A(j, k);
            }
        }
    }
}

///用m1解m2
void DTSTRF_SS_v2(BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    for (int i = 0; i < BLOCK_SIDE; i++) {
        for (int j = m2->row_pointers[i]; j < m2->row_pointers[i + 1]; ++j) {
            int b_col = m2->col_indices[j];
            ELE_TYPE l = B(i, b_col) /= A(b_col, b_col);
            for (int k = m1->row_pointers[b_col] + 1; k < m1->row_pointers[b_col + 1]; k++) {
                int a_col = m1->col_indices[k];
                B(i, a_col) -= l * A(b_col, a_col);
            }
        }
    }
}

//U CSR L CSC
void DTSTRF_SS(BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    for (int i = 0; i < BLOCK_SIDE; i++) {
        for (int j = m2->col_pointers[i]; j < m2->col_pointers[i + 1]; ++j) {
            int b_row = m2->row_indices[j];
            ELE_TYPE l = B(b_row, i) /= A(i, i);
            for (int k = m1->row_pointers[i] + 1; k < m1->row_pointers[i + 1]; k++) {
                int a_col = m1->col_indices[k];
                B(b_row, a_col) -= l * A(i, a_col);
            }
        }
    }
}

//U稀疏
void DTSTRF_SD(BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    for (int i = 0; i < BLOCK_SIDE; i++) {
        for (int j = 0; j < BLOCK_SIDE; ++j) {
            B(j, i) /= A(i, i);
        }
        for (int k = m1->row_pointers[i] + 1; k < m1->row_pointers[i + 1]; k++) {
            int a_col = m1->col_indices[k];
            ELE_TYPE a_i_a_col = A(i, a_col);
            for (int j = 0; j < BLOCK_SIDE; ++j) {
                B(j, a_col) -= B(j, i) * a_i_a_col;
            }
        }
    }
}

//U稠密
void DTSTRF_DS(BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    for (int i = 0; i < BLOCK_SIDE; i++) {
        for (int j = m2->col_pointers[i]; j < m2->col_pointers[i + 1]; ++j) {
            int b_row = m2->row_indices[j];
            ELE_TYPE l = B(b_row, i) /= A(i, i);
            for (int k = i + 1; k < BLOCK_SIDE; k++) {
                B(b_row, k) -= l * A(i, k);
            }
        }
    }
}


int main() {
    int row_ptr[] = {0, 3, 5, 7, 8, 9}; // 每行非零元素的起始位置
    int col_indices[] = {0, 2, 4, 1, 3, 2, 4, 3, 4}; // 非零元素的列索引

    int col_ptr[] = {0, 1, 2, 4, 6, 9};
    int row_indices[] = {0, 1, 0, 2, 1, 3, 0, 2, 4};

    ELE_TYPE values[] = {
        1, 0, 6, 0, 7,
        0, 2, 0, 8, 0,
        0, 0, 3, 0, 9,
        0, 0, 0, 4, 0,
        0, 0, 0, 0, 5
        //1, 1, 2, 2, 3
        //0, 1, 2, 4, 6, 9
    }; // 非零元素的值
    BlockMatrix m1 = {
        .l_nnz = 6,
        .row_pointers = row_ptr,
        .col_indices = col_indices,
        .values = values,
    };
    int row_ptr2[] = {0, 2, 5, 8, 9, 11}; // 每行非零元素的起始位置
    int col_indices2[] = {1, 3, 1, 3, 4, 0, 2, 4, 3, 1, 3}; // 非零元素的列索引
    int col_ptr2[] = {0, 1, 4, 5, 9, 11};
    int row_indices2[] = {
        2,
        0, 1, 4,
        2,
        0, 1, 3, 4,
        1, 2
    };
    ELE_TYPE values2[] = {
        0, 1, 0, 2, 0, //0
        0, 3, 0, 0, 4, //1
        5, 0, 6, 0, 0, //2
        0, 0, 0, 7, 0, //3
        0, 8, 0, 0, 0, //4
        //1, 3, 1, 2, 1
        //0, 1, 4, 5, 7, 8
    }; // 非零元素的值

    BlockMatrix m2 = {
        .u_nnz = 8,
        .col_pointers = col_ptr2,
        .row_indices = row_indices2,
        .values = values2,
    };
    ELE_TYPE *original_values2 = (ELE_TYPE *) calloc(BLOCK_SIDE * BLOCK_SIDE, sizeof(ELE_TYPE));
    memcpy(original_values2, values2, BLOCK_SIDE * BLOCK_SIDE * sizeof(ELE_TYPE));

    // DTSTRF_DD(values, values2);
    // DTSTRF_SS(&m1, &m2);
    // DTSTRF_SD(&m1,&m2);
    DTSTRF_DS(&m1, &m2);
    print_dense_matrix(values2,BLOCK_SIDE);


    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                BLOCK_SIDE, BLOCK_SIDE, BLOCK_SIDE, 1.0,
                m2.values, BLOCK_SIDE, m1.values, BLOCK_SIDE, -1.0, original_values2, BLOCK_SIDE);
    printf("C:\n");
    print_dense_matrix(original_values2,BLOCK_SIDE);

    return EXIT_SUCCESS;
}
