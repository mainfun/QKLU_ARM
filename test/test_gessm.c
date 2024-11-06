//
// Created by mainf on 2024/10/12.
//
#include <layer_matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <cblas.h>

#define ELE_TYPE double
#define INDEX_TYPE int
#define BLOCK_SIDE 5
#define A(i,j) a[(i)*BLOCK_SIDE+(j)]
#define B(i,j) b[(i)*BLOCK_SIDE+(j)]
#define C(i,j) c[(i)*BLOCK_SIDE+(j)]
#define L(i,j) lx[(i)*BLOCK_SIDE+(j)]
#define U(i,j) ux[(i)*BLOCK_SIDE+(j)]
#define DA(i,j) a[(i)*BLOCK_SIDE+(j)]
#define DB(i,j) b[(i)*BLOCK_SIDE+(j)]
#define DC(i,j) c[(i)*BLOCK_SIDE+(j)]
#define DL(i,j) lx[(i)*BLOCK_SIDE+(j)]
#define DU(i,j) ux[(i)*BLOCK_SIDE+(j)]

///用m1解m2 L是CSC的 U csr
void DGESSM_SS(const int *Lp_start, const int *Lp_end, BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *b = m1->values;
    ELE_TYPE *a = m2->values;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int p = Lp_start[i]; p < Lp_end[i]; p++) {
            int j = m1->row_indices[p];
            ELE_TYPE l = B(j, i);
            //printf("第%d行消%d行\n", i, j);
            for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                int c = m2->col_indices[k];
                //printf("%lf -= %lf * %lf\n", A(j, c), l, A(i, c));
                A(j, c) -= l * A(i, c);
            }
        }
    }
}

///用m1解m2 L是CSR的 U csr
void DGESSM_SS2(const int *Lp_start, const int *Lp_end, BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *a = m2->values;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int p = Lp_start[i]; p < Lp_end[i]; p++) {
            int j = m1->col_indices[p];
            ELE_TYPE l = m1->values[i * BLOCK_SIDE + j];
            //printf("第%d行消%d行\n", j, i);
            for (int k = m2->row_pointers[j]; k < m2->row_pointers[j + 1]; k++) {
                int c = m2->col_indices[k];
                //printf("%lf -= %lf * %lf\n", A(i, c), l, A(j, c));
                A(i, c) -= l * A(j, c);
            }
        }
    }
}

void DGESSM_DD(const ELE_TYPE *a, ELE_TYPE *b) {
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int j = i - 1; j < BLOCK_SIDE; ++j) {
            ELE_TYPE l = DA(j, i);
            for (int c = 0; c < BLOCK_SIDE; ++c) {
                DB(j, c) -= l * DB(i, c);
            }
        }
    }
}

///稀疏L解稠密矩阵
void DGESSM_DS(const ELE_TYPE *lx, BlockMatrix *m1) {
    ELE_TYPE *a = m1->values;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int j = i; j < BLOCK_SIDE; j++) {
            ELE_TYPE l = L(j, i);
            //printf("第%d行消%d行\n", i, j);
            for (int k = m1->row_pointers[i]; k < m1->row_pointers[i + 1]; k++) {
                int c = m1->col_indices[k];
                //printf("%lf -= %lf * %lf\n", A(j, c), l, A(i, c));
                A(j, c) -= l * A(i, c);
            }
        }
    }
}

void DGESSM_SD(const int *Lp_start, const int *Lp_end, BlockMatrix *m1, ELE_TYPE *a) {
    ELE_TYPE *lx = m1->values;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int p = Lp_start[i]; p < Lp_end[i]; p++) {
            int j = m1->row_indices[p];
            ELE_TYPE l = L(j, i);
            //printf("第%d行消%d行\n", i, j);
            for (int k = 0; k < BLOCK_SIDE; k++) {
                //printf("%lf -= %lf * %lf\n", A(j, c), l, A(i, c));
                A(j, k) -= l * A(i, k);
            }
        }
    }
}

int main() {
    int row_ptr[] = {0, 0, 1, 2, 3, 5}; // 每行非零元素的起始位置
    int col_indices[] = {0, 0, 1, 0, 2}; // 非零元素的列索引

    int col_ptr[] = {0, 3, 4, 5, 5, 5}; // 每行非零元素的起始位置
    int row_indices[] = {1, 2, 4, 3, 4}; // 非零元素的列索引

    ELE_TYPE values[] = {
        0, 0, 0, 0, 0,
        1, 0, 0, 0, 0,
        2, 0, 0, 0, 0,
        0, 4, 0, 0, 0,
        5, 0, 6, 0, 0
    }; // 非零元素的值
    BlockMatrix m1 = {
        .l_nnz = 6,
        .col_pointers = col_ptr,
        .row_indices = row_indices,
        // .col_pointers = col_ptr,
        // .row_indices = row_indices,
        .values = values,
    };
    int row_ptr2[] = {0, 2, 5, 9, 12, 16}; // 每行非零元素的起始位置
    int col_indices2[] = {1, 3, 1, 3, 4, 0, 1, 2, 3, 1, 3, 4, 0, 1, 2, 3}; // 非零元素的列索引
    ELE_TYPE values2[] = {
        0, 1, 0, 2, 0,
        0, 3, 0, 0, 4,
        5, 0, 6, 0, 0,
        0, 0, 0, 7, 0,
        0, 8, 0, 0, 0
    }; // 非零元素的值

    BlockMatrix m2 = {
        .u_nnz = 8,
        .row_pointers = row_ptr2,
        .col_indices = col_indices2,
        .values = values2,
    };
    ELE_TYPE *c = (ELE_TYPE *) calloc(BLOCK_SIDE * BLOCK_SIDE, sizeof(ELE_TYPE));
    // DGESSM_SS(m1.col_pointers, m1.col_pointers + 1, &m1, &m2);
    // DGESSM_SS2(m1.row_pointers, m1.row_pointers + 1, &m1, &m2);
    // DGESSM_DD(m1.values, m2.values);
    DGESSM_SD(m1.col_pointers, m1.col_pointers + 1, &m1, m2.values);
    // DGESSM_DS(m1.values, &m2);
    printf("解：\n");
    print_dense_matrix(m2.values,BLOCK_SIDE);

    memcpy(c, values2, BLOCK_SIDE * BLOCK_SIDE * sizeof(ELE_TYPE));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                BLOCK_SIDE, BLOCK_SIDE, BLOCK_SIDE, 1.0,
                m1.values, BLOCK_SIDE, m2.values, BLOCK_SIDE, 1.0, c, BLOCK_SIDE);
    printf("还原\n");
    print_dense_matrix(c,BLOCK_SIDE);
    return EXIT_SUCCESS;
}
