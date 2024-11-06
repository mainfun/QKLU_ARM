#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <layer_matrix.h>
#include <cblas.h>

#define ELE_TYPE double
#define INDEX_TYPE long long
#define BLOCK_SIDE 5
#define A(i,j) a[(i)*BLOCK_SIDE+(j)]
#define B(i,j) b[(i)*BLOCK_SIDE+(j)]
#define C(i,j) c[(i)*BLOCK_SIDE+(j)]

///核心是SpGEMM m3+=m1*m2,m1 m2都是CSR的
void SpGEMM_v2(BlockMatrix *m1, BlockMatrix *m2, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    ELE_TYPE *c = m3->array2D;
    //c+=a*b
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int j = m1->row_pointers[i]; j < m1->row_pointers[i + 1]; j++) {
            const int a_col = m1->col_indices[j];
            const ELE_TYPE a_v = a[j];
            for (int k = m2->row_pointers[a_col]; k < m2->row_pointers[a_col + 1]; k++) {
                const int b_col = m2->col_indices[k];
                const ELE_TYPE b_v = b[k];
                C(i, b_col) += a_v * b_v;
            }
        }
    }
}

///核心是SpGEMM m3+=m1*m2，m1 CSC, m2 CSR
void SpGEMM(BlockMatrix *m1, BlockMatrix *m2, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    ELE_TYPE *c = m3->values;
    int *offset3 = m3->offset;
    //c+=a*b
    if (m3->format == SPARSE) {
        for (int i = 0; i < BLOCK_SIDE; ++i) {
            for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
                const int j = m1->row_indices[p];
                const ELE_TYPE l = A(j, i);
                for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                    const int b_col = m2->col_indices[k];
                    C(j, b_col) -= l * B(i, b_col);
                }
            }
        }
    } else {
        for (int i = 0; i < BLOCK_SIDE; ++i) {
            for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
                const int j = m1->row_indices[p];
                const ELE_TYPE l = A(j, i);
                for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                    const int b_col = m2->col_indices[k];
                    C(j, b_col) -= l * B(i, b_col);
                }
            }
        }
    }
}

void SpMM_DS(const ELE_TYPE *a, const BlockMatrix *m2, BlockMatrix *m3) {
    const ELE_TYPE *b = m2->values;
    ELE_TYPE *c = m3->values;
    int *offset3 = m3->offset;
    //c+=a*b
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int j = 0; j < BLOCK_SIDE; j++) {
            const ELE_TYPE l = A(j, i);
            for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                const int b_col = m2->col_indices[k];
                C(j, b_col) -= l * B(i, b_col);
            }
        }
    }
}

void SpMM_SD(const BlockMatrix *m1, const ELE_TYPE *b, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *c = m3->values;
    int *offset3 = m3->offset;
    //c+=a*b
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
            const int j = m1->row_indices[p];
            const ELE_TYPE l = A(j, i);
            for (int b_col = 0; b_col < BLOCK_SIDE; b_col++) {
                C(j, b_col) -= l * B(i, b_col);
            }
        }
    }
}

//C=α(A×B)+βC
void GEMM(const ELE_TYPE *a, const ELE_TYPE *b, ELE_TYPE *c) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                BLOCK_SIDE, BLOCK_SIDE, BLOCK_SIDE, -1.0, a, BLOCK_SIDE, b, BLOCK_SIDE, 1.0, c, BLOCK_SIDE);
}


int main() {
    printf("BLAS version: %s\n", OPENBLAS_VERSION);
    int row_ptr[] = {0, 3, 5, 8, 9, 11}; // 每行非零元素的起始位置
    int col_indices[] = {0, 2, 4, 1, 3, 0, 2, 3, 1, 0, 2}; // 非零元素的列索引

    int col_ptr[] = {0, 3, 5, 8, 10, 11};
    int row_idx[] = {
        0, 2, 4,
        1, 3,
        0, 2, 4,
        1, 2,
        0
    };

    ELE_TYPE arr[] = {
        1, 0, 2, 0, 3, //0
        0, 4, 0, 5, 0, //1
        6, 0, 7, 8, 0, //2
        0, 9, 0, 0, 0, //3
        10, 0, 11, 0, 0 //4
        //0,3, 5, 8, 10,11
    }; // 非零元素的值
    //ELE_TYPE values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    BlockMatrix m1 = {
        .l_nnz = 11,
        .values = arr,
        .col_pointers = col_ptr,
        .row_indices = row_idx,
    };

    int row_ptr2[] = {0, 2, 4, 6, 7, 8}; // 每行非零元素的起始位置
    int col_indices2[] = {1, 3, 1, 4, 0, 2, 3, 1}; // 非零元素的列索引
    ELE_TYPE arr2[] = {
        0, 1, 0, 2, 0,
        0, 3, 0, 0, 4,
        5, 0, 6, 0, 0,
        0, 0, 0, 7, 0,
        0, 8, 0, 0, 0
    }; // 非零元素的值
    ELE_TYPE values2[] = {1, 2, 3, 4, 5, 6, 7, 8};
    BlockMatrix m2 = {
        .u_nnz = 8,
        .values = arr2,
        //.array2D = arr2,
        .row_pointers = row_ptr2,
        .col_indices = col_indices2,
    };
    ELE_TYPE *c = (ELE_TYPE *) calloc(BLOCK_SIDE * BLOCK_SIDE, sizeof(ELE_TYPE));
    // SpGEMM_v2(&m1, &m2, c);
    // GEMM(m1.values, m2.values, c);
    //SpMM_DS(m1.values, &m2, c);
    SpMM_SD(&m1, m2.values, c);

    printf("m3:\n");
    print_dense_matrix(c,BLOCK_SIDE);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                BLOCK_SIDE, BLOCK_SIDE, BLOCK_SIDE, 1.0, arr, BLOCK_SIDE, arr2, BLOCK_SIDE, 1.0, c, BLOCK_SIDE);

    printf("误差:\n");
    print_dense_matrix(c,BLOCK_SIDE);
    return 0;
}
