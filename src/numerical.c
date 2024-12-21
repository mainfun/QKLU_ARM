#include <numerical.h>
#include <preprocess.h>
#include <cblas.h>
#include <omp.h>

long long sp_fma_count = 0;
long long dense_fma_count = 0;
long long other_fma_count = 0;
double gemm_time = 0;
double spgemm_time = 0;

void SpGEMM_v0(BlockMatrix *m1, BlockMatrix *m2, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    ELE_TYPE *c = m3->values;
    int *offset1 = m1->offset;
    int *offset2 = m2->offset;
    int *offset3 = m3->offset;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
            const int j = m1->row_indices[p];
            const ELE_TYPE l = A(j, i);
            ELE_TYPE *b_ptr = b + offset2[i];
            ELE_TYPE *c_ptr = c + offset3[j];
            for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                const int b_col = m2->col_indices[k];
                c_ptr[b_col] -= l * b_ptr[b_col];
                sp_fma_count++;
            }
        }
    }
}

long long count_dense_row = 0;
long long count_sparse_row = 0;

void SpGEMM_v3(BlockMatrix *m1, BlockMatrix *m2, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    ELE_TYPE *c = m3->values;
    int *offset1 = m1->offset;
    int *offset2 = m2->offset;
    int *offset3 = m3->offset;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        if (m2->row_pointers[i] == m2->row_pointers[i + 1]) continue;
        int start_col = m2->col_indices[m2->row_pointers[i]];
        int end_col = m2->col_indices[m2->row_pointers[i + 1] - 1];
        if (end_col - start_col < 1.5 * (m2->row_pointers[i + 1] - m2->row_pointers[i])) {
            for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
                count_dense_row += m2->row_pointers[i + 1] - m2->row_pointers[i];
                const int j = m1->row_indices[p];
                const ELE_TYPE l = A(j, i);
                ELE_TYPE *b_ptr = b + offset2[i];
                ELE_TYPE *c_ptr = c + offset3[j];
                for (int k = start_col; k <= end_col; k++) c_ptr[k] -= l * b_ptr[k];
            }
        } else {
            for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
                count_sparse_row += m2->row_pointers[i + 1] - m2->row_pointers[i];
                const int j = m1->row_indices[p];
                const ELE_TYPE l = A(j, i);
                ELE_TYPE *b_ptr = b + offset2[i];
                ELE_TYPE *c_ptr = c + offset3[j];
                for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                    const int b_col = m2->col_indices[k];
                    c_ptr[b_col] -= l * b_ptr[b_col];
                }
            }
        }
    }
}

void SpGEMM_v4(BlockMatrix *m1, BlockMatrix *m2, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    ELE_TYPE *c = m3->values;
    int *offset1 = m1->offset;
    int *offset2 = m2->offset;
    int *offset3 = m3->offset;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        if (m2->row_pointers[i] == m2->row_pointers[i + 1]) continue;
        int start_col = m2->col_indices[m2->row_pointers[i]];
        int end_col = m2->col_indices[m2->row_pointers[i + 1] - 1];
        if (end_col - start_col < 1.5 * (m2->row_pointers[i + 1] - m2->row_pointers[i])) {
            for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
                count_dense_row += m2->row_pointers[i + 1] - m2->row_pointers[i];
                const int j = m1->row_indices[p];
                const ELE_TYPE l = A(j, i);
                ELE_TYPE *b_ptr = b + offset2[i];
                ELE_TYPE *c_ptr = c + offset3[j];
                for (int k = start_col; k <= end_col; k++) c_ptr[k] -= l * b_ptr[k];
            }
        } else {
            int p = m1->col_pointers[i];
            for (; p <= m1->col_pointers[i + 1] - 6; p += 6) {
                const int j0 = m1->row_indices[p + 0];
                const int j1 = m1->row_indices[p + 1];
                const int j2 = m1->row_indices[p + 2];
                const int j3 = m1->row_indices[p + 3];
                const int j4 = m1->row_indices[p + 4];
                const int j5 = m1->row_indices[p + 5];
                const ELE_TYPE l0 = A(j0, i);
                const ELE_TYPE l1 = A(j1, i);
                const ELE_TYPE l2 = A(j2, i);
                const ELE_TYPE l3 = A(j3, i);
                const ELE_TYPE l4 = A(j4, i);
                const ELE_TYPE l5 = A(j5, i);
                const ELE_TYPE *b_ptr = b + offset2[i];
                ELE_TYPE *c_ptr0 = c + offset3[j0];
                ELE_TYPE *c_ptr1 = c + offset3[j1];
                ELE_TYPE *c_ptr2 = c + offset3[j2];
                ELE_TYPE *c_ptr3 = c + offset3[j3];
                ELE_TYPE *c_ptr4 = c + offset3[j4];
                ELE_TYPE *c_ptr5 = c + offset3[j5];
                for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                    const int b_col = m2->col_indices[k];
                    const ELE_TYPE pv = b_ptr[b_col];
                    c_ptr0[b_col] -= l0 * pv;
                    c_ptr1[b_col] -= l1 * pv;
                    c_ptr2[b_col] -= l2 * pv;
                    c_ptr3[b_col] -= l3 * pv;
                    c_ptr4[b_col] -= l4 * pv;
                    c_ptr5[b_col] -= l5 * pv;
                }
            }
            for (; p <= m1->col_pointers[i + 1] - 2; p += 2) {
                const int j = m1->row_indices[p];
                const int j1 = m1->row_indices[p + 1];
                const ELE_TYPE l = A(j, i);
                const ELE_TYPE l1 = A(j1, i);
                ELE_TYPE *b_ptr = b + offset2[i];
                ELE_TYPE *c_ptr = c + offset3[j];
                ELE_TYPE *c_ptr1 = c + offset3[j1];
                for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                    const int b_col = m2->col_indices[k];
                    const ELE_TYPE pv = b_ptr[b_col];
                    c_ptr[b_col] -= l * pv;
                    c_ptr1[b_col] -= l1 * pv;
                }
            }
            for (; p < m1->col_pointers[i + 1]; p++) {
                const int j = m1->row_indices[p];
                const ELE_TYPE l = A(j, i);
                ELE_TYPE *b_ptr = b + offset2[i];
                ELE_TYPE *c_ptr = c + offset3[j];
                for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                    const int b_col = m2->col_indices[k];
                    c_ptr[b_col] -= l * b_ptr[b_col];
                }
            }
        }
    }
}


void SpGEMM_v1(BlockMatrix *m1, BlockMatrix *m2, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    ELE_TYPE *c = m3->values;
    int *offset1 = m1->offset;
    int *offset2 = m2->offset;
    int *offset3 = m3->offset;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        if (m2->row_pointers[i] == m2->row_pointers[i + 1]) continue;
        int p = m1->col_pointers[i];
        for (; p <= m1->col_pointers[i + 1] - 6; p += 6) {
            const int j0 = m1->row_indices[p + 0];
            const int j1 = m1->row_indices[p + 1];
            const int j2 = m1->row_indices[p + 2];
            const int j3 = m1->row_indices[p + 3];
            const int j4 = m1->row_indices[p + 4];
            const int j5 = m1->row_indices[p + 5];
            const ELE_TYPE l0 = A(j0, i);
            const ELE_TYPE l1 = A(j1, i);
            const ELE_TYPE l2 = A(j2, i);
            const ELE_TYPE l3 = A(j3, i);
            const ELE_TYPE l4 = A(j4, i);
            const ELE_TYPE l5 = A(j5, i);
            const ELE_TYPE *b_ptr = b + offset2[i];
            ELE_TYPE *c_ptr0 = c + offset3[j0];
            ELE_TYPE *c_ptr1 = c + offset3[j1];
            ELE_TYPE *c_ptr2 = c + offset3[j2];
            ELE_TYPE *c_ptr3 = c + offset3[j3];
            ELE_TYPE *c_ptr4 = c + offset3[j4];
            ELE_TYPE *c_ptr5 = c + offset3[j5];
            for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                const int b_col = m2->col_indices[k];
                const ELE_TYPE pv = b_ptr[b_col];
                c_ptr0[b_col] -= l0 * pv;
                c_ptr1[b_col] -= l1 * pv;
                c_ptr2[b_col] -= l2 * pv;
                c_ptr3[b_col] -= l3 * pv;
                c_ptr4[b_col] -= l4 * pv;
                c_ptr5[b_col] -= l5 * pv;
                sp_fma_count += 6;
            }
        }
        for (; p <= m1->col_pointers[i + 1] - 2; p += 2) {
            const int j = m1->row_indices[p];
            const int j1 = m1->row_indices[p + 1];
            const ELE_TYPE l = A(j, i);
            const ELE_TYPE l1 = A(j1, i);
            ELE_TYPE *b_ptr = b + offset2[i];
            ELE_TYPE *c_ptr = c + offset3[j];
            ELE_TYPE *c_ptr1 = c + offset3[j1];
            for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                const int b_col = m2->col_indices[k];
                const ELE_TYPE pv = b_ptr[b_col];
                c_ptr[b_col] -= l * pv;
                c_ptr1[b_col] -= l1 * pv;
                sp_fma_count += 2;
            }
        }
        for (; p < m1->col_pointers[i + 1]; p++) {
            const int j = m1->row_indices[p];
            const ELE_TYPE l = A(j, i);
            ELE_TYPE *b_ptr = b + offset2[i];
            ELE_TYPE *c_ptr = c + offset3[j];
            for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                const int b_col = m2->col_indices[k];
                c_ptr[b_col] -= l * b_ptr[b_col];
                sp_fma_count++;
            }
        }
    }
}

void LUS(const int *Lp_start, const int *Lp_end, const int *Li,
         const int *Up_start, const int *Up_end, const int *Ui,
         ELE_TYPE *a, int n, const int *offset1) {
    for (int i = 0; i < n; ++i) {
        for (int p = Lp_start[i]; p < Lp_end[i]; p++) {
            int j = Li[p];
            //printf("第%d行消第%d行\n", i, j);
            ELE_TYPE l = A(j, i) /= A(i, i);
            ELE_TYPE *pivot_row_ptr = a + offset1[i];
            ELE_TYPE *eli_row_ptr = a + offset1[j];
            for (int k = Up_start[i] + 1; k < Up_end[i]; k++) {
                int c = Ui[k];
                //printf("eli_row_ptr[%d]-= %lf * pivot_row_ptr[%d]\n", c, l, c);
                eli_row_ptr[c] -= l * pivot_row_ptr[c];
                other_fma_count++;
            }
        }
    }
}

void LUS_v1(BlockMatrix *m1) {
    int *offset1 = m1->offset;
    ELE_TYPE *a = m1->values;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        int p = m1->col_pointers[i];
        for (; p <= m1->col_pointers[i + 1] - 2; p += 2) {
            int j0 = m1->row_indices[p + 0];
            int j1 = m1->row_indices[p + 1];
            // int j2 = m1->row_indices[p + 2];
            // int j3 = m1->row_indices[p + 3];
            // int j4 = m1->row_indices[p + 4];
            // int j5 = m1->row_indices[p + 5];
            ELE_TYPE l0 = A(j0, i) /= A(i, i);
            ELE_TYPE l1 = A(j1, i) /= A(i, i);
            // ELE_TYPE l2 = A(j2, i) /= A(i, i);
            // ELE_TYPE l3 = A(j3, i) /= A(i, i);
            // ELE_TYPE l4 = A(j4, i) /= A(i, i);
            // ELE_TYPE l5 = A(j5, i) /= A(i, i);
            ELE_TYPE *pivot_row_ptr = a + offset1[i];
            ELE_TYPE *e0 = a + offset1[j0];
            ELE_TYPE *e1 = a + offset1[j1];
            // ELE_TYPE *e2 = a + offset1[j2];
            // ELE_TYPE *e3 = a + offset1[j3];
            // ELE_TYPE *e4 = a + offset1[j4];
            // ELE_TYPE *e5 = a + offset1[j5];
            for (int k = m1->row_pointers[i] + 1; k < m1->row_pointers[i + 1]; k++) {
                int c = m1->col_indices[k];
                const ELE_TYPE pv = pivot_row_ptr[c];
                e0[c] -= l0 * pv;
                e1[c] -= l1 * pv;
                // e2[c] -= l2 * pv;
                // e3[c] -= l3 * pv;
                // e4[c] -= l4 * pv;
                // e5[c] -= l5 * pv;
            }
        }
        //余量
        for (; p < m1->col_pointers[i + 1]; p++) {
            int j = m1->row_indices[p];
            ELE_TYPE l = A(j, i) /= A(i, i);
            ELE_TYPE *pivot_row_ptr = a + offset1[i];
            ELE_TYPE *eli_row_ptr = a + offset1[j];
            for (int k = m1->row_pointers[i] + 1; k < m1->row_pointers[i + 1]; k++) {
                int c = m1->col_indices[k];
                eli_row_ptr[c] -= l * pivot_row_ptr[c];
            }
        }
    }
}

void SpGEMM_v2(BlockMatrix *m1, BlockMatrix *m2, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    ELE_TYPE *c = m3->values;
    int *offset1 = m1->offset;
    int *offset2 = m2->offset;
    int *offset3 = m3->offset;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        int p = m1->col_pointers[i];
        for (; p <= m1->col_pointers[i + 1] - 6; p += 6) {
            const int j0 = m1->row_indices[p + 0];
            const int j1 = m1->row_indices[p + 1];
            const int j2 = m1->row_indices[p + 2];
            const int j3 = m1->row_indices[p + 3];
            const int j4 = m1->row_indices[p + 4];
            const int j5 = m1->row_indices[p + 5];
            const ELE_TYPE l0 = A(j0, i);
            const ELE_TYPE l1 = A(j1, i);
            const ELE_TYPE l2 = A(j2, i);
            const ELE_TYPE l3 = A(j3, i);
            const ELE_TYPE l4 = A(j4, i);
            const ELE_TYPE l5 = A(j5, i);
            const ELE_TYPE *b_ptr = b + offset2[i];
            ELE_TYPE *c_ptr0 = c + offset3[j0];
            ELE_TYPE *c_ptr1 = c + offset3[j1];
            ELE_TYPE *c_ptr2 = c + offset3[j2];
            ELE_TYPE *c_ptr3 = c + offset3[j3];
            ELE_TYPE *c_ptr4 = c + offset3[j4];
            ELE_TYPE *c_ptr5 = c + offset3[j5];
            for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                const int b_col = m2->col_indices[k];
                const ELE_TYPE pv = b_ptr[b_col];
                c_ptr0[b_col] -= l0 * pv;
                c_ptr1[b_col] -= l1 * pv;
                c_ptr2[b_col] -= l2 * pv;
                c_ptr3[b_col] -= l3 * pv;
                c_ptr4[b_col] -= l4 * pv;
                c_ptr5[b_col] -= l5 * pv;
            }
        }
        for (; p <= m1->col_pointers[i + 1] - 2; p += 2) {
            const int j = m1->row_indices[p];
            const int j1 = m1->row_indices[p + 1];
            const ELE_TYPE l = A(j, i);
            const ELE_TYPE l1 = A(j1, i);
            ELE_TYPE *b_ptr = b + offset2[i];
            ELE_TYPE *c_ptr = c + offset3[j];
            ELE_TYPE *c_ptr1 = c + offset3[j1];
            for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                const int b_col = m2->col_indices[k];
                const ELE_TYPE pv = b_ptr[b_col];
                c_ptr[b_col] -= l * pv;
                c_ptr1[b_col] -= l1 * pv;
            }
        }
        for (; p < m1->col_pointers[i + 1]; p++) {
            const int j = m1->row_indices[p];
            const ELE_TYPE l = A(j, i);
            ELE_TYPE *b_ptr = b + offset2[i];
            ELE_TYPE *c_ptr = c + offset3[j];
            for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                const int b_col = m2->col_indices[k];
                c_ptr[b_col] -= l * b_ptr[b_col];
            }
        }
    }
}

void DTSTRF_SS_v1(BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    int *offset1 = m1->offset;
    int *offset2 = m2->offset;
    for (int i = 0; i < BLOCK_SIDE; i++) {
        int j = m2->col_pointers[i];
        for (; j <= m2->col_pointers[i + 1] - 4; j += 4) {
            const int b_row0 = m2->row_indices[j + 0];
            const int b_row1 = m2->row_indices[j + 1];
            const int b_row2 = m2->row_indices[j + 2];
            const int b_row3 = m2->row_indices[j + 3];
            const ELE_TYPE l0 = B(b_row0, i) /= A(i, i);
            const ELE_TYPE l1 = B(b_row1, i) /= A(i, i);
            const ELE_TYPE l2 = B(b_row2, i) /= A(i, i);
            const ELE_TYPE l3 = B(b_row3, i) /= A(i, i);
            const ELE_TYPE *a_ptr = a + offset1[i];
            ELE_TYPE *b_ptr0 = b + offset2[b_row0];
            ELE_TYPE *b_ptr1 = b + offset2[b_row1];
            ELE_TYPE *b_ptr2 = b + offset2[b_row2];
            ELE_TYPE *b_ptr3 = b + offset2[b_row3];
            for (int k = m1->row_pointers[i] + 1; k < m1->row_pointers[i + 1]; k++) {
                const int a_col = m1->col_indices[k];
                const ELE_TYPE pv = a_ptr[a_col];
                b_ptr0[a_col] -= l0 * pv;
                b_ptr1[a_col] -= l1 * pv;
                b_ptr2[a_col] -= l2 * pv;
                b_ptr3[a_col] -= l3 * pv;
            }
        }
        for (; j < m2->col_pointers[i + 1]; ++j) {
            int b_row = m2->row_indices[j];
            ELE_TYPE l = B(b_row, i) /= A(i, i);
            for (int k = m1->row_pointers[i] + 1; k < m1->row_pointers[i + 1]; k++) {
                int a_col = m1->col_indices[k];
                B(b_row, a_col) -= l * A(i, a_col);
            }
        }
    }
}

void LUD(ELE_TYPE *a) {
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        ELE_TYPE pivot = DA(i, i);
        for (int j = i + 1; j < BLOCK_SIDE; ++j) {
            ELE_TYPE l = DA(j, i) /= pivot;
            if (l == 0) continue;
            for (int k = i + 1; k < BLOCK_SIDE; ++k) {
                DA(j, k) -= l * DA(i, k);
                other_fma_count++;
            }
        }
    }
}

//U CSR L CSC
void DTSTRF_SS(BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    int *offset1 = m1->offset;
    int *offset2 = m2->offset;
    for (int i = 0; i < BLOCK_SIDE; i++) {
        for (int j = m2->col_pointers[i]; j < m2->col_pointers[i + 1]; ++j) {
            const int b_row = m2->row_indices[j];
            const ELE_TYPE l = B(b_row, i) /= A(i, i);
            for (int k = m1->row_pointers[i] + 1; k < m1->row_pointers[i + 1]; k++) {
                const int a_col = m1->col_indices[k];
                B(b_row, a_col) -= l * A(i, a_col);
                other_fma_count++;
            }
        }
    }
}

// void DTSTRF_DD(const ELE_TYPE *a, ELE_TYPE *b) {
//     for (int i = 0; i < BLOCK_SIDE; i++) {
//         for (int j = 0; j < BLOCK_SIDE; j++) {
//             //第j行消第i行
//             ELE_TYPE l = DB(i, j) /= DA(j, j);
//             if (l == 0) continue;
//             for (int k = j + 1; k < BLOCK_SIDE; k++) {
//                 DB(i, k) -= l * DA(j, k);
//             }
//         }
//     }
// }
//
// //U稀疏
// void DTSTRF_SD(const BlockMatrix *m1, ELE_TYPE *b) {
//     ELE_TYPE *a = m1->values;
//     int *offset1 = m1->offset;
//     for (int i = 0; i < BLOCK_SIDE; i++) {
//         for (int j = 0; j < BLOCK_SIDE; ++j) {
//             DB(j, i) /= A(i, i);
//         }
//         for (int k = m1->row_pointers[i] + 1; k < m1->row_pointers[i + 1]; k++) {
//             int a_col = m1->col_indices[k];
//             ELE_TYPE a_i_a_col = A(i, a_col);
//             for (int j = 0; j < BLOCK_SIDE; ++j) {
//                 DB(j, a_col) -= DB(j, i) * a_i_a_col;
//             }
//         }
//     }
// }

// //U稠密
// void DTSTRF_DS(const ELE_TYPE *a, BlockMatrix *m2) {
//     ELE_TYPE *b = m2->values;
//     int *offset2 = m2->offset;
//     for (int i = 0; i < BLOCK_SIDE; i++) {
//         for (int j = m2->col_pointers[i]; j < m2->col_pointers[i + 1]; ++j) {
//             int b_row = m2->row_indices[j];
//             ELE_TYPE l = B(b_row, i) /= DA(i, i);
//             for (int k = i + 1; k < BLOCK_SIDE; k++) {
//                 B(b_row, k) -= l * DA(i, k);
//             }
//         }
//     }
// }

///用m1解m2 L是CSC的 U csr
void DGESSM_SS(BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *lx = m1->values;
    int *offset_l = m1->offset;
    ELE_TYPE *a = m2->values;
    int *offset1 = m2->offset;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
            int j = m1->row_indices[p];
            ELE_TYPE l = L(j, i);
            //printf("第%d行消%d行\n", i, j);
            for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                int c = m2->col_indices[k];
                //printf("%lf -= %lf * %lf\n", A(j, c), l, A(i, c));
                A(j, c) -= l * A(i, c);
                other_fma_count++;
            }
        }
    }
}

// void DGESSM_DD(const ELE_TYPE *a, ELE_TYPE *b) {
//     for (int i = 0; i < BLOCK_SIDE; ++i) {
//         for (int j = i - 1; j < BLOCK_SIDE; ++j) {
//             ELE_TYPE l = DA(j, i);
//             if (l == 0) continue;
//             for (int c = 0; c < BLOCK_SIDE; ++c) {
//                 DB(j, c) -= l * DB(i, c);
//             }
//         }
//     }
// }
//
// ///稀疏L解稠密矩阵
// void DGESSM_DS(const ELE_TYPE *lx, BlockMatrix *m1) {
//     ELE_TYPE *a = m1->values;
//     int *offset1 = m1->offset;
//     for (int i = 0; i < BLOCK_SIDE; ++i) {
//         for (int j = i; j < BLOCK_SIDE; j++) {
//             ELE_TYPE l = DL(j, i);
//             if (l == 0) continue;
//             //printf("第%d行消%d行\n", i, j);
//             for (int k = m1->row_pointers[i]; k < m1->row_pointers[i + 1]; k++) {
//                 int c = m1->col_indices[k];
//                 //printf("%lf -= %lf * %lf\n", A(j, c), l, A(i, c));
//                 A(j, c) -= l * A(i, c);
//             }
//         }
//     }
// }

// void DGESSM_SD(BlockMatrix *m1, ELE_TYPE *a) {
//     ELE_TYPE *lx = m1->values;
//     int *offset_l = m1->offset;
//     for (int i = 0; i < BLOCK_SIDE; ++i) {
//         for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
//             int j = m1->row_indices[p];
//             ELE_TYPE l = L(j, i);
//             //printf("第%d行消%d行\n", i, j);
//             for (int k = 0; k < BLOCK_SIDE; k++) {
//                 //printf("%lf -= %lf * %lf\n", A(j, c), l, A(i, c));
//                 DA(j, k) -= l * DA(i, k);
//             }
//         }
//     }
// }

///核心是SpGEMM m3+=m1*m2，m1 CSC, m2 CSR
void SpGEMM(BlockMatrix *m1, BlockMatrix *m2, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    ELE_TYPE *c = m3->values;
    int *offset1 = m1->offset;
    int *offset2 = m2->offset;
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
                    DC(j, b_col) -= l * B(i, b_col);
                }
            }
        }
    }
}

void SpMM_DS(const ELE_TYPE *a, const BlockMatrix *m2, BlockMatrix *m3) {
    const ELE_TYPE *b = m2->values;
    ELE_TYPE *c = m3->values;
    int *offset2 = m2->offset;
    int *offset3 = m3->offset;
    //c+=a*b
    if (m3->format == SPARSE) {
        for (int i = 0; i < BLOCK_SIDE; ++i) {
            for (int j = 0; j < BLOCK_SIDE; j++) {
                const ELE_TYPE l = DA(j, i);
                if (l == 0)continue;
                for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                    const int b_col = m2->col_indices[k];
                    C(j, b_col) -= l * B(i, b_col);
                }
            }
        }
    } else {
        for (int i = 0; i < BLOCK_SIDE; ++i) {
            for (int j = 0; j < BLOCK_SIDE; j++) {
                const ELE_TYPE l = DA(j, i);
                for (int k = m2->row_pointers[i]; k < m2->row_pointers[i + 1]; k++) {
                    const int b_col = m2->col_indices[k];
                    DC(j, b_col) -= l * B(i, b_col);
                }
            }
        }
    }
}

void SpMM_SD_v1(const BlockMatrix *m1, ELE_TYPE *b, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *c = m3->values;
    int *offset1 = m1->offset;
    int *offset3 = m3->offset;
    //c+=a*b
    if (m3->format == SPARSE) {
        for (int i = 0; i < BLOCK_SIDE; ++i) {
            for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
                const int j = m1->row_indices[p];
                const ELE_TYPE l = A(j, i);
                ELE_TYPE *c_j_ptr = c + offset3[j];
                ELE_TYPE *b_i_ptr = b + i * BLOCK_SIDE;
                for (int b_col = 0; b_col < BLOCK_SIDE; b_col++) {
                    c_j_ptr[b_col] -= l * b_i_ptr[b_col];
                }
            }
        }
    }
}

void SpMM_SD(const BlockMatrix *m1, const ELE_TYPE *b, BlockMatrix *m3) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *c = m3->values;
    int *offset1 = m1->offset;
    int *offset3 = m3->offset;
    //c+=a*b
    if (m3->format == SPARSE) {
        for (int i = 0; i < BLOCK_SIDE; ++i) {
            for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
                const int j = m1->row_indices[p];
                const ELE_TYPE l = A(j, i);
                for (int b_col = 0; b_col < BLOCK_SIDE; b_col++) {
                    C(j, b_col) -= l * DB(i, b_col);
                }
            }
        }
    } else {
        for (int i = 0; i < BLOCK_SIDE; ++i) {
            for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
                const int j = m1->row_indices[p];
                const ELE_TYPE l = A(j, i);
                for (int b_col = 0; b_col < BLOCK_SIDE; b_col++) {
                    DC(j, b_col) -= l * DB(i, b_col);
                }
            }
        }
    }
}

//C=α(A×B)+βC
void GEMM(const ELE_TYPE *a, const ELE_TYPE *b, ELE_TYPE *c) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                BLOCK_SIDE, BLOCK_SIDE, BLOCK_SIDE, -1.0, a, BLOCK_SIDE, b, BLOCK_SIDE, 1.0, c, BLOCK_SIDE);
}


void sample_factor(BlockMatrix *m) {
    if (m->format == DENSE) {
        LUD(m->values);
    } else {
        LUS(m->col_pointers, m->col_pointers + 1, m->row_indices,
            m->row_pointers, m->row_pointers + 1, m->col_indices,
            m->values,BLOCK_SIDE, m->offset);
        // LUS_v1(m);
    }
}


// void upper_solving(BlockMatrix *m1, BlockMatrix *m2) {
//     // if (m1->format == DENSE) {
//     //     if (m2->format == DENSE) {
//     //         DTSTRF_DD(m1->values, m2->values);
//     //     } else {
//     //         DTSTRF_DS(m1->values, m2);
//     //     }
//     // } else {
//     //     if (m2->format == DENSE) {
//     //         DTSTRF_SD(m1, m2->values);
//     //     } else {
//     //         DTSTRF_SS(m1, m2);
//     //     }
//     // }
// }

// void lower_solving(BlockMatrix *m1, BlockMatrix *m2) {
//     if (m1->format == DENSE) {
//         if (m2->format == DENSE) {
//             DGESSM_DD(m1->values, m2->values);
//         } else {
//             DGESSM_DS(m1->values, m2);
//         }
//     } else {
//         if (m2->format == DENSE) {
//             DGESSM_SD(m1, m2->values);
//         } else {
//             DGESSM_SS(m1, m2);
//         }
//     }
// }

void add_matrix_DS_csc(BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    int *offset1 = m1->offset;
    int *offset2 = m2->offset;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
            const int r = m1->row_indices[p];
            b[offset2[r] + i] -= a[offset1[r] + i];
        }
    }
}

void add_matrix_DS_csr(BlockMatrix *m1, BlockMatrix *m2) {
    ELE_TYPE *a = m1->values;
    ELE_TYPE *b = m2->values;
    int *offset1 = m1->offset;
    int *offset2 = m2->offset;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int p = m1->row_pointers[i]; p < m1->row_pointers[i + 1]; p++) {
            const int c = m1->col_indices[p];
            b[offset2[i] + c] -= a[offset1[i] + c];
        }
    }
}


///m3+=m1*m2
void schur_complement(BlockMatrix *m1, BlockMatrix *m2, BlockMatrix *m3) {
    // if (m1->format == DENSE) {
    //     if (m2->format == DENSE && m3->format == DENSE) {
    //         GEMM(m1->values, m2->values, m3->values);
    //     } else {
    //         SpMM_DS(m1->values, m2, m3);
    //     }
    // } else {
    //     if (m2->format == DENSE) {
    //         SpMM_SD(m1, m2->values, m3);
    //     } else {
    //         SpGEMM_v3(m1, m2, m3);
    //     }
    // }
    if (m1->format == DENSE && m2->format == DENSE && m3->format == DENSE) {
        double t = omp_get_wtime();
        GEMM(m1->values, m2->values, m3->values);
        gemm_time += omp_get_wtime() - t;
        dense_fma_count += BLOCK_SIDE * BLOCK_SIDE * BLOCK_SIDE;
    } else {
        double t = omp_get_wtime();
        SpGEMM_v1(m1, m2, m3);
        spgemm_time += omp_get_wtime() - t;
    }
}

void block_factor_up_looking_v1(L2Matrix *l2) {
    int u_nnz_0_500 = 0;
    int u_nnz_500_2000 = 0;
    int u_nnz_2000_4000 = 0;
    int u_nnz_4000 = 0;
    double mm_time = 0;
    for (int i = 00; i < l2->num_row_block; ++i) {
        for (int p = l2->row_pointers[i]; p < l2->diag_index[i]; p++) {
            int j = l2->col_indices[p];
            //上三角解
            BlockMatrix *l_bm = get_block(l2, i, j);
            BlockMatrix *diag_block = get_diag_block(l2, j);
            DTSTRF_SS(diag_block, l_bm);
            //schur a(i,c)-=l(i,j) * u(j,c)
            for (int k = l2->diag_index[j] + 1; k < l2->row_pointers[j + 1]; k++) {
                int col_idx = l2->col_indices[k];
                BlockMatrix *bm = get_block(l2, i, col_idx);
                if (bm != NULL) {
                    double t = omp_get_wtime();
                    BlockMatrix *u_bm = get_block(l2, j, col_idx);
                    schur_complement(l_bm, u_bm, bm);
                    mm_time += omp_get_wtime() - t;
                }
            }
        }
        //LU factor
        BlockMatrix *diag_block = get_diag_block(l2, i);
        sample_factor(diag_block);
        //下三角解
        for (int p = l2->diag_index[i] + 1; p < l2->row_pointers[i + 1]; p++) {
            const int c = l2->col_indices[p];
            DGESSM_SS(diag_block, get_block(l2, i, c));
        }
    }
    printf("\nu_nnz: %d,%d,%d,%d\n\n", u_nnz_0_500, u_nnz_500_2000, u_nnz_2000_4000, u_nnz_4000);
    LOG_INFO("DGESSM elapsed time: %lf ms", (mm_time) * 1000.0);
}

///块间右看法
void block_factor_right_looking_v0(L2Matrix *l2) {
    for (int i = 0; i < l2->num_row_block; ++i) {
        //LU分解 A[i][i]
        BlockMatrix *diag_block = get_diag_block(l2, i);
        sample_factor(diag_block);
        for (int j = l2->diag_index[i] + 1; j < l2->row_pointers[i + 1]; ++j) {
            //L[i][i]解A[i][]
            int c = l2->col_indices[j];
            DGESSM_SS(diag_block, get_block(l2, i, c));
        }
        //L CSC
        for (int j = l2->diag_index_csc[i] + 1; j < l2->col_pointers[i + 1]; ++j) {
            int r = l2->row_indices[j];
            BlockMatrix *bm_l = get_block(l2, r, i);
            //U[i][i]解A[][i]
            DTSTRF_SS(diag_block, bm_l);
            for (int k = l2->diag_index[i] + 1; k < l2->row_pointers[i + 1]; ++k) {
                int c = l2->col_indices[k];
                BlockMatrix *bm_A = get_block(l2, r, c);
                if (bm_A != NULL) schur_complement(bm_l, get_block(l2, i, c), bm_A);
            }
        }
    }
}

void block_factor_right_looking_omp_v0(L2Matrix *l2) {
    #pragma omp parallel
    #pragma omp single nowait
    for (int i = 0; i < l2->num_row_block; ++i) {
        //LU分解 A[i][i]
        BlockMatrix *diag_block = get_diag_block(l2, i);
        #pragma omp task depend(inout:*diag_block)
        sample_factor(diag_block);
        for (int j = l2->diag_index[i] + 1; j < l2->row_pointers[i + 1]; ++j) {
            //L[i][i]解A[i][]
            int c = l2->col_indices[j];
            BlockMatrix *u_bm=get_block(l2, i, c);
            #pragma omp task depend(in:*diag_block) depend(inout:*u_bm)
            DGESSM_SS(diag_block, u_bm);
        }
        //L CSC
        for (int j = l2->diag_index_csc[i] + 1; j < l2->col_pointers[i + 1]; ++j) {
            int r = l2->row_indices[j];
            BlockMatrix *bm_l = get_block(l2, r, i);
            //U[i][i]解A[][i]
            #pragma omp task depend(in:*diag_block) depend(inout:*bm_l)
            DTSTRF_SS(diag_block, bm_l);
            for (int k = l2->diag_index[i] + 1; k < l2->row_pointers[i + 1]; ++k) {
                int c = l2->col_indices[k];
                BlockMatrix *bm_A = get_block(l2, r, c);
                BlockMatrix *u_bm=get_block(l2, i, c);
                if (bm_A != NULL) {
                    #pragma omp task depend(in:*bm_l) depend(in:*u_bm) depend(inout:*bm_A)
                    schur_complement(bm_l, u_bm, bm_A);
                }
            }
        }
    }
}

void block_factor(L2Matrix *l2) {
    LOG_INFO("分解开始");
    openblas_set_num_threads(1);
    double factor_time = omp_get_wtime();
    block_factor_right_looking_v0(l2);
    LOG_INFO("分解 elapsed time: %lf ms", (omp_get_wtime() - factor_time) * 1000.0);
    LOG_INFO("count_dense_row=%lld, count_sparse_row=%lld", count_dense_row, count_sparse_row);
    LOG_INFO("sp_fma_count=%lld, dense_fma_count=%lld, other_fma_count=%lld", sp_fma_count, dense_fma_count,
             other_fma_count);
    LOG_INFO("spgemm time: %lf ms. gemm time: %lf ms", spgemm_time * 1000.0, gemm_time * 1000.0);
}

void block_parallel_factor(L2Matrix *l2) {
    LOG_INFO("分解开始");
    double factor_time = omp_get_wtime();
    LOG_INFO("分解 elapsed time: %lf ms", (omp_get_wtime() - factor_time) * 1000.0);
}
