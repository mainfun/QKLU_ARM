#include <numerical.h>
#include <preprocess.h>
#include <cblas.h>

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
            int b_row = m2->row_indices[j];
            ELE_TYPE l = B(b_row, i) /= A(i, i);
            for (int k = m1->row_pointers[i] + 1; k < m1->row_pointers[i + 1]; k++) {
                int a_col = m1->col_indices[k];
                B(b_row, a_col) -= l * A(i, a_col);
            }
        }
    }
}

void DTSTRF_DD(const ELE_TYPE *a, ELE_TYPE *b) {
    for (int i = 0; i < BLOCK_SIDE; i++) {
        for (int j = 0; j < BLOCK_SIDE; j++) {
            //第j行消第i行
            ELE_TYPE l = DB(i, j) /= DA(j, j);
            if (l == 0) continue;
            for (int k = j + 1; k < BLOCK_SIDE; k++) {
                DB(i, k) -= l * DA(j, k);
            }
        }
    }
}

//U稀疏
void DTSTRF_SD(const BlockMatrix *m1, ELE_TYPE *b) {
    ELE_TYPE *a = m1->values;
    int *offset1 = m1->offset;
    for (int i = 0; i < BLOCK_SIDE; i++) {
        for (int j = 0; j < BLOCK_SIDE; ++j) {
            DB(j, i) /= A(i, i);
        }
        for (int k = m1->row_pointers[i] + 1; k < m1->row_pointers[i + 1]; k++) {
            int a_col = m1->col_indices[k];
            ELE_TYPE a_i_a_col = A(i, a_col);
            for (int j = 0; j < BLOCK_SIDE; ++j) {
                DB(j, a_col) -= DB(j, i) * a_i_a_col;
            }
        }
    }
}

//U稠密
void DTSTRF_DS(const ELE_TYPE *a, BlockMatrix *m2) {
    ELE_TYPE *b = m2->values;
    int *offset2 = m2->offset;
    for (int i = 0; i < BLOCK_SIDE; i++) {
        for (int j = m2->col_pointers[i]; j < m2->col_pointers[i + 1]; ++j) {
            int b_row = m2->row_indices[j];
            ELE_TYPE l = B(b_row, i) /= DA(i, i);
            for (int k = i + 1; k < BLOCK_SIDE; k++) {
                B(b_row, k) -= l * DA(i, k);
            }
        }
    }
}

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
            }
        }
    }
}

void DGESSM_DD(const ELE_TYPE *a, ELE_TYPE *b) {
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int j = i - 1; j < BLOCK_SIDE; ++j) {
            ELE_TYPE l = DA(j, i);
            if (l == 0) continue;
            for (int c = 0; c < BLOCK_SIDE; ++c) {
                DB(j, c) -= l * DB(i, c);
            }
        }
    }
}

///稀疏L解稠密矩阵
void DGESSM_DS(const ELE_TYPE *lx, BlockMatrix *m1) {
    ELE_TYPE *a = m1->values;
    int *offset1 = m1->offset;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int j = i; j < BLOCK_SIDE; j++) {
            ELE_TYPE l = DL(j, i);
            if (l == 0) continue;
            //printf("第%d行消%d行\n", i, j);
            for (int k = m1->row_pointers[i]; k < m1->row_pointers[i + 1]; k++) {
                int c = m1->col_indices[k];
                //printf("%lf -= %lf * %lf\n", A(j, c), l, A(i, c));
                A(j, c) -= l * A(i, c);
            }
        }
    }
}

void DGESSM_SD(BlockMatrix *m1, ELE_TYPE *a) {
    ELE_TYPE *lx = m1->values;
    int *offset_l = m1->offset;
    for (int i = 0; i < BLOCK_SIDE; ++i) {
        for (int p = m1->col_pointers[i]; p < m1->col_pointers[i + 1]; p++) {
            int j = m1->row_indices[p];
            ELE_TYPE l = L(j, i);
            //printf("第%d行消%d行\n", i, j);
            for (int k = 0; k < BLOCK_SIDE; k++) {
                //printf("%lf -= %lf * %lf\n", A(j, c), l, A(i, c));
                DA(j, k) -= l * DA(i, k);
            }
        }
    }
}

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


void lu_factor(BlockMatrix *m) {
    if (m->format == DENSE) {
        LUD(m->values);
    } else {
        LUS(m->col_pointers, m->col_pointers + 1, m->row_indices,
            m->row_pointers, m->row_pointers + 1, m->col_indices,
            m->values,BLOCK_SIDE, m->offset);
    }
}


void upper_solving(BlockMatrix *m1, BlockMatrix *m2) {
    if (m1->format == DENSE) {
        if (m2->format == DENSE) {
            DTSTRF_DD(m1->values, m2->values);
        } else {
            DTSTRF_DS(m1->values, m2);
        }
    } else {
        if (m2->format == DENSE) {
            DTSTRF_SD(m1, m2->values);
        } else {
            DTSTRF_SS(m1, m2);
        }
    }
}

void lower_solving(BlockMatrix *m1, BlockMatrix *m2) {
    if (m1->format == DENSE) {
        if (m2->format == DENSE) {
            DGESSM_DD(m1->values, m2->values);
        } else {
            DGESSM_DS(m1->values, m2);
        }
    } else {
        if (m2->format == DENSE) {
            DGESSM_SD(m1, m2->values);
        } else {
            DGESSM_SS(m1, m2);
        }
    }
}

///m3+=m1*m2
void schur_complement(BlockMatrix *m1, BlockMatrix *m2, BlockMatrix *m3) {
    if (m1->format == DENSE) {
        if (m2->format == DENSE) {
            GEMM(m1->values, m2->values, m3->values);
        } else {
            SpMM_DS(m1->values, m2, m3);
        }
    } else {
        if (m2->format == DENSE) {
            SpMM_SD(m1, m2->values, m3);
        } else {
            SpGEMM(m1, m2, m3);
        }
    }
}

/**
 * @details 块间是上看法，块内用右看法，CPU单线程运行
*/
void block_factor_up_looking_v0(L2Matrix *l2) {
    for (int i = 0; i < l2->num_row_block; ++i) {
        for (int p = l2->row_pointers[i]; p < l2->diag_index[i]; p++) {
            int j = l2->col_indices[p];
            //上三角解
            BlockMatrix *l_bm = get_block(l2, i, j);
            BlockMatrix *diag_block = get_diag_block(l2, j);
            upper_solving(diag_block, l_bm);
            //schur a(i,c)-=l(i,j) * u(j,c)
            for (int k = l2->diag_index[j] + 1; k < l2->row_pointers[j + 1]; k++) {
                int col_idx = l2->col_indices[k];
                BlockMatrix *bm = get_block(l2, i, col_idx);
                if (bm != NULL) {
                    schur_complement(l_bm, get_block(l2, j, col_idx), bm);
                }
            }
        }
        //LU factor
        BlockMatrix *diag_block = get_diag_block(l2, i);
        lu_factor(diag_block);
        //下三角解
        for (int p = l2->diag_index[i] + 1; p < l2->row_pointers[i + 1]; p++) {
            const int c = l2->col_indices[p];
            lower_solving(diag_block, get_block(l2, i, c));
        }
    }
}

///块间右看法
void block_factor_right_looking_v0(L2Matrix *l2) {
    for (int i = 0; i < l2->num_row_block; ++i) {
        //LU分解 A[i][i]
        BlockMatrix *diag_block = get_diag_block(l2, i);
        lu_factor(diag_block);
        for (int j = l2->diag_index[i] + 1; j < l2->row_pointers[i + 1]; ++j) {
            //L[i][i]解A[i][]
            int c = l2->col_indices[j];
            lower_solving(diag_block, get_block(l2, i, c));
        }
        //L CSC
        for (int j = l2->diag_index_csc[i] + 1; j < l2->col_pointers[i + 1]; ++j) {
            int r = l2->row_indices[j];
            BlockMatrix *bm_l = get_block(l2, r, i);
            //U[i][i]解A[][i]
            upper_solving(diag_block, bm_l);
            for (int k = l2->diag_index[i] + 1; k < l2->row_pointers[i + 1]; ++k) {
                int c = l2->col_indices[j];
                BlockMatrix *bm_A = get_block(l2, r, c);
                if (bm_A != NULL) schur_complement(bm_l, get_block(l2, i, c), bm_A);
            }
        }
    }
}

void block_factor(L2Matrix *l2) {
    LOG_INFO("分解开始");
    clock_t factor_time = clock();
    block_factor_up_looking_v0(l2);
    LOG_INFO("分解 elapsed time: %lf ms", ((double) (clock() - factor_time)) / CLOCKS_PER_SEC * 1000.0);

}
