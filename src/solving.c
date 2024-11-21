//
// Created by mainf on 2024/8/17.
//

#include "solving.h"

#include <layer_matrix.h>

#include "check.h"

void lower_solver(const SparseMatrix *L, ELE_TYPE *x, const ELE_TYPE *b) {
    INDEX_TYPE *Lp = L->col_pointers;
    INDEX_TYPE *Li = L->row_indices;
    ELE_TYPE *Lx = L->csc_values;
    INDEX_TYPE n = L->num_row;
    for (INDEX_TYPE i = 0; i < n; ++i) {
        x[i] = 0;
    }
    for (INDEX_TYPE i = 0; i < n; ++i) {
        x[i] = b[i];
    }
    for (INDEX_TYPE j = 0; j < n; ++j) {
        if (x[j] == 0) continue;
        for (INDEX_TYPE i = Lp[j]; i < Lp[j + 1]; i++) {
            x[Li[i]] -= Lx[i] * x[j];
        }
    }
}

void lower_solver_csr(const SparseMatrix *L, ELE_TYPE *x, const ELE_TYPE *b) {
    INDEX_TYPE *Lp = L->row_pointers;
    INDEX_TYPE *Li = L->col_indices;
    ELE_TYPE *Lx = L->csr_values;
    INDEX_TYPE n = L->num_row;
    for (INDEX_TYPE i = 0; i < n; ++i) {
        x[i] = b[i];
    }
    for (INDEX_TYPE i = 0; i < n; ++i) {
        for (INDEX_TYPE j = Lp[i]; j < Lp[i + 1]; j++) {
            x[i] -= Lx[j] * x[Li[j]];
        }
    }
}

void lower_solver_block(const L2Matrix *l2, ELE_TYPE *x, const ELE_TYPE *b, INDEX_TYPE n) {
    for (INDEX_TYPE i = 0; i < n; ++i) {
        x[i] = b[i];
    }
    for (int i = 0; i < l2->num_col_block; ++i) {
        //Spmv
        for (int j = l2->row_pointers[i]; j < l2->diag_index[i + 1]; j++) {
            BlockMatrix m = *l2->block_matrices[j];
            //块内CSR
            for (int ii = 0; ii < BLOCK_SIDE; ++ii) {
                for (int jj = m.row_pointers[ii]; jj < m.row_pointers[ii + 1]; jj++) {
                    x[ii] -= m.values[jj] * x[m.col_indices[jj]];
                }
            }
        }
        //解对角块
        BlockMatrix bm = *get_diag_block(l2, i);
        for (int ii = 0; ii < BLOCK_SIDE; ++ii) {
            for (int jj = bm.row_pointers[ii]; jj < l2->diag_index[ii]; jj++) {
                x[ii] -= bm.values[jj] * x[bm.col_indices[jj]];
            }
        }
    }
}

void lower_solver_block_v2(const L2Matrix *l2, ELE_TYPE *x, const ELE_TYPE *b) {
    int n = BLOCK_SIDE;
    for (int i = 0; i < l2->num_col_block; ++i) {
        //Spmv
        for (int j = l2->row_pointers[i]; j < l2->diag_index[i] + 1; j++) {
            const int block_col_idx = l2->col_indices[j];
            const BlockMatrix *bm = get_block(l2, i, block_col_idx);
            //块内CSR
            if (bm->format == SPARSE || bm->format == DENSE) {
                for (int c = 0; c < n; ++c) {
                    INDEX_TYPE big_col = (INDEX_TYPE) block_col_idx * n + c;
                    INDEX_TYPE big_ow = (INDEX_TYPE) i * n; //big_row = big_ow + r
                    //if (x[big_col] == 0) continue;
                    for (int jj = bm->col_pointers[c]; jj < bm->col_pointers[c + 1]; jj++) {
                        const int r = bm->row_indices[jj];
                        x[big_ow + r] -= bm->values[bm->offset[r] + c] * x[big_col];
                    }
                }
            } else {
                //对角
                INDEX_TYPE big_ow = (INDEX_TYPE) i * n; //big_row = big_ow + r
                if (block_col_idx == i) {
                    for (int r = 0; r < n; r++) {
                        for (int c = 0; c < r; c++) {
                            INDEX_TYPE big_col = (INDEX_TYPE) block_col_idx * n + c;
                            x[big_ow + r] -= bm->values[r * n + c] * x[big_col];
                        }
                    }
                } else {
                    for (int r = 0; r < n; r++) {
                        for (int c = 0; c < n; ++c) {
                            INDEX_TYPE big_col = (INDEX_TYPE) block_col_idx * n + c;
                            x[big_ow + r] -= bm->values[r * n + c] * x[big_col];
                        }
                    }
                }
            }
        }
    }
}

void upper_solver_block(const L2Matrix *l2, ELE_TYPE *x, const ELE_TYPE *b) {
    int n = BLOCK_SIDE;
    for (int i = l2->num_col_block - 1; i >= 0; --i) {
        for (int j = l2->diag_index[i] + 1; j < l2->row_pointers[i + 1]; j++) {
            const int block_col_idx = l2->col_indices[j];
            const BlockMatrix *bm = get_block(l2, i, block_col_idx);
            if (bm->format == SPARSE || bm->format == DENSE) {
                //SPMV
                for (int ii = 0; ii < BLOCK_SIDE; ++ii) {
                    INDEX_TYPE big_row = (INDEX_TYPE) i * n + ii;
                    INDEX_TYPE big_ol = (INDEX_TYPE) block_col_idx * n;
                    for (int jj = bm->row_pointers[ii]; jj < bm->row_pointers[ii + 1]; jj++) {
                        const int c = bm->col_indices[jj];
                        x[big_row] -= bm->values[bm->offset[ii] + c] * x[big_ol + c];
                    }
                    //x[big_row] /= bm->values[bm->offset[ii] + bm->row_pointers[ii]];
                }
            }
        }
        //解对角块
        BlockMatrix *bm = get_diag_block(l2, i);
        n=bm->side;
        if (bm->format == SPARSE || bm->format == DENSE) {
            const int block_col_idx = i;
            for (int ii = bm->side-1; ii >= 0; --ii) {
                INDEX_TYPE big_row = (INDEX_TYPE) i * n + ii;
                //if(big_row>=36057)continue;
                INDEX_TYPE big_ol = (INDEX_TYPE) block_col_idx * n;
                for (int jj = bm->row_pointers[ii] + 1; jj < bm->row_pointers[ii + 1]; jj++) {
                    const int c = bm->col_indices[jj];
                    x[big_row] -= bm->values[bm->offset[ii] + c] * x[big_ol + c];
                }
                x[big_row] /= bm->values[bm->offset[ii] + bm->col_indices[bm->row_pointers[ii]]];
            }
        } else {
        }
    }
}

/**
 * 解UX=b，U是CSC格式
 */
void upper_solver_csc(const SparseMatrix *U, ELE_TYPE *x, const ELE_TYPE *b) {
    INDEX_TYPE *Up = U->col_pointers;
    INDEX_TYPE *Ui = U->row_indices;
    ELE_TYPE *Ux = U->csc_values;
    INDEX_TYPE n = U->num_row;
    for (INDEX_TYPE i = 0; i < n; ++i) {
        x[i] = 0;
    }
    for (INDEX_TYPE i = 0; i < n; ++i) {
        x[i] = b[i];
    }
    printf("\n");
    for (INDEX_TYPE j = n - 1; j >= 0; --j) {
        //        if (x[j] == 0) continue;
        //if (fabs(Ux[Up[j + 1] - 1]) < 1e-8) Ux[Up[j + 1] - 1] = 1e-8;
        x[j] /= Ux[Up[j + 1] - 1];
        //        printf("\t\tX[%d]= %lf/%lf\n", j, x[j], Ux[Up[j + 1] - 1]);
        INDEX_TYPE i;
        for (i = Up[j]; i < Up[j + 1] - 1; i++) {
            x[Ui[i]] -= Ux[i] * x[j];
            //            printf("X[%d]-=U[%d][%d]*X[%d] (%.3lf*%.3lf=%.3lf)\n",Ui[i],Ui[i],j,j,Ux[i],x[j],Ux[i] * x[j]);
        }
    }
}

/**
 * 解UX=b，U是CSR格式
 */
void upper_solver_csr(const SparseMatrix *U, ELE_TYPE *x, const ELE_TYPE *b) {
    INDEX_TYPE *Up = U->row_pointers;
    INDEX_TYPE *Ui = U->col_indices;
    ELE_TYPE *Ux = U->csr_values;
    INDEX_TYPE n = U->num_row;
    for (INDEX_TYPE i = 0; i < n; ++i) {
        x[i] = b[i];
    }
    for (INDEX_TYPE i = n - 1; i >= 0; --i) {
        for (INDEX_TYPE j = Up[i] + 1; j < Up[i + 1]; j++) {
            x[i] -= Ux[j] * x[Ui[j]];
        }
        x[i] /= Ux[Up[i]];
    }
}

/**
 * b'= Pc Pr Dr b
 * 解y：L'y=b'
 */
void forward_substitution(const SparseMatrix *L, const INDEX_TYPE *reorder_iperm,
                          const INDEX_TYPE *reorder_perm, const ELE_TYPE *Dr,
                          const INDEX_TYPE *mc64_perm, const INDEX_TYPE *mc64_iperm,
                          const ELE_TYPE *b, ELE_TYPE *y, INDEX_TYPE n) {
    ELE_TYPE *b_prime = (ELE_TYPE *) malloc(n * sizeof(ELE_TYPE));
    for (INDEX_TYPE i = 0; i < n; i++) {
        b_prime[reorder_iperm[mc64_perm[i]]] = b[i] * Dr[i];
    }
    lower_solver_csr(L, y, b_prime);
    free(b_prime);
}

void forward_substitution_block(const L2Matrix *L, const INDEX_TYPE *reorder_iperm,
                                const INDEX_TYPE *reorder_perm, const ELE_TYPE *Dr,
                                const INDEX_TYPE *mc64_perm, const INDEX_TYPE *mc64_iperm,
                                const ELE_TYPE *b, ELE_TYPE *y, INDEX_TYPE n) {
    ELE_TYPE *b_prime = (ELE_TYPE *) malloc(n * sizeof(ELE_TYPE));
    for (INDEX_TYPE i = 0; i < n; i++) {
        b_prime[reorder_iperm[mc64_perm[i]]] = b[i] * Dr[i];
    }
    for (INDEX_TYPE i = 0; i < n; ++i) {
        y[i] = b_prime[i];
    }
    lower_solver_block_v2(L, y, b_prime);
    free(b_prime);
}

/**
 * U'Z=Y
 * 解x = D_c P^T Z
 */
void backward_substitution(const SparseMatrix *U, const INDEX_TYPE *reorder_perm,
                           const INDEX_TYPE *reorder_iperm, const ELE_TYPE *Dc,
                           ELE_TYPE *x, const ELE_TYPE *y, INDEX_TYPE n) {
    ELE_TYPE *z = (ELE_TYPE *) malloc(n * sizeof(ELE_TYPE));
    upper_solver_csr(U, z, y);
    //check_upper_solver(U, z, y);
    for (INDEX_TYPE i = 0; i < n; i++) {
        x[i] = z[reorder_iperm[i]] * Dc[i];
    }
    free(z);
}

void backward_substitution_block(const L2Matrix *U, const INDEX_TYPE *reorder_perm,
                                 const INDEX_TYPE *reorder_iperm, const ELE_TYPE *Dc,
                                 ELE_TYPE *x, const ELE_TYPE *y, INDEX_TYPE n) {
    ELE_TYPE *z = (ELE_TYPE *) malloc(n * sizeof(ELE_TYPE));
    for (INDEX_TYPE i = 0; i < n; ++i) {
        z[i] = y[i];
    }
    upper_solver_block(U, z, y);
    for (INDEX_TYPE i = 0; i < n; i++) {
        x[i] = z[reorder_iperm[i]] * Dc[i];
    }
    free(z);
}

// void solve(SparseMatrix *A, ELE_TYPE *x, const ELE_TYPE *b) {
//     INDEX_TYPE n = A->num_row;
//     PreprocessInfo *info = init_preprocess_info();
//     preprocess(A, info, true, true, true);
//     factor(A, info);
//     //check_lu(A, a->L, a->U);
//     ELE_TYPE *y = (ELE_TYPE *) malloc(n * sizeof(ELE_TYPE));
//     forward_substitution(info->L, info->mc64_iperm, info->mc64_perm,
//                          info->Dr, info->reorder_perm, info->reorder_iperm, b, y, n);
//     backward_substitution(info->U, info->reorder_perm, info->reorder_iperm, info->Dc, x, y, n);
// }
