//
// Created by mainf on 2024/8/17.
//

#include "solving.h"
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
void upper_solver(const SparseMatrix *U, ELE_TYPE *x, const ELE_TYPE *b) {
    INDEX_TYPE *Up = U->row_pointers;
    INDEX_TYPE *Ui = U->col_indices;
    ELE_TYPE *Ux = U->csr_values;
    INDEX_TYPE n = U->num_row;
    for (INDEX_TYPE i = 0; i < n; ++i) {
        x[i] = b[i];
//        if (b[i] > 1e5 || b[i] < -1e5) {
//            printf("%lld:  %lg\n", i, b[i]);
//        }
    }
    for (INDEX_TYPE j = n - 1; j >= 0; --j) {
        //if (x[j] == 0) continue;
        for (INDEX_TYPE i = Up[j] + 1; i < Up[j + 1]; i++) {
            ELE_TYPE t = x[j];
            x[j] -= Ux[i] * x[Ui[i]];
            //printf("%lld:  %lf = %lf-%lf*%lf\n", j, x[j], t, Ux[i], x[Ui[i]]);
            //if (j < n - 30) exit(0);
            if (isnan(x[j]) || isinf(x[j]) || x[j] > 1e5 || x[j] < -1e5) {
                printf("%lld:  %lf-%lf*%lf is NaN.\n", j, t, Ux[i], x[Ui[i]]);
                exit(0);
            }
        }
        ELE_TYPE t = x[j];
        x[j] /= Ux[Up[j]];
//        if (isnan(x[j])) {
//            printf("%lld:  %lf/%lf is NaN.\n", j, t, Ux[Up[j]]);
//        }
//        if (isinf(x[j])) {
//            printf("%lld:  %lf/%lf is inf.\n", j, t, Ux[Up[j]]);
//        }
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
//    ELE_TYPE *bd = (ELE_TYPE *) malloc(n * sizeof(ELE_TYPE));
//    for (INDEX_TYPE i = 0; i < n; i++) {
//        bd[i] = b[i] * Dr[i];
//    }
//    ELE_TYPE *temp = (ELE_TYPE *) malloc(n * sizeof(ELE_TYPE));
//    for (INDEX_TYPE i = 0; i < n; i++) {
//        temp[i] = bd[mc64_perm[i]];
//    }
//    for (INDEX_TYPE i = 0; i < n; i++) {
//        b_prime[i] = temp[reorder_perm[i]];
//    }

    lower_solver(L, y, b_prime);

//    for (int i = 0; i < n; ++i) {
//        printf("%lf\n",y[i]);
//    }
    //check_lower_solver(L, y, b_prime);
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
    upper_solver_csc(U, z, y);
    //check_upper_solver(U, z, y);
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