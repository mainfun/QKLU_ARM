//
// Created by mainf on 2024/11/29.
//
#include <omp.h>
// #include <algorithm>
#include "symbol_calc.h"

void prune_graph_size(const INDEX_TYPE *Lp_start, const INDEX_TYPE *Lp_end, const INDEX_TYPE *Li,
                      const INDEX_TYPE *Up_start, const INDEX_TYPE *Up_end, const INDEX_TYPE *Ui,
                      INDEX_TYPE *lp_prune, INDEX_TYPE *up_prune,
                      INDEX_TYPE n) {
    // lp_prune[0] = 0;
    // up_prune[0] = 0;
    // for (INDEX_TYPE i = 1; i < n + 1; ++i) {
    //     if (Lp_start[i] == Lp_end[i]) {
    //         lp_prune[i] = lp_prune[i - 1];
    //     } else {
    //         INDEX_TYPE idx = Li[Lp_start[i]];
    //         lp_prune[i] = i - idx + lp_prune[i - 1];
    //     }
    //     //printf("%lld,", lp_prune[i]);
    // }
    // for (INDEX_TYPE i = 1; i < n+1; ++i) {
    //     if (Up_start[i] == Up_end[i]) {
    //         up_prune[i] = up_prune[i - 1];
    //     } else {
    //         INDEX_TYPE idx = Ui[Up_start[i]];
    //         up_prune[i] = i - idx + up_prune[i - 1];
    //     }
    //     //printf("%lld,", up_prune[i]);
    // }
    for (int i = 0; i < n+1; ++i) {
        up_prune[i]=i*10000;
        lp_prune[i]=i*10000;
    }
}

// Lp, Li: 下三角矩阵的列指针和行索引 (CSC格式)
// Up, Ui: 上三角矩阵的行指针和列索引 (CSR格式)
// n: 矩阵的大小
void prune_graph(const INDEX_TYPE *Lp_start, const INDEX_TYPE *Lp_end, const INDEX_TYPE *Li,
                 const INDEX_TYPE *Up_start, const INDEX_TYPE *Up_end, const INDEX_TYPE *Ui,
                 INDEX_TYPE *li_prune, INDEX_TYPE *li_top,
                 INDEX_TYPE *ui_prune, INDEX_TYPE *ui_top,
                 INDEX_TYPE n, INDEX_TYPE i) {
    INDEX_TYPE j = Lp_start[i], k = Up_start[i] + 1;
    INDEX_TYPE row_idx = Li[j], col_idx = Ui[k];
    while (j < Lp_end[i] && k < Up_end[i]) {
        if (col_idx == row_idx) {
            li_prune[li_top[row_idx]++] = i;
            ui_prune[ui_top[col_idx]++] = i;
            return;
        } else if (col_idx < row_idx) {
            ui_prune[ui_top[col_idx]++] = i;
            col_idx = Ui[++k];
        } else {
            li_prune[li_top[row_idx]++] = i;
            row_idx = Li[++j];
        }
    }
    for (; j < Lp_end[i]; ++j) {
        row_idx = Li[j];
        li_prune[li_top[row_idx]++] = i;
    }
    for (; k < Up_end[i]; ++k) {
        col_idx = Ui[k];
        ui_prune[ui_top[col_idx]++] = i;
    }
}
int compare(const void *a, const void *b) {
    return (*(INDEX_TYPE*)a - *(INDEX_TYPE*)b);
}

void static_symbol_calc(SparseMatrix *A, INDEX_TYPE *l_count,INDEX_TYPE *u_count,
                        INDEX_TYPE *Rp_l, INDEX_TYPE *Ri_l,
                        INDEX_TYPE *Rp_u, INDEX_TYPE *Ri_u) {
    double time = omp_get_wtime();
    INDEX_TYPE count=0;
    //--------------------------------初始化--------------------------------
    INDEX_TYPE n = A->num_row;
    INDEX_TYPE *Li = A->row_indices;
    INDEX_TYPE *Ui = A->col_indices;
    INDEX_TYPE *Lp_start = get_diag_index_v2(A->col_pointers, A->row_indices, n);
    INDEX_TYPE *Up_start = get_diag_index_v2(A->row_pointers, A->col_indices, n);
    INDEX_TYPE *Lp_end = A->col_pointers + 1;
    INDEX_TYPE *Up_end = A->row_pointers + 1;
    INDEX_TYPE *lp_prune = (INDEX_TYPE *) lu_calloc((n + 1), sizeof(INDEX_TYPE));
    INDEX_TYPE *up_prune = (INDEX_TYPE *) lu_calloc((n + 1), sizeof(INDEX_TYPE));
    prune_graph_size(A->row_pointers, Up_start, Ui,
                     A->col_pointers, Lp_start, Li, lp_prune, up_prune, n);
    INDEX_TYPE *li_top = (INDEX_TYPE *) lu_malloc(lp_prune[n] * sizeof(INDEX_TYPE)); //li_prune_每行大小
    INDEX_TYPE *ui_top = (INDEX_TYPE *) lu_malloc(up_prune[n] * sizeof(INDEX_TYPE)); //ui_prune_每列大小
    for (INDEX_TYPE i = 0; i < n; ++i) {
        li_top[i] = lp_prune[i];
        ui_top[i] = up_prune[i];
    }
    INDEX_TYPE *li_prune = (INDEX_TYPE *) lu_malloc(lp_prune[n] * sizeof(INDEX_TYPE));
    INDEX_TYPE *ui_prune = (INDEX_TYPE *) lu_malloc(up_prune[n] * sizeof(INDEX_TYPE));
    INDEX_TYPE *mark_l = (INDEX_TYPE *) lu_calloc(lp_prune[n], sizeof(INDEX_TYPE));
    INDEX_TYPE *mark_u = (INDEX_TYPE *) lu_calloc(up_prune[n], sizeof(INDEX_TYPE));
    *l_count = 0, *u_count = 0;
    Rp_l[0] = Rp_u[0] = 0;
    LOG_INFO("malloc elapsed time: %lf ms", (omp_get_wtime() - time) * 1000.0);
    //--------------------------------end 初始化--------------------------------
    for (INDEX_TYPE i = 0; i < n; ++i) {
        //复制L
        for (INDEX_TYPE j = Lp_start[i] + 1; j < Lp_end[i]; ++j) {
            INDEX_TYPE row_idx = Li[j];
            Ri_l[(*l_count)++] = row_idx;
            mark_l[row_idx] = i;
            count++;
        }
        //复制U
        for (INDEX_TYPE j = Up_start[i]; j < Up_end[i]; ++j) {
            INDEX_TYPE col_idx = Ui[j];
            Ri_u[(*u_count)++] = col_idx;
            mark_u[col_idx] = i;
            count++;
        }
        //L并依赖
        // printf("%lld L : ", i);
        for (INDEX_TYPE j = up_prune[i]; j < ui_top[i]; ++j) {
            INDEX_TYPE col_idx = ui_prune[j]; //依赖列
            // printf("%lld, ", col_idx);
            //列求并集
            for (INDEX_TYPE k = Rp_l[col_idx]; k < Rp_l[col_idx + 1]; ++k) {
                INDEX_TYPE row_idx = Ri_l[k];
                //printf("(%lld)", row_idx);
                count++;
                if (row_idx <= i) continue;
                if (mark_l[row_idx] != i) {
                    mark_l[row_idx] = i;
                    Ri_l[(*l_count)++] = row_idx;
                    // printf("(%lld)", row_idx);
                }
            }
        }
        // printf("\n");
        Rp_l[i + 1] = *l_count;
        // printf("L[%lld]:", i);
        // for (INDEX_TYPE k = Rp_l[i]; k < Rp_l[i + 1]; ++k) {
        //     INDEX_TYPE row_idx = Ri_l[k];
        //     printf("(%lld)", row_idx);
        // }
        // printf("\n");
        //U并依赖
        // printf("%lld U : ", i);
        for (INDEX_TYPE j = lp_prune[i]; j < li_top[i]; ++j) {
            INDEX_TYPE row_idx = li_prune[j]; //依赖行
            // printf("%lld, ", row_idx);
            //行求并集
            for (INDEX_TYPE k = Rp_u[row_idx]; k < Rp_u[row_idx + 1]; ++k) {
                INDEX_TYPE col_idx = Ri_u[k];
                count++;
                if (col_idx < i) continue;
                if (mark_u[col_idx] != i) {
                    mark_u[col_idx] = i;
                    Ri_u[(*u_count)++] = col_idx;
                    // printf("(%lld)", col_idx);
                }
            }
        }
        // printf("\n");
        Rp_u[i + 1] = *u_count;
        // printf("U[%lld]:", i);
        // for (INDEX_TYPE k = Rp_u[i]; k < Rp_u[i + 1]; ++k) {
        //     INDEX_TYPE col_idx = Ri_u[k];
        //     printf("(%lld)", col_idx);
        // }
        // printf("\n");
        prune_graph(Rp_l, Rp_l + 1, Ri_l, Rp_u, Rp_u + 1, Ri_u,
                    li_prune, li_top, ui_prune, ui_top, n, i);
    }

    for (INDEX_TYPE i = 0; i < n; ++i) {
        qsort(&Ri_l[Rp_l[i]], Rp_l[i + 1] - Rp_l[i], sizeof(INDEX_TYPE), compare);
    }
    for (INDEX_TYPE i = 0; i < n; ++i) {
        qsort(&Ri_u[Rp_u[i]], Rp_u[i + 1] - Rp_u[i], sizeof(INDEX_TYPE), compare);
    }
    //--------------------------------free--------------------------------
    lu_free(mark_l);
    lu_free(mark_u);
    lu_free(li_top);
    lu_free(ui_top);
    lu_free(li_prune);
    lu_free(ui_prune);
    lu_free(lp_prune);
    lu_free(up_prune);
    LOG_INFO("static_symb_calc elapsed time: %lf ms", (omp_get_wtime() - time) * 1000.0);
    LOG_DEBUG("symb calc count:%lld",count);
}
