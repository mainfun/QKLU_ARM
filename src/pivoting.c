//
// Created by mainf on 2024/5/25.
//
#include "pivoting.h"

#define MAX(A, B) (((A) > (B)) ? (A) : (B))
#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define ABS(A) ((A > 0) ? (A) : (-A))

#define MC64_FLAG (-5)

void mc64Dd(INDEX_TYPE col, INDEX_TYPE n, INDEX_TYPE *queue, const ELE_TYPE *row_scale_value, INDEX_TYPE *save_tmp) {
    INDEX_TYPE loc = save_tmp[col];
    ELE_TYPE rsv_min = row_scale_value[col];

    for (INDEX_TYPE i = 0; i < n; i++) {
        if (loc <= 0) {
            break;
        }
        INDEX_TYPE tmp_loc = (loc - 1) / 2;
        INDEX_TYPE index_loc = queue[tmp_loc];
        if (rsv_min >= row_scale_value[index_loc]) {
            break;
        }
        queue[loc] = index_loc;
        save_tmp[index_loc] = loc;
        loc = tmp_loc;
    }

    queue[loc] = col;
    save_tmp[col] = loc;
}

void mc64Ed(INDEX_TYPE *queue_length, INDEX_TYPE n, INDEX_TYPE *queue, const ELE_TYPE *row_scale_value, INDEX_TYPE *save_tmp) {
    INDEX_TYPE loc = 0;
    (*queue_length)--;
    INDEX_TYPE now_queue_length = *queue_length;
    ELE_TYPE rsv_min = row_scale_value[queue[now_queue_length]];

    for (INDEX_TYPE i = 0; i < n; i++) {
        INDEX_TYPE tmp_loc = (loc + 1) * 2 - 1;
        if (tmp_loc > now_queue_length) {
            break;
        }
        ELE_TYPE rsv_now = row_scale_value[queue[tmp_loc]];
        if (tmp_loc < now_queue_length) {
            ELE_TYPE rsv_after = row_scale_value[queue[tmp_loc + 1]];
            if (rsv_now > rsv_after) {
                rsv_now = rsv_after;
                tmp_loc++;
            }
        }
        if (rsv_min <= rsv_now) {
            break;
        }
        queue[loc] = queue[tmp_loc];
        save_tmp[queue[loc]] = loc;
        loc = tmp_loc;
    }

    queue[loc] = queue[now_queue_length];
    save_tmp[queue[now_queue_length]] = loc;
}

void mc64Fd(INDEX_TYPE loc_origin,
            INDEX_TYPE *queue_length,
            INDEX_TYPE n,
            INDEX_TYPE *queue,
            const ELE_TYPE *row_scale_value,
            INDEX_TYPE *save_tmp) {
    (*queue_length)--;
    INDEX_TYPE now_queue_length = *queue_length;

    if (loc_origin == now_queue_length) {
        return;
    }

    INDEX_TYPE loc = loc_origin;
    ELE_TYPE rsv_min = row_scale_value[queue[now_queue_length]];

    // begin mc64dd
    for (INDEX_TYPE i = 0; i < n; ++i) {
        if (loc <= 0) {
            break;
        }
        INDEX_TYPE tmp_loc = (loc - 1) / 2;
        INDEX_TYPE index_loc = queue[tmp_loc];
        if (rsv_min >= row_scale_value[index_loc]) {
            break;
        }
        queue[loc] = index_loc;
        save_tmp[index_loc] = loc;
        loc = tmp_loc;
    }

    queue[loc] = queue[now_queue_length];
    save_tmp[queue[now_queue_length]] = loc;

    // begin mc64ed
    for (INDEX_TYPE i = 0; i < n; i++) {
        INDEX_TYPE tmp_loc = (loc + 1) * 2 - 1;
        if (tmp_loc > now_queue_length) {
            break;
        }
        ELE_TYPE rsv_now = row_scale_value[queue[tmp_loc]];
        if (tmp_loc < now_queue_length) {
            ELE_TYPE rsv_after = row_scale_value[queue[tmp_loc + 1]];
            if (rsv_now > rsv_after) {
                rsv_now = rsv_after;
                tmp_loc++;
            }
        }
        if (rsv_min <= rsv_now) {
            break;
        }
        queue[loc] = queue[tmp_loc];
        save_tmp[queue[loc]] = loc;
        loc = tmp_loc;
    }

    queue[loc] = queue[now_queue_length];
    save_tmp[queue[now_queue_length]] = loc;
}

void static_pivoting(SparseMatrix *S, INDEX_TYPE **perm, INDEX_TYPE **iperm,
                     ELE_TYPE **row_scale, ELE_TYPE **col_scale) {
    clock_t start_time=clock();
    INDEX_TYPE n = S->num_row;
    INDEX_TYPE nnz = S->nnz;

    INDEX_TYPE finish_flag = 0;

    INDEX_TYPE *rowptr = S->row_pointers;
    INDEX_TYPE *colidx = S->col_indices;
    ELE_TYPE *val = S->csr_values;

    INDEX_TYPE *col_perm = (INDEX_TYPE *) malloc(sizeof(INDEX_TYPE) * n);
    INDEX_TYPE *row_perm = (INDEX_TYPE *) malloc(sizeof(INDEX_TYPE) * n);
    INDEX_TYPE *rowptr_tmp = (INDEX_TYPE *) malloc(sizeof(INDEX_TYPE) * n);
    INDEX_TYPE *queue = (INDEX_TYPE *) malloc(sizeof(INDEX_TYPE) * n);
    INDEX_TYPE *save_tmp = (INDEX_TYPE *) malloc(sizeof(INDEX_TYPE) * n);
    INDEX_TYPE *ans_queue = (INDEX_TYPE *) malloc(sizeof(INDEX_TYPE) * n);

    ELE_TYPE *fabs_value = (ELE_TYPE *) malloc(sizeof(ELE_TYPE) * nnz);
    ELE_TYPE *max_value = (ELE_TYPE *) malloc(sizeof(ELE_TYPE) * n);
    ELE_TYPE *col_scale_value = (ELE_TYPE *) malloc(sizeof(ELE_TYPE) * n);
    ELE_TYPE *row_scale_value = (ELE_TYPE *) malloc(sizeof(ELE_TYPE) * n);

    for (INDEX_TYPE i = 0; i < n; i++) {
        if (rowptr[i] >= rowptr[i + 1]) {
            LOG_ERROR("error exit zero row in matrix");
        }
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        col_perm[i] = -1;
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        row_perm[i] = -1;
    }
    for (INDEX_TYPE i = 0; i < n; i++) {
        col_scale_value[i] = DBL_MAX;
    }
    for (INDEX_TYPE i = 0; i < n; i++) {
        row_scale_value[i] = 0.0;
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        rowptr_tmp[i] = rowptr[i];
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        save_tmp[i] = -1;
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        if (rowptr[i] >= rowptr[i + 1]) {
            LOG_ERROR("error exit zero row in matrix");
        }
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        max_value[i] = 0.0;//行最大值
        for (INDEX_TYPE j = rowptr[i]; j < rowptr[i + 1]; j++) {
            fabs_value[j] = ABS(val[j]);
            max_value[i] = MAX(fabs_value[j], max_value[i]);
        }

        ELE_TYPE now_row_max = max_value[i];

        if ((now_row_max != 0.0)) {
            now_row_max = log(now_row_max);
        } else {
            LOG_ERROR("now_row_max = 0");
        }

        for (INDEX_TYPE j = rowptr[i]; j < rowptr[i + 1]; j++) {
            if ((fabs_value[j] != 0.0)) {
                fabs_value[j] = now_row_max - log(fabs_value[j]);
            } else {
                fabs_value[j] = DBL_MAX / 5.0;
            }
        }
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        for (INDEX_TYPE j = rowptr[i]; j < rowptr[i + 1]; j++) {
            INDEX_TYPE col = colidx[j];
            if (fabs_value[j] <= col_scale_value[col]) {
                col_scale_value[col] = fabs_value[j];
                col_perm[col] = i;
                save_tmp[col] = j;
            }
        }
    }
    for (INDEX_TYPE i = 0; i < n; i++) {
        if (col_perm[i] >= 0) {
            if (row_perm[col_perm[i]] < 0) {
                finish_flag++;
                col_perm[i] = col_perm[i];
                row_perm[col_perm[i]] = save_tmp[i];
            } else {
                col_perm[i] = -1;
            }
        }
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        if (finish_flag == n) {
            break;
        }
        if (row_perm[i] < 0) {
            ELE_TYPE col_max = DBL_MAX;
            INDEX_TYPE save_col = -1;
            INDEX_TYPE save_index = -1;
            for (INDEX_TYPE j = rowptr[i]; j < rowptr[i + 1]; j++) {
                INDEX_TYPE col = colidx[j];
                ELE_TYPE now_value = fabs_value[j] - col_scale_value[col];
                if (now_value > col_max) {
                    // nothing to do
                } else if ((now_value >= col_max) && (now_value != DBL_MAX)) {
                    if ((col_perm[col] < 0) && (col_perm[save_col] >= 0)) {
                        col_max = now_value;
                        save_col = col;
                        save_index = j;
                    }
                } else {
                    col_max = now_value;
                    save_col = col;
                    save_index = j;
                }
            }

            row_scale_value[i] = col_max;
            if (col_perm[save_col] < 0) {
                row_perm[i] = save_index;
                col_perm[save_col] = i;
                rowptr_tmp[i] = save_index + 1;
                finish_flag++;
            } else {
                INDEX_TYPE break_flag = 0;
                for (INDEX_TYPE j = save_index; j < rowptr[i + 1]; j++) {
                    INDEX_TYPE col = colidx[j];
                    if ((fabs_value[j] - col_scale_value[col]) <= col_max) {
                        INDEX_TYPE now_col = col_perm[col];
                        if (rowptr_tmp[now_col] < rowptr[now_col + 1]) {
                            for (INDEX_TYPE k = rowptr_tmp[now_col]; k < rowptr[now_col + 1]; k++) {
                                INDEX_TYPE tmp_col = colidx[k];
                                if (col_perm[tmp_col] < 0) {
                                    if ((fabs_value[k] - col_scale_value[tmp_col]) <= row_scale_value[now_col]) {
                                        row_perm[now_col] = k;
                                        col_perm[tmp_col] = now_col;
                                        rowptr_tmp[now_col] = k + 1;
                                        break_flag = 1;
                                        break;
                                    }
                                }
                            }
                            if (break_flag == 1) {
                                row_perm[i] = j;
                                col_perm[col] = i;
                                rowptr_tmp[i] = j + 1;
                                finish_flag++;
                                break;
                            }
                            rowptr_tmp[now_col] = rowptr[now_col + 1];
                        }
                    }
                }
            }
        }
    }
    if (finish_flag != n) {
        for (INDEX_TYPE i = 0; i < n; i++) {
            row_scale_value[i] = DBL_MAX;
        }

        for (INDEX_TYPE i = 0; i < n; i++) {
            save_tmp[i] = -1;
        }
    }

    for (INDEX_TYPE now_row = 0; now_row < n; now_row++) {
        if (finish_flag == n) {
            break;
        }
        if (row_perm[now_row] < 0) {
            INDEX_TYPE row = now_row;
            INDEX_TYPE queue_length = 0;
            INDEX_TYPE low = n;
            INDEX_TYPE top = n;
            INDEX_TYPE save_index = -1;
            INDEX_TYPE save_row = -1;

            rowptr_tmp[row] = MC64_FLAG;
            ELE_TYPE min_cost = DBL_MAX;
            ELE_TYPE sum_cost = DBL_MAX;

            for (INDEX_TYPE k = rowptr[row]; k < rowptr[row + 1]; k++) {
                INDEX_TYPE col = colidx[k];
                ELE_TYPE now_value = fabs_value[k] - col_scale_value[col];
                if (now_value < sum_cost) {
                    if (col_perm[col] < 0) {
                        sum_cost = now_value;
                        save_index = k;
                        save_row = row;
                    } else {
                        min_cost = MIN(now_value, min_cost);
                        row_scale_value[col] = now_value;
                        queue[queue_length++] = k;
                    }
                }
            }

            INDEX_TYPE now_queue_length = queue_length;
            queue_length = 0;

            for (INDEX_TYPE k = 0; k < now_queue_length; k++) {
                INDEX_TYPE queue_index = queue[k];
                INDEX_TYPE col = colidx[queue_index];
                if (row_scale_value[col] >= sum_cost) {
                    row_scale_value[col] = DBL_MAX;
                } else {
                    if (row_scale_value[col] <= min_cost) {
                        low--;
                        queue[low] = col;
                        save_tmp[col] = low;
                    } else {
                        save_tmp[col] = queue_length++;
                        mc64Dd(col, n, queue, row_scale_value, save_tmp);
                    }
                    INDEX_TYPE now_col = col_perm[col];
                    ans_queue[now_col] = queue_index;
                    rowptr_tmp[now_col] = row;
                }
            }

            for (INDEX_TYPE k = 0; k < finish_flag; k++) {
                if (low == top) {

                    if ((queue_length == 0) || (row_scale_value[queue[0]] >= sum_cost)) {
                        break;
                    }
                    INDEX_TYPE col = queue[0];
                    min_cost = row_scale_value[col];
                    do {
                        mc64Ed(&queue_length, n, queue, row_scale_value, save_tmp);
                        queue[--low] = col;
                        save_tmp[col] = low;
                        if (queue_length == 0) {
                            break;
                        }
                        col = queue[0];
                    } while (row_scale_value[col] <= min_cost);
                }
                INDEX_TYPE now_queue_length = queue[top - 1];
                if (row_scale_value[now_queue_length] >= sum_cost) {
                    break;
                }
                top--;
                row = col_perm[now_queue_length];
                ELE_TYPE row_sum_max = row_scale_value[now_queue_length] - fabs_value[row_perm[row]] +
                                       col_scale_value[now_queue_length];

                for (INDEX_TYPE k = rowptr[row]; k < rowptr[row + 1]; k++) {
                    INDEX_TYPE col = colidx[k];
                    if (save_tmp[col] < top) {
                        ELE_TYPE now_value = row_sum_max + fabs_value[k] - col_scale_value[col];
                        if (now_value < sum_cost) {
                            if (col_perm[col] < 0) {
                                sum_cost = now_value;
                                save_index = k;
                                save_row = row;
                            } else {
                                if ((row_scale_value[col] > now_value) && (save_tmp[col] < low)) {
                                    row_scale_value[col] = now_value;
                                    if (now_value <= min_cost) {
                                        if (save_tmp[col] >= 0) {
                                            mc64Fd(save_tmp[col], &queue_length, n, queue, row_scale_value, save_tmp);
                                        }
                                        low--;
                                        queue[low] = col;
                                        save_tmp[col] = low;
                                    } else {
                                        if (save_tmp[col] < 0) {
                                            save_tmp[col] = queue_length++;
                                        }
                                        mc64Dd(col, n, queue, row_scale_value, save_tmp);
                                    }

                                    INDEX_TYPE now_col = col_perm[col];
                                    ans_queue[now_col] = k;
                                    rowptr_tmp[now_col] = row;
                                }
                            }
                        }
                    }
                }
            }

            if (sum_cost != DBL_MAX) {
                finish_flag++;
                col_perm[colidx[save_index]] = save_row;
                row_perm[save_row] = save_index;
                row = save_row;

                for (INDEX_TYPE k = 0; k < finish_flag; k++) {
                    INDEX_TYPE now_rowptr_tmp = rowptr_tmp[row];
                    if (now_rowptr_tmp == MC64_FLAG) {
                        break;
                    }
                    INDEX_TYPE col = colidx[ans_queue[row]];
                    col_perm[col] = now_rowptr_tmp;
                    row_perm[now_rowptr_tmp] = ans_queue[row];
                    row = now_rowptr_tmp;
                }

                for (INDEX_TYPE k = top; k < n; k++) {
                    INDEX_TYPE col = queue[k];
                    col_scale_value[col] = col_scale_value[col] + row_scale_value[col] - sum_cost;
                }
            }

            for (INDEX_TYPE k = low; k < n; k++) {
                INDEX_TYPE col = queue[k];
                row_scale_value[col] = DBL_MAX;
                save_tmp[col] = -1;
            }

            for (INDEX_TYPE k = 0; k < queue_length; k++) {
                INDEX_TYPE col = queue[k];
                row_scale_value[col] = DBL_MAX;
                save_tmp[col] = -1;
            }
        }
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        INDEX_TYPE now_col = row_perm[i];
        if (now_col >= 0) {
            row_scale_value[i] = fabs_value[now_col] - col_scale_value[colidx[now_col]];
        } else {
            row_scale_value[i] = 0.0;
        }
        if (col_perm[i] < 0) {
            col_scale_value[i] = 0.0;
        }
    }

    if (finish_flag == n) {
        for (INDEX_TYPE i = 0; i < n; i++) {
            row_perm[col_perm[i]] = i;

            ELE_TYPE a = max_value[i];
            if (a != 0.0) {
                row_scale_value[i] -= log(a);
            } else {
                row_scale_value[i] = 0.0;
            }
        }
    } else {
        for (INDEX_TYPE i = 0; i < n; i++) {
            row_perm[i] = -1;
        }

        INDEX_TYPE ans_queue_length = 0;
        for (INDEX_TYPE i = 0; i < n; i++) {
            if (col_perm[i] < 0) {
                ans_queue[ans_queue_length++] = i;
            } else {
                row_perm[col_perm[i]] = i;
            }
        }
        ans_queue_length = 0;
        for (INDEX_TYPE i = 0; i < n; i++) {
            if (row_perm[i] < 0) {
                col_perm[ans_queue[ans_queue_length++]] = i;
            }
        }
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        col_scale_value[i] = exp(col_scale_value[i]);
    }

    for (INDEX_TYPE i = 0; i < n; i++) {
        row_scale_value[i] = exp(row_scale_value[i]);
    }

    free(rowptr_tmp);
    free(queue);
    free(save_tmp);
    free(ans_queue);

    free(fabs_value);
    free(max_value);

    *perm = row_perm;
    *iperm = col_perm;

    *col_scale = col_scale_value;
    *row_scale = row_scale_value;
    LOG_INFO("static pivoting elapsed time: %lf ms",((double) (clock() - start_time)) / CLOCKS_PER_SEC * 1000.0);
}

#undef MAX
#undef MIN
#undef ABS