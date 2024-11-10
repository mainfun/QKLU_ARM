//
// Created by mainf on 2024/10/9.
//

#include "layer_matrix.h"
#include <symbolic_analysis.h>
#include <base/base_math.h>

#define DENSE_THRESHOLD (0.5 * BLOCK_SIDE * BLOCK_SIDE)

BlockMatrix *init_BlockMatrix() {
    BlockMatrix *sub_matrix = (BlockMatrix *) lu_calloc(1, sizeof(BlockMatrix));
    return sub_matrix;
}

L2Matrix *init_L2Matrix() {
    L2Matrix *matrix = (L2Matrix *) lu_calloc(1, sizeof(L2Matrix));
    return matrix;
}

void free_BlockMatrix(BlockMatrix *matrix) {
    if (matrix) {
        lu_free(matrix->row_indices);
        lu_free(matrix->row_pointers);
        lu_free(matrix->col_indices);
        lu_free(matrix->values);
        lu_free(matrix);
    }
}

void free_L2Matrix(L2Matrix *m) {
}

BlockMatrix *get_diag_block(const L2Matrix *l2, const int i) {
    return l2->block_matrices[i * l2->num_row_block + i];
}

BlockMatrix *get_block(const L2Matrix *l2, const int i, const int j) {
    return l2->block_matrices[i * l2->num_row_block + j];
}

void print_csr(BlockMatrix *mat, int n, int d_size) {
    printf("%s\n", mat->format==SPARSE?"SPARSE":"DENSE");
    printf("%d x %d\n", mat->num_row, mat->num_col);
    printf("Row ptr: ");
    for (int i = 0; i <= n; i++) {
        printf("%d ", mat->row_pointers[i]);
    }
    printf("\nCol indices: ");
    for (int i = 0; i < n; i++) {
        for (int j = mat->row_pointers[i]; j < mat->row_pointers[i + 1]; j++) {
            printf("%d ", mat->col_indices[j]);
        }
        printf(",");
    }
    printf("\nCSR values: ");
    for (int i = 0; i < d_size; i++) {
        printf("%lf ", mat->values[i]);
    }
    printf("\nOffset: ");
    for (int i = 0; i < n; i++) {
        printf("%d ", mat->offset[i]);
    }
    printf("\n");
}

void print_csc(BlockMatrix *mat, int n, int d_size) {
    printf("%s\n", mat->format==SPARSE?"SPARSE":"DENSE");
    printf("%d x %d\n", mat->num_row, mat->num_col);
    printf("Col ptr: ");
    for (int i = 0; i <= n; i++) {
        printf("%d ", mat->col_pointers[i]);
    }
    printf("\nRow indices: ");
    for (int i = 0; i < n; i++) {
        for (int j = mat->col_pointers[i]; j < mat->col_pointers[i + 1]; j++) {
            printf("%d ", mat->row_indices[j]);
        }
        printf(",");
    }
    printf("\nCSC values: ");
    for (int i = 0; i < d_size; i++) {
        printf("%lf ", mat->values[i]);
    }
    printf("\nOffset: ");
    for (int i = 0; i < n; i++) {
        printf("%d ", mat->offset[i]);
    }
    printf("\n");
}

void print_block(L2Matrix *mat) {
    for (int i = 0; i < mat->num_row_block; i++) {
        for (int j = 0; j < mat->num_col_block; j++) {
            printf("\nSubmatrix (%d, %d):\n", i, j);
            BlockMatrix *bm = mat->block_matrices[i * mat->num_col_block + j];
            //if (bm->format == SPARSE) {
                if (i >= j)
                    print_csc(bm, mat->block_width, bm->d_size);
                if (i <= j)
                    print_csr(bm, mat->block_width, bm->d_size);
            //}
        }
    }
}

typedef struct {
    int *row_pointers;
    int *col_indices;
    int nnz;
} temp_csr_l;

int *get_offset(const int *row_pointers, const int *col_indices, int BLOCK_SIDE, int *d_size) {
    //计算offset
    int *offset = (int *) lu_malloc(BLOCK_SIDE * sizeof(int));
    for (int k = 0; k < BLOCK_SIDE; ++k) {
        int start = row_pointers[k];
        int end = row_pointers[k + 1];
        offset[k] = *d_size;
        //如果本行非空
        if (end != start) {
            *d_size += col_indices[end - 1] - col_indices[start] + 1;
        }
    }
    for (int r = 0; r < BLOCK_SIDE; ++r) {
        //todo 空行的offset是冗余的
        offset[r] -= col_indices[row_pointers[r]];
    }
    return offset;
}

int *get_offset_diag(const int *row_ptr_l, const int *col_idx_l,
                     const int *row_ptr_u, const int *col_idx_u,
                     int BLOCK_SIDE, int *d_size) {
    //计算offset
    int *offset = (int *) lu_malloc(BLOCK_SIDE * sizeof(int));
    for (int k = 0; k < BLOCK_SIDE; ++k) {
        int start_l = row_ptr_l[k];
        int end_l = row_ptr_l[k + 1];
        //int start_u = row_ptr_u[k];
        int end_u = row_ptr_u[k + 1];
        offset[k] = *d_size;
        //(U永远非空，L可能是空行)
        int col1 = start_l == end_l ? k : col_idx_l[start_l];
        int col2 = col_idx_u[end_u - 1];
        *d_size += col2 - col1 + 1;
    }
    for (int r = 0; r < BLOCK_SIDE; ++r) {
        //如果L空行
        if (row_ptr_l[r] == row_ptr_l[r + 1]) {
            offset[r] -= r;
        } else {
            offset[r] -= col_idx_l[row_ptr_l[r]];
        }
    }
    return offset;
}

//symbolic_blocking A
//calc offset for each block
//numerical_blocking A
void symbolic_blocking(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai,
                       L2Matrix *l2, const int BLOCK_SIDE, INDEX_TYPE n) {
}

void numerical_blocking() {
}

void csr2L2Matrix(const INDEX_TYPE *Lp, const INDEX_TYPE *Li,
                  const INDEX_TYPE *Up, const INDEX_TYPE *Ui,
                  const INDEX_TYPE *Ap, const INDEX_TYPE *Ai, const ELE_TYPE *Ax,
                  L2Matrix *l2, const int BLOCK_SIDE, INDEX_TYPE n) {
    clock_t block_time = clock();
    //---------------------get block nnz---------------------
    const int block_height = BLOCK_SIDE;
    const int block_width = BLOCK_SIDE;
    int num_row_block = CEIL_DIV(n, block_height);
    int num_col_block = CEIL_DIV(n, block_width);
    LOG_DEBUG("num_row_block: %d, num_col_block: %d", num_row_block, num_col_block);

    int *block_nnz_arr = (int *) lu_calloc(num_col_block * num_row_block, sizeof(int));
    for (INDEX_TYPE i = 0; i < n; ++i) {
        //i / block_height 是块的行号，Ai[j] / block_width是块的列号
        INDEX_TYPE block_row_ptr = (i / block_height) * num_col_block;
        for (INDEX_TYPE j = Lp[i]; j < Lp[i + 1]; ++j) {
            block_nnz_arr[block_row_ptr + Li[j] / block_width]++;
        }
        for (INDEX_TYPE j = Up[i]; j < Up[i + 1]; ++j) {
            block_nnz_arr[block_row_ptr + Ui[j] / block_width]++;
        }
    }
    int block_count = 0;
    int dense_count = 0;
    for (int i = 0; i < num_row_block * num_col_block; ++i) {
        if (block_nnz_arr[i] > 0) {
            block_count++;
            if (block_nnz_arr[i] > DENSE_THRESHOLD) dense_count++;
            //printf("%d ",block_nnz_arr[i]);
        }
    }
    LOG_DEBUG("block_count: %d", block_count);
    LOG_DEBUG("dense_count: %d", dense_count);
    if (block_count < 500) {
        LOG_DEBUG("block nnz: \n");
        for (int i = 0; i < num_row_block; ++i) {
            for (int j = 0; j < num_col_block; ++j) {
                printf("%4d ", block_nnz_arr[i * num_col_block + j]);
            }
            printf("\n");
        }
    }
    l2->num_row_block = num_row_block;
    l2->num_col_block = num_col_block;
    l2->block_width = block_width;
    l2->block_height = block_height;
    l2->block_count = block_count;
    //--------------------end get block nnz--------------------
    //分配BlockMatrix
    BlockMatrix *block_matrices = (BlockMatrix *) lu_calloc(block_count, sizeof(BlockMatrix));
    BlockMatrix **block_matrices_ptr = (BlockMatrix **)
            lu_calloc(num_col_block * num_row_block, sizeof(BlockMatrix *));
    temp_csr_l *csr_l_arr = lu_calloc(block_count, sizeof(temp_csr_l));
    temp_csr_l **csr_l_ptr = (temp_csr_l **)
            lu_calloc(num_col_block * num_row_block, sizeof(temp_csr_l *));
    l2->block_matrices = block_matrices_ptr;
    /**---------------分配L2Matrix---------------**/
    int block_index = 0;
    for (int i = 0; i < num_row_block; i++) {
        //U
        for (int j = 0; j < num_col_block; ++j) {
            const int INDEX = i * num_row_block + j;
            int nnz = block_nnz_arr[INDEX];
            if (nnz > 0) {
                BlockMatrix *bm = &block_matrices[block_index];
                bm->format = nnz > DENSE_THRESHOLD ? DENSE : SPARSE;
                //L
                if (j <= i) {
                    bm->col_pointers = (int *) lu_calloc((block_width + 1), sizeof(int));
                    bm->row_indices = (int *) lu_malloc(nnz * sizeof(int));

                    temp_csr_l *t = &csr_l_arr[block_index];
                    t->row_pointers = (int *) lu_calloc((block_height + 1), sizeof(int));
                    t->col_indices = (int *) lu_malloc(nnz * sizeof(int));
                    csr_l_ptr[INDEX] = t;
                }
                //U
                if (j >= i) {
                    bm->row_pointers = (int *) lu_calloc((block_height + 1), sizeof(int));
                    bm->col_indices = (int *) lu_malloc(nnz * sizeof(int));
                }
                block_matrices_ptr[INDEX] = bm;
                block_index++;
            }
        }
    }
    /**--------------end分配空间---------------**/
    /**--------------计算非零索引---------------**/
    for (INDEX_TYPE i = 0; i < n; i++) {
        INDEX_TYPE row_block_idx = i / block_height;
        int local_row = (int) (i % block_height);
        for (INDEX_TYPE row_nz_idx = Lp[i]; row_nz_idx < Lp[i + 1]; row_nz_idx++) {
            INDEX_TYPE col_idx = Li[row_nz_idx];
            INDEX_TYPE col_block_idx = col_idx / block_width;
            // 获取该元素对应的子块矩阵
            const INDEX_TYPE INDEX = row_block_idx * num_col_block + col_block_idx;
            BlockMatrix *bm = block_matrices_ptr[INDEX];
            //printf("(%lld,%lld)", row_block_idx, col_block_idx);
            // 计算该子块中的行和列相对位置
            int local_col = (int) (col_idx % block_width);
            // 添加元素到子块中

            csr_l_ptr[INDEX]->row_pointers[local_row + 1]++;
            csr_l_ptr[INDEX]->col_indices[csr_l_ptr[INDEX]->nnz] = local_col;
            csr_l_ptr[INDEX]->nnz++;
        }
        for (INDEX_TYPE row_nz_idx = Up[i]; row_nz_idx < Up[i + 1]; row_nz_idx++) {
            INDEX_TYPE col_idx = Ui[row_nz_idx];
            INDEX_TYPE col_block_idx = col_idx / block_width;
            BlockMatrix *bm = block_matrices_ptr[row_block_idx * num_col_block + col_block_idx];
            // 计算该子块中的行和列相对位置
            int local_col = (int) (col_idx % block_width);
            // 添加元素到子块中

            bm->row_pointers[local_row + 1]++;
            bm->col_indices[bm->u_nnz] = local_col;
            bm->u_nnz++;
        }
    }
    long long sum_d_size = 0;
    for (int i = 0; i < num_row_block; ++i) {
        // 更新每个子块的行列指针数组
        for (int j = 0; j < num_col_block; ++j) {
            BlockMatrix *bm = block_matrices_ptr[i * num_col_block + j];
            if (bm == NULL) continue;

            int d_size = 0;
            const temp_csr_l *t = csr_l_ptr[i * num_col_block + j];
            //L
            if (i > j) {
                for (int r = 1; r <= BLOCK_SIDE; r++) {
                    t->row_pointers[r] += t->row_pointers[r - 1];
                }
                bm->l_nnz = t->nnz;
                csr2csc_pattern_v2(t->row_pointers, t->col_indices,
                                   bm->col_pointers, bm->row_indices, t->nnz, BLOCK_SIDE);
                if (bm->format == SPARSE)
                    bm->offset = get_offset(t->row_pointers, t->col_indices, BLOCK_SIDE, &d_size);
            }
            //U
            if (i < j) {
                for (int r = 1; r <= BLOCK_SIDE; r++) {
                    bm->row_pointers[r] += bm->row_pointers[r - 1];
                }
                if (bm->format == SPARSE)
                    bm->offset = get_offset(bm->row_pointers, bm->col_indices, BLOCK_SIDE, &d_size);
            }
            //对角
            if (i == j) {
                bm->l_nnz = t->nnz;
                for (int r = 1; r <= BLOCK_SIDE; r++) {
                    t->row_pointers[r] += t->row_pointers[r - 1];
                    bm->row_pointers[r] += bm->row_pointers[r - 1];
                }
                csr2csc_pattern_v2(t->row_pointers, t->col_indices,
                                   bm->col_pointers, bm->row_indices, t->nnz, BLOCK_SIDE);
                if (bm->format == SPARSE)
                    bm->offset = get_offset_diag(t->row_pointers, t->col_indices,
                                                 bm->row_pointers, bm->col_indices,
                                                 BLOCK_SIDE, &d_size);
                // d_size = BLOCK_SIDE * BLOCK_SIDE;
                // bm->offset = lu_calloc(BLOCK_SIDE * BLOCK_SIDE, sizeof(int));
            }
            //printf("(%d,%d) ", d_size,block_nnz_arr[i*num_row_block+j]);
            if (bm->format == SPARSE) {
                sum_d_size += d_size;
                bm->d_size = d_size;
            }

            if (bm->format == DENSE) { //稠密结构
                int *offset = (int *) lu_malloc(BLOCK_SIDE * sizeof(int));
                for (int k = 0; k < BLOCK_SIDE; ++k) {
                    offset[k] = k * BLOCK_SIDE;
                }
                bm->offset = offset;
                sum_d_size += BLOCK_SIDE * BLOCK_SIDE;
                bm->d_size = BLOCK_SIDE * BLOCK_SIDE;
            }
        }
    }
    block_matrices_ptr[18731 * num_row_block + 18731]->d_size = 1000;
    sum_d_size += 2000;
    LOG_INFO("\n the values size = %lld MB. (sum_d_size=%lld)", sum_d_size / 1024 / 1024 * 8, sum_d_size);
    if (sum_d_size <= 0) {
        LOG_ERROR("sum_d_size<=0");
    }
    lu_free(csr_l_ptr);
    lu_free(csr_l_arr);
    /**--------------end 计算非零索引---------------**/
    //-----------------------块索引-----------------------
    l2->row_pointers = (int *) lu_malloc((num_row_block + 1) * sizeof(int));
    l2->col_indices = (int *) lu_malloc(block_count * sizeof(int));
    l2->row_pointers[0] = 0;
    int row_count = 0;
    for (int i = 0; i < num_row_block; ++i) {
        for (int j = 0; j < num_col_block; ++j) {
            if (block_nnz_arr[i * num_row_block + j] != 0) {
                l2->col_indices[row_count++] = j;
            }
        }
        l2->row_pointers[i + 1] = row_count;
    }
    l2->diag_index = get_diag_index_i(l2->row_pointers, l2->col_indices, l2->num_row_block);
    lu_free(block_nnz_arr);
    //-----------------------end块索引-----------------------
    // 值数组初始化
    ELE_TYPE *values = (ELE_TYPE *) lu_calloc(sum_d_size, sizeof(ELE_TYPE *));
    long long values_index = 0;
    for (int i = 0; i < num_row_block; ++i) {
        for (int j = l2->row_pointers[i]; j < l2->row_pointers[i + 1]; ++j) {
            BlockMatrix *bm = block_matrices_ptr[i * num_row_block + l2->col_indices[j]];
            if (bm != NULL) {
                if (bm->d_size < 0) {
                    LOG_DEBUG("bm->d_size<0 at (%d,%d)", i, l2->col_indices[j]);
                }
                bm->values = &values[values_index];
                values_index += bm->d_size;
            }
        }
    }
    // 遍历原始矩阵并填充子块值
    for (INDEX_TYPE i = 0; i < n; i++) {
        INDEX_TYPE row_block_idx = i / block_height;
        for (INDEX_TYPE row_nz_idx = Ap[i]; row_nz_idx < Ap[i + 1]; row_nz_idx++) {
            INDEX_TYPE col_idx = Ai[row_nz_idx];
            INDEX_TYPE col_block_idx = col_idx / block_width;
            // 获取该元素对应的子块矩阵
            BlockMatrix *sub_matrix = block_matrices_ptr[row_block_idx * num_col_block + col_block_idx];
            // 计算该子块中的行和列相对位置
            INDEX_TYPE local_row = i % block_height;
            INDEX_TYPE local_col = col_idx % block_width;
            // 添加元素到子块中
            //if (sub_matrix->format == SPARSE) {
            sub_matrix->values[sub_matrix->offset[local_row] + local_col] = Ax[row_nz_idx];
            // } else {
            //     sub_matrix->values[local_row * block_width + local_col] = Ax[row_nz_idx];
            // }
        }
    }
    if (l2->block_count < 100)print_block(l2);
    LOG_INFO("blocking elapsed time: %lf ms", ((double) (clock() - block_time)) / CLOCKS_PER_SEC * 1000.0);
}
