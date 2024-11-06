//
// Created by mainf on 2024/5/6.
//
#include <algorithm>
#include "symbolic_analysis.h"
#include "base/matrix.h"

void symbolic_calc_sym(SparseMatrix *matrix, PreprocessInfo *info) {
    clock_t symbolic_calc_time = clock();
    LOG_INFO("symbolic calc start");
    //CSR->CSC
    if (matrix->col_pointers == nullptr) matrix->col_pointers = matrix->row_pointers;
    if (matrix->row_indices == nullptr) matrix->row_indices = matrix->col_indices;
    //CSC->CSR
    if (matrix->row_pointers == nullptr) matrix->row_pointers = matrix->col_pointers;
    if (matrix->col_indices == nullptr) matrix->col_indices = matrix->row_indices;
    get_diag_index(matrix, info);
    create_row_etree_sym(matrix, info);
    get_l_col_patterns(matrix, info);
    get_l_row_patterns(matrix, info);

    for (INDEX_TYPE i = 0; i < matrix->num_row; ++i) {
        INDEX_TYPE col_start = info->L->col_pointers[i];
        INDEX_TYPE col_end = info->L->col_pointers[i + 1];
        std::sort(&(info->L->row_indices[col_start]), &(info->L->row_indices[col_end]));
    }
    for (INDEX_TYPE i = 0; i < matrix->num_row; ++i) {
        INDEX_TYPE start = info->L->row_pointers[i];
        INDEX_TYPE end = info->L->row_pointers[i + 1];
        std::sort(&(info->L->col_indices[start]), &(info->L->col_indices[end]));
    }
    info->U = init_sparse_matrix(info->L->num_row, info->L->num_col, info->L->nnz);
    info->U->row_pointers = info->L->col_pointers;
    info->U->col_indices = info->L->row_indices;
    info->U->csr_values = info->L->csc_values;
    LOG_INFO("symbolic calc elapsed time: %lf ms", ((double) (clock() - symbolic_calc_time)) / CLOCKS_PER_SEC * 1000.0);\
}

INDEX_TYPE binary_search_l(const INDEX_TYPE *array, INDEX_TYPE start, INDEX_TYPE end, INDEX_TYPE target) {
    INDEX_TYPE mid;
    while (start <= end) {
        mid = start + (end - start) / 2;
        if (array[mid] == target)
            return mid;
        else if (array[mid] < target)
            start = mid + 1;
        else
            end = mid - 1;
    }
    LOG_WARN("without diag element!");
    return -1; // 表示未找到
}

int binary_search_i(const int *array, int start, int end, int target) {
    int mid;
    while (start <= end) {
        mid = start + (end - start) / 2;
        if (array[mid] == target)
            return mid;
        else if (array[mid] < target)
            start = mid + 1;
        else
            end = mid - 1;
    }
    LOG_ERROR("without diag element!");
    return -1; // 表示未找到
}

///求对角元素索引
void get_diag_index(SparseMatrix *A, PreprocessInfo *preprocess_info) {
    if (A->row_indices == nullptr || A->col_indices == nullptr ||
        A->row_pointers == nullptr || A->col_pointers == nullptr) {
        LOG_ERROR("CSR or CSC data is null");
    }
    auto *diag_index_csr = (INDEX_TYPE *) lu_malloc(A->num_row * sizeof(INDEX_TYPE));
    auto *diag_index_csc = (INDEX_TYPE *) lu_malloc(A->num_col * sizeof(INDEX_TYPE));
    //#pragma omp parallel for
    for (INDEX_TYPE i = 0; i < A->num_row; i++) {
        diag_index_csc[i] = binary_search_l(A->row_indices, A->col_pointers[i], A->col_pointers[i + 1] - 1, i);
        diag_index_csr[i] = binary_search_l(A->col_indices, A->row_pointers[i], A->row_pointers[i + 1] - 1, i);
    }
    preprocess_info->diag_index_csc = diag_index_csc;
    preprocess_info->diag_index_csr = diag_index_csr;
}

///求对角元素索引
INDEX_TYPE *get_diag_index_v2(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai, INDEX_TYPE n) {
    INDEX_TYPE *diag_index = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    //#pragma omp parallel for
    for (INDEX_TYPE i = 0; i < n; i++) {
        diag_index[i] = binary_search_l(Ai, Ap[i], Ap[i + 1] - 1, i);
    }
    return diag_index;
}

///求对角元素索引
int *get_diag_index_i(const int *Ap, const int *Ai, int n) {
    int *diag_index = (int *) lu_malloc(n * sizeof(int));
    //#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        diag_index[i] = binary_search_i(Ai, Ap[i], Ap[i + 1] - 1, i);
    }
    return diag_index;
}

//行复制原则，用CSR格式的下三角生成
void create_row_etree_sym(SparseMatrix *csr_mat, PreprocessInfo *preprocess_info) {
    INDEX_TYPE n = csr_mat->num_row;
    Forest *forest = init_forest(n);
    auto *ancestor = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    INDEX_TYPE root;
    for (INDEX_TYPE r = 0; r < n; ++r) {
        ancestor[r] = -1;
        for (INDEX_TYPE j = csr_mat->row_pointers[r]; j < preprocess_info->diag_index_csr[r]; ++j) {
            INDEX_TYPE c = csr_mat->col_indices[j];
            root = c;
            // 如果有祖先，并且祖先不是本节点，设置祖先为本节点。
            // 例如：8找到1，7，且有ancestor[1] = 7。
            // 修改为ancestor[1] = 8;root=7;
            while (ancestor[root] != -1 && ancestor[root] != r) {
                INDEX_TYPE l = ancestor[root];
                ancestor[root] = r;
                root = l;
            }
            //ancestor[7] = 8;
            //parent[7] = 8;
            if (ancestor[root] == -1) {
                ancestor[root] = r;
                //添加
                add_child(forest, r, root);
                //printf("父：%lld  子：%lld\n",r+1,root+1);
            }
        }
    }
    lu_free(ancestor);
    preprocess_info->row_etree = forest;
}

void get_l_col_patterns(const SparseMatrix *csr_mat, PreprocessInfo *preprocess_info) {
    INDEX_TYPE n = csr_mat->num_row;
    preprocess_info->L = init_csc_matrix(csr_mat->num_row, csr_mat->num_col, long(preprocess_info->max_nz) + n);
    Forest *forest = preprocess_info->row_etree;
    INDEX_TYPE *Lp = preprocess_info->L->col_pointers;
    INDEX_TYPE *Li = preprocess_info->L->row_indices;
    auto *mark = (INDEX_TYPE *) lu_calloc(n, sizeof(INDEX_TYPE));
    INDEX_TYPE l_col_count = 0;
    Lp[0] = 0;
    for (INDEX_TYPE i = 0; i < n; ++i) {
        //复制
        for (INDEX_TYPE j = preprocess_info->diag_index_csc[i] + 1; j < csr_mat->col_pointers[i + 1]; ++j) {
            INDEX_TYPE row = csr_mat->row_indices[j];
            Li[l_col_count++] = row;
            mark[row] = i;
        }
        //查孩子
        INDEX_TYPE child = forest->first_child[i];
        while (child != -1) {
            //printf("Child of node %lld: %lld\n", i + 1, child + 1);
            //求并集
            for (INDEX_TYPE k = Lp[child + 1] - 1; k >= Lp[child]; k--) {
                INDEX_TYPE row = Li[k];
                if (row <= i) continue;
                //printf("访问：%lld，", row);
                if (mark[row] != i) {
                    mark[row] = i;
                    //todo 有序插入(归并树)
                    Li[l_col_count++] = row;
                }
            }
            child = forest->next_sibling[child];
        }
        Lp[i + 1] = l_col_count;
    }
    preprocess_info->L->nnz = l_col_count;
    lu_free(mark);
}

void get_l_row_patterns(const SparseMatrix *csr_mat, PreprocessInfo *info) {
    INDEX_TYPE n = csr_mat->num_row;
    auto max_nz = (INDEX_TYPE) info->max_nz;
    SparseMatrix *L;
    if (info->L == NULL) {
        L = init_csr_matrix(csr_mat->num_row, csr_mat->num_col, max_nz);
    } else {
        L = info->L;
        L->row_pointers = (INDEX_TYPE *) lu_malloc((csr_mat->num_row + 1) * sizeof(INDEX_TYPE));
        L->col_indices = (INDEX_TYPE *) lu_malloc(max_nz * sizeof(INDEX_TYPE));
        L->csr_values = (ELE_TYPE *) lu_malloc(max_nz * sizeof(ELE_TYPE));
    }
    INDEX_TYPE *parent = info->row_etree->parent;
    INDEX_TYPE *Lp = L->row_pointers;
    INDEX_TYPE *Li = L->col_indices;
    auto *mark = (INDEX_TYPE *) lu_calloc(n, sizeof(INDEX_TYPE));
    Lp[0] = 0;
    INDEX_TYPE l_row_count = 0;
    for (INDEX_TYPE r = 0; r < n; ++r) {
        mark[r] = r; //标记是否被添加过
        for (INDEX_TYPE k = csr_mat->row_pointers[r]; k < info->diag_index_csr[r]; ++k) {
            INDEX_TYPE j = csr_mat->col_indices[k];
            while (mark[j] != r) {
                mark[j] = r;
                Li[l_row_count++] = j;
                j = parent[j];
            }
        }
        Lp[r + 1] = l_row_count;
    }
    free(mark);
}

/**-------------------------------------------------------------------------**/
/**---------------------------------非对称-----------------------------------**/
/**-------------------------------------------------------------------------**/
INDEX_TYPE pruneL(INDEX_TYPE jcol,INDEX_TYPE *U_r_idx,INDEX_TYPE *U_c_ptr,INDEX_TYPE *L_r_idx,INDEX_TYPE *L_c_ptr,
                  INDEX_TYPE *work_space,INDEX_TYPE *prune_space) {
    if (jcol == 0)
        return 0;

    INDEX_TYPE min = U_c_ptr[jcol];
    INDEX_TYPE max = U_c_ptr[jcol + 1];
    INDEX_TYPE cmin, cmax, crow, doprune;
    doprune = 0;

    for (INDEX_TYPE i = min; i < max; i++) //遍历U的一列
    {
        doprune = 0;
        crow = U_r_idx[i]; //U的行号
        cmin = L_c_ptr[crow]; //前端矩阵的列
        cmax = prune_space[crow]; //修剪后的列尾
        for (INDEX_TYPE j = cmin; j < cmax; j++) {
            if (L_r_idx[j] == jcol) {
                doprune = 1;
                break;
            }
        }

        if (doprune == 1) {
            for (INDEX_TYPE j = cmin; j < cmax; j++) {
                INDEX_TYPE ccrow = L_r_idx[j];
                if (ccrow > jcol && work_space[ccrow] == jcol) //加入L的结构，但又>jcol，进行移除
                {
                    INDEX_TYPE temp = L_r_idx[cmax - 1];
                    L_r_idx[cmax - 1] = L_r_idx[j];
                    L_r_idx[j] = temp;
                    cmax--;
                    j--;
                }
            }
        }
        prune_space[crow] = cmax;
    }
    return 0;
}

void fill_in_2_no_sort_pruneL(INDEX_TYPE n, INDEX_TYPE nnz,
                              INDEX_TYPE *ai, INDEX_TYPE *ap,
                              INDEX_TYPE **L_rowpointer, INDEX_TYPE **L_columnindex,
                              INDEX_TYPE **U_rowpointer, INDEX_TYPE **U_columnindex) {
    INDEX_TYPE relloc_zjq = nnz;
    INDEX_TYPE *U_r_idx = (INDEX_TYPE *) lu_malloc(relloc_zjq * sizeof(INDEX_TYPE)); //exclude diagonal
    INDEX_TYPE *U_c_ptr = (INDEX_TYPE *) lu_malloc((n + 1) * sizeof(INDEX_TYPE));
    INDEX_TYPE *L_r_idx = (INDEX_TYPE *) lu_malloc(relloc_zjq * sizeof(INDEX_TYPE)); //include diagonal
    INDEX_TYPE *L_c_ptr = (INDEX_TYPE *) lu_malloc((n + 1) * sizeof(INDEX_TYPE));
    U_c_ptr[0] = 0;
    L_c_ptr[0] = 0;

    INDEX_TYPE *parent = (INDEX_TYPE *) lu_malloc((n + 1) * sizeof(INDEX_TYPE)); //for dfs
    INDEX_TYPE *xplore = (INDEX_TYPE *) lu_malloc((n + 1) * sizeof(INDEX_TYPE));

    INDEX_TYPE *work_space = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE)); //use this to avoid sorting
    INDEX_TYPE *prune_space = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));

    INDEX_TYPE U_maxsize = relloc_zjq;
    INDEX_TYPE L_maxsize = relloc_zjq;

    INDEX_TYPE U_size = 0;
    INDEX_TYPE L_size = 0; //record quantity

    INDEX_TYPE row = -1;
    INDEX_TYPE oldrow = -1;
    INDEX_TYPE xdfs = -1;
    INDEX_TYPE maxdfs = -1;
    INDEX_TYPE kchild = -1;
    INDEX_TYPE kpar = -1;

    for (INDEX_TYPE k = 0; k < n; k++) {
        work_space[k] = -1; //avoid conflict
        parent[k] = -1;
        xplore[k] = 0;
    }

    //列
    for (INDEX_TYPE i = 0; i < n; i++) {
        //列大小
        INDEX_TYPE n_rows = ap[i + 1] - ap[i];

        for (INDEX_TYPE k = 0; k < n_rows; k++) {
            row = (ai + ap[i])[k]; //行号
            if (work_space[row] == i)
                continue;

            work_space[row] = i;
            if (row >= i) //L部分
            {
                L_r_idx[L_size] = row;
                L_size++;

                if (L_size >= L_maxsize - 100) {
                    L_r_idx = (INDEX_TYPE *) lu_realloc(L_r_idx, (L_maxsize + relloc_zjq) * sizeof(INDEX_TYPE));
                    L_maxsize = L_maxsize + nnz;
                }
            } else {
                U_r_idx[U_size] = row;
                U_size++;
                if (U_size >= U_maxsize - 100) {
                    U_r_idx = (INDEX_TYPE *) realloc(U_r_idx, (U_maxsize + relloc_zjq) * sizeof(INDEX_TYPE));
                    U_maxsize = U_maxsize + nnz;
                }
                //U得L的结构
                //do dfs
                oldrow = -1;
                parent[row] = oldrow;
                xdfs = L_c_ptr[row]; //xdfs ~ maxdfs
                maxdfs = prune_space[row]; //prune
                do {
                    /* code */
                    while (xdfs < maxdfs) {
                        kchild = L_r_idx[xdfs];
                        xdfs++;

                        if (work_space[kchild] != i) {
                            work_space[kchild] = i;
                            if (kchild >= i) {
                                L_r_idx[L_size] = kchild;
                                L_size++;
                                if (L_size >= L_maxsize - 100) {
                                    L_r_idx = (INDEX_TYPE *) realloc(
                                        L_r_idx, (L_maxsize + relloc_zjq) * sizeof(INDEX_TYPE));
                                    L_maxsize = L_maxsize + nnz;
                                }
                            } else {
                                U_r_idx[U_size] = kchild;
                                U_size++;
                                if (U_size >= U_maxsize - 100) {
                                    U_r_idx = (INDEX_TYPE *) realloc(
                                        U_r_idx, (U_maxsize + relloc_zjq) * sizeof(INDEX_TYPE));
                                    U_maxsize = U_maxsize + nnz;
                                }

                                xplore[row] = xdfs;
                                oldrow = row;
                                row = kchild;
                                parent[row] = oldrow;
                                xdfs = L_c_ptr[row];
                                maxdfs = prune_space[row]; //prune
                            }
                        }
                    }

                    kpar = parent[row];

                    if (kpar == -1)
                        break;

                    row = kpar;
                    xdfs = xplore[row];
                    maxdfs = prune_space[row];
                } while (kpar != -1);
            }
        }
        U_c_ptr[i + 1] = U_size;
        L_c_ptr[i + 1] = L_size;
        prune_space[i] = L_size;

        pruneL(i, U_r_idx, U_c_ptr, L_r_idx, L_c_ptr, work_space, prune_space);
    }

    //L
    for (INDEX_TYPE i = 0; i < n; ++i) {
        std::sort(&L_r_idx[L_c_ptr[i]], &L_r_idx[L_c_ptr[i + 1]]);
    }
    //U
    for (INDEX_TYPE i = 0; i < n; ++i) {
        std::sort(&U_r_idx[U_c_ptr[i]], &U_r_idx[U_c_ptr[i + 1]]);
    }

    lu_free(parent);
    lu_free(xplore);
    lu_free(work_space);
    lu_free(prune_space);

    *L_rowpointer = L_c_ptr;
    *L_columnindex = L_r_idx;
    *U_rowpointer = U_c_ptr;
    *U_columnindex = U_r_idx;
}

void csr2csc_pattern(SparseMatrix *m) {
    // 初始化
    m->col_pointers = (INDEX_TYPE *) lu_calloc(m->num_col + 1, sizeof(INDEX_TYPE));
    m->row_indices = (INDEX_TYPE *) lu_malloc(m->nnz * sizeof(INDEX_TYPE));
    // 计算每列的非零元素个数
    for (INDEX_TYPE i = 0; i < m->nnz; i++) {
        m->col_pointers[m->col_indices[i] + 1]++;
    }
    // 转换为列指针的累加形式
    for (INDEX_TYPE i = 0; i < m->num_col; i++) {
        m->col_pointers[i + 1] += m->col_pointers[i];
    }
    for (INDEX_TYPE i = 0; i < m->num_row; i++) {
        for (INDEX_TYPE j = m->row_pointers[i]; j < m->row_pointers[i + 1]; j++) {
            INDEX_TYPE col = m->col_indices[j];
            INDEX_TYPE dst = m->col_pointers[col];
            m->row_indices[dst] = i;
            m->col_pointers[col]++;
        }
    }
    // 修正列指针
    for (INDEX_TYPE i = m->num_col; i > 0; i--) {
        m->col_pointers[i] = m->col_pointers[i - 1];
    }
    m->col_pointers[0] = 0;
}

void csc2csr_pattern(SparseMatrix *m) {
    // 初始化
    m->row_pointers = (INDEX_TYPE *) lu_calloc(m->num_row + 1, sizeof(INDEX_TYPE));
    m->col_indices = (INDEX_TYPE *) lu_malloc(m->nnz * sizeof(INDEX_TYPE));
    // 计算每列的非零元素个数
    for (INDEX_TYPE i = 0; i < m->nnz; i++) {
        m->row_pointers[m->row_indices[i] + 1]++;
    }
    // 转换为累加形式
    for (INDEX_TYPE i = 0; i < m->num_row; i++) {
        m->row_pointers[i + 1] += m->row_pointers[i];
    }
    for (INDEX_TYPE i = 0; i < m->num_col; i++) {
        for (INDEX_TYPE j = m->col_pointers[i]; j < m->col_pointers[i + 1]; j++) {
            INDEX_TYPE col = m->row_indices[j];
            INDEX_TYPE dst = m->row_pointers[col];
            m->col_indices[dst] = i;
            m->row_pointers[col]++;
        }
    }
    // 修正列指针
    for (INDEX_TYPE i = m->num_row; i > 0; i--) {
        m->row_pointers[i] = m->row_pointers[i - 1];
    }
    m->row_pointers[0] = 0;
}

void csr2csc_pattern_v2(const int *row_pointers, const int *col_indices,
                        int *col_pointers, int *row_indices,
                        int nnz, int n) {
    // 计算每列的非零元素个数
    for (int i = 0; i < nnz; i++) {
        col_pointers[col_indices[i] + 1]++;
    }
    // 转换为列指针的累加形式
    for (int i = 0; i < n; i++) {
        col_pointers[i + 1] += col_pointers[i];
    }
    for (int i = 0; i < n; i++) {
        for (int j = row_pointers[i]; j < row_pointers[i + 1]; j++) {
            int col = col_indices[j];
            int dst = col_pointers[col];
            row_indices[dst] = i;
            col_pointers[col]++;
        }
    }
    // 修正列指针
    for (int i = n; i > 0; i--) {
        col_pointers[i] = col_pointers[i - 1];
    }
    col_pointers[0] = 0;
}

void csc2csr_pattern_v2(int *row_pointers, int *col_indices,
                        const int *col_pointers, const int *row_indices,
                        int nnz, int n) {
    // 计算每列的非零元素个数
    for (int i = 0; i < nnz; i++) {
        row_pointers[row_indices[i] + 1]++;
    }
    // 转换为累加形式
    for (int i = 0; i < n; i++) {
        row_pointers[i + 1] += row_pointers[i];
    }
    for (int i = 0; i < n; i++) {
        for (int j = col_pointers[i]; j < col_pointers[i + 1]; j++) {
            int col = row_indices[j];
            int dst = row_pointers[col];
            col_indices[dst] = i;
            row_pointers[col]++;
        }
    }
    // 修正列指针
    for (int i = n; i > 0; i--) {
        row_pointers[i] = row_pointers[i - 1];
    }
    row_pointers[0] = 0;
}

/**
 * 方
 */
void csc2csr(SparseMatrix *m) {
    // 初始化
    m->row_pointers = (INDEX_TYPE *) lu_calloc(m->num_row + 1, sizeof(INDEX_TYPE));
    m->col_indices = (INDEX_TYPE *) lu_malloc(m->nnz * sizeof(INDEX_TYPE));
    m->csr_values = (ELE_TYPE *) lu_malloc(m->nnz * sizeof(ELE_TYPE));
    // 计算每列的非零元素个数
    for (INDEX_TYPE i = 0; i < m->nnz; i++) {
        m->row_pointers[m->row_indices[i] + 1]++;
    }
    // 转换为累加形式
    for (INDEX_TYPE i = 0; i < m->num_row; i++) {
        m->row_pointers[i + 1] += m->row_pointers[i];
    }
    for (INDEX_TYPE i = 0; i < m->num_col; i++) {
        for (INDEX_TYPE j = m->col_pointers[i]; j < m->col_pointers[i + 1]; j++) {
            INDEX_TYPE col = m->row_indices[j];
            INDEX_TYPE dst = m->row_pointers[col];
            ELE_TYPE v = m->csc_values[j];
            m->col_indices[dst] = i;
            m->csr_values[dst] = v;
            m->row_pointers[col]++;
        }
    }
    // 修正列指针
    for (INDEX_TYPE i = m->num_row; i > 0; i--) {
        m->row_pointers[i] = m->row_pointers[i - 1];
    }
    m->row_pointers[0] = 0;
}

void csr2csc(SparseMatrix *m) {
    // 初始化
    m->col_pointers = (INDEX_TYPE *) lu_calloc(m->num_row + 1, sizeof(INDEX_TYPE));
    m->row_indices = (INDEX_TYPE *) lu_malloc(m->nnz * sizeof(INDEX_TYPE));
    m->csc_values = (ELE_TYPE *) lu_malloc(m->nnz * sizeof(ELE_TYPE));
    // 计算每列的非零元素个数
    for (INDEX_TYPE i = 0; i < m->nnz; i++) {
        m->col_pointers[m->col_indices[i] + 1]++;
    }
    // 转换为列指针的累加形式
    for (INDEX_TYPE i = 0; i < m->num_col; i++) {
        m->col_pointers[i + 1] += m->col_pointers[i];
    }
    for (INDEX_TYPE i = 0; i < m->num_row; i++) {
        for (INDEX_TYPE j = m->row_pointers[i]; j < m->row_pointers[i + 1]; j++) {
            INDEX_TYPE col = m->col_indices[j];
            INDEX_TYPE dst = m->col_pointers[col];
            ELE_TYPE v = m->csr_values[j];
            m->row_indices[dst] = i;
            m->csc_values[dst] = v;
            m->col_pointers[col]++;
        }
    }
    // 修正列指针
    for (INDEX_TYPE i = m->num_col; i > 0; i--) {
        m->col_pointers[i] = m->col_pointers[i - 1];
    }
    m->col_pointers[0] = 0;
}


void get_pattern_asym(SparseMatrix *m, PreprocessInfo *info) {
    //csc l, csr u
    csr2csc_pattern(m);
    get_diag_index(m, info);
    INDEX_TYPE *sym_pointers = (INDEX_TYPE *) lu_calloc(m->num_row, sizeof(INDEX_TYPE)); //对称点
    for (INDEX_TYPE i = 0; i < m->num_row; i++) {
        //找对称点
        INDEX_TYPE l_j = info->diag_index_csc[i];
        for (INDEX_TYPE u_j = info->diag_index_csr[i]; u_j < m->row_pointers[i + 1]; u_j++) {
            INDEX_TYPE col = m->col_indices[u_j];
            INDEX_TYPE row = m->row_indices[l_j];
            if (row == col) {
                sym_pointers[i] = row;
                break;
            } else if (row > col) {
                l_j++;
            }
        }
        //求并集
    }
}
