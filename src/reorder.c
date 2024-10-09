//
// Created by mainf on 2024/4/15.
//
#include "reorder.h"
#include <amd.h>
#include "base/sort.h"
#include "symbolic_analysis.h"

void apply_permutation(const SparseMatrix *A, SparseMatrix *R, const INDEX_TYPE iperm[]) {
    const clock_t start_time = clock();
    const INDEX_TYPE *Ap = A->row_pointers;
    const INDEX_TYPE *Ai = A->col_indices;
    const ELE_TYPE *Ax = A->csr_values;
    const INDEX_TYPE n = A->num_row;

    //求PAP^T
    //计算Rp
    for (INDEX_TYPE i = 0; i < n; i++) {
        const INDEX_TYPE row_num = Ap[i + 1] - Ap[i];
        R->row_pointers[iperm[i] + 1] = row_num;
    }
    R->row_pointers[0] = 0;
    for (INDEX_TYPE i = 0; i < n; i++) {
        R->row_pointers[i + 1] += R->row_pointers[i];
    }
    //计算Ri Rx
    for (INDEX_TYPE i = 0; i < n; i++) {
        const INDEX_TYPE r = iperm[i];
        const INDEX_TYPE r_start = R->row_pointers[r];
        const INDEX_TYPE r_end = R->row_pointers[r + 1];
        INDEX_TYPE a_index = Ap[i];
        for (INDEX_TYPE j = r_start; j < r_end; j++) {
            R->col_indices[j] = iperm[Ai[a_index]];
            R->csr_values[j] = Ax[a_index++];
        }
        const INDEX_TYPE num = r_end - r_start;
        //todo 排序优化
        bubbleSort(R->col_indices + r_start, R->csr_values + r_start, num);
    }
    LOG_INFO("AMD apply permutation elapsed time: %lf ms", ((double) (clock() - start_time)) / CLOCKS_PER_SEC * 1000.0);
}

INDEX_TYPE *reorder_csr_amd(const SparseMatrix *matrix, PreprocessInfo *preprocess_info) {
    clock_t start_time = clock();
    if (matrix == NULL) LOG_ERROR("matrix is NULL");
    const INDEX_TYPE n = matrix->num_row;
    const INDEX_TYPE *Ap = matrix->row_pointers;
    const INDEX_TYPE *Ai = matrix->col_indices;
    INDEX_TYPE *perm = (INDEX_TYPE *) lu_calloc(n, sizeof(INDEX_TYPE));
    INDEX_TYPE *iperm = (INDEX_TYPE *) lu_malloc(n * sizeof(INDEX_TYPE));
    // 调用AMD进行重排序
    double Info[AMD_INFO];
    INDEX_TYPE status;
    if(sizeof(INDEX_TYPE) == sizeof(double)) {
        status = amd_l_order(n, (int64_t*)Ap, (int64_t*)Ai, (int64_t*)perm, NULL, Info);
    }else {
        status = amd_order((int)n, (int*)Ap, (int*)Ai, (int*)perm, NULL, Info);
    }
    //转置矩阵
    for (int i = 0; i < n; ++i) {
        iperm[perm[i]] = i;
    }
    if (status == AMD_OK) {
        LOG_DEBUG("AMD ordering succeeded!\n");
    } else {
        LOG_DEBUG("AMD ordering failed! Status code: %lld\n", status);
        if (Info[AMD_STATUS] == AMD_INVALID) {
            LOG_ERROR("AMD ordering: Invalid matrix (jumbled or duplicate entries or wrong dimensions)\n");
        } else if (Info[AMD_STATUS] == AMD_OUT_OF_MEMORY) {
            LOG_ERROR("AMD ordering: Out of memory\n");
        } else {
            LOG_ERROR("AMD ordering: Unknown failure\n");
        }
    }
    amd_l_info(Info);
    preprocess_info->max_nz = Info[AMD_LNZ];
    preprocess_info->reorder_perm = perm;
    preprocess_info->reorder_iperm = iperm;
    LOG_DEBUG("\nL max_nz:%f .exclude dialog\n", Info[AMD_LNZ]);
    LOG_INFO("AMD reorder elapsed time: %lf ms", ((double) (clock() - start_time)) / CLOCKS_PER_SEC * 1000.0);
    return perm;
}