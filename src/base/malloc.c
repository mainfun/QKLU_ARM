//
// Created by mainf on 2024/8/14.
//
#include "malloc.h"

void lu_free(void *ptr) {
    if (ptr != NULL) free(ptr), ptr = NULL;
}

void *aligned_malloc(size_t size) {
    void *ptr;
    if (posix_memalign(&ptr, 64, size) != 0) {
        return NULL; // 分配失败
    }
    return ptr;
}

void *lu_malloc(size_t size) {
    void *ptr = malloc(size);
    if (ptr == NULL)
        LOG_ERROR("内存分配失败\n");
    return ptr;
}

void *lu_calloc(size_t n, size_t size) {
    void *ptr = calloc(n, size);
    if (ptr == NULL)
        LOG_ERROR("内存分配失败\n");
    return ptr;
}

void *lu_realloc(void *ptr, size_t size) {
    void *tmp = realloc(ptr, size);
    if (tmp == NULL) {
        LOG_ERROR("重新分配内存失败\n");
    }
    return tmp;
}
