//
// Created by mainf on 2024/4/15.
//
#ifndef QKLU_MALLOC_H
#define QKLU_MALLOC_H
#ifdef __cplusplus
extern "C" {
#endif

#include "stdlib.h"
#include "log.h"

void lu_free(void *p);

void *lu_malloc(size_t size);

void *lu_calloc(size_t n, size_t size);

void *lu_realloc(void *ptr, size_t size);

#ifdef __cplusplus
}
#endif
#endif //QKLU_MALLOC_H
