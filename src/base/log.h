//
// Created by mainf on 2024/9/15.
//

#ifndef QKLU_LOG_H
#define QKLU_LOG_H
#ifdef __cplusplus
extern "C" {
#endif
#include <time.h>
#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#define SHOW_DEBUG

#ifdef SHOW_DEBUG
#define LOG_DEBUG(message, ...) fprintf(stdout, "[DEBUG]:" message "\n", ##__VA_ARGS__)
#else
#define LOG_DEBUG(message, ...) (void)0
#endif

// 定义不同级别的日志宏
#define LOG_WARN(message, ...) fprintf(stdout, "\n[WARN]:---------" message "\n", ##__VA_ARGS__)

#define LOG_INFO(message, ...) fprintf(stdout, "[INFO]:" message "\n", ##__VA_ARGS__)

#define LOG_ERROR(message, ...) do { \
        time_t now = time(NULL); \
        char *date = ctime(&now); \
        date[strlen(date) - 1] = '\0';  /* 移除换行符 */ \
        fprintf(stderr, "[ERROR] %s:%d %s --- " message "\n", __FILE__, __LINE__, date, ##__VA_ARGS__); \
        exit(EXIT_FAILURE);                             \
    } while (0)


#ifdef __cplusplus
}
#endif
#endif //QKLU_LOG_H
