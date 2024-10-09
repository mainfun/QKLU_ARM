//
// Created by mainf on 2024/8/19.
//
#ifndef QKLU_IDENTIFY_CHIP_H
#define QKLU_IDENTIFY_CHIP_H
#ifdef __cplusplus
extern "C" {
#endif //__cplusplus
#include <stdbool.h>
// 定义枚举类型表示不同的架构
typedef enum {
    ARCH_UNKNOWN,
    ARCH_ARM32,
    ARCH_ARM64,
    ARCH_X86_32,
    ARCH_X86_64
} CPU_ARCHITECTURE;

CPU_ARCHITECTURE identify_architecture();

const char *get_cpu_architecture_name(CPU_ARCHITECTURE arch);

bool supports_neon128();

bool supports_avx512();

bool supports_sse_sse2();

void log_cpu_info();

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //QKLU_IDENTIFY_CHIP_H
