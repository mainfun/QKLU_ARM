//
// Created by mainf on 2024/8/19.
//

#include <stdbool.h>
#include "identify_chip.h"
#include "stdio.h"
#include "log.h"

// 识别架构的函数，返回枚举类型
CPU_ARCHITECTURE identify_architecture() {
#ifdef __aarch64__
    return ARCH_ARM64;
#elif defined(__arm__)
    return ARCH_ARM32;
#elif defined(__x86_64__) || defined(_M_X64)
    return ARCH_X86_64;
#elif defined(__i386) || defined(_M_IX86)
    return ARCH_X86_32;
#else
    return ARCH_UNKNOWN;
#endif
}

// 输出架构信息的函数
const char* get_cpu_architecture_name(CPU_ARCHITECTURE arch) {
    switch (arch) {
        case ARCH_ARM32:
            return "32-bit ARM architecture (ARM)";
        case ARCH_ARM64:
            return "64-bit ARM architecture (AArch64)";
        case ARCH_X86_32:
            return "32-bit x86 architecture (x86)";
        case ARCH_X86_64:
            return "64-bit x86 architecture (x86_64)";
        default:
            return "Unknown architecture";
    }
}

// 检查Neon是否支持的SIMD128
bool supports_neon128() {
#if defined(__ARM_NEON)
    return 1;
#else
    return 0;
#endif
}

// 定义cpuid指令函数 仅x86平台
void cpuid(int registers[4], int eax, int ecx) {
    #if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    __asm__ volatile ("cpuid"
            : "=a" (registers[0]), "=b" (registers[1]), "=c" (registers[2]), "=d" (registers[3])
            : "a" (eax), "c" (ecx));
    #endif
}

// 检查是否支持AVX512的函数
bool supports_avx512() {
    int registers[4];
    // 调用cpuid指令，使用EAX=7和ECX=0获取扩展特性信息
    cpuid(registers, 7, 0);
    // AVX512F支持在EBX寄存器的第16位
    return (registers[1] & (1 << 16)) != 0;
}

// 检查是否支持128位SIMD (SSE, SSE2) 的函数
bool supports_sse_sse2() {
    int registers[4];

    // 使用 cpuid，eax=1 来获取特性信息
    cpuid(registers, 1, 0);

    // 检查 edx 寄存器中的位
    // SSE 支持在 edx 寄存器的第 25 位 (bit 25)
    // SSE2 支持在 edx 寄存器的第 26 位 (bit 26)
    bool sse_supported = (registers[3] & (1 << 25)) != 0;
    bool sse2_supported = (registers[3] & (1 << 26)) != 0;

    return sse_supported && sse2_supported;
}

void log_cpu_info(){
    CPU_ARCHITECTURE arch = identify_architecture();
    LOG_INFO("Detected architecture: %s", get_cpu_architecture_name(arch));
    switch (arch) {
        case ARCH_ARM32:;
        case ARCH_ARM64:
            #ifndef __ARM_NEON
            LOG_INFO("not support neon 128");
            #else
            LOG_INFO("support neon 128");
            #endif
            break;
        case ARCH_X86_32:;
        case ARCH_X86_64:
            if(!supports_avx512())
                LOG_INFO("not support avx 512");
            else
                LOG_INFO("support avx 512");
            if(!supports_sse_sse2())
                LOG_INFO("not support sse sse2");
            else
                LOG_INFO("support sse sse2");
            break;
        case ARCH_UNKNOWN:
            LOG_INFO("ARCH_UNKNOWN");
            break;
    }
}