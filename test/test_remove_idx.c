#include <stdio.h>
#include <stdbool.h>

#define MAX_SIZE 100

// 函数声明
int simd_remove_idx(int *arr, int start, int end);

int main() {
    int arr[MAX_SIZE] = {1, 2, 3, 4, 5, 7, 8, 9};
    int start = 0; // 开始索引
    int end = 8;   // 结束索引，当前数组大小

    printf("Original array: ");
    for (int i = start; i < end; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");

    // 处理数组
    int new_size = simd_remove_idx(arr, start, end);

    printf("Processed array: ");
    for (int i = start; i < new_size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");

    return 0;
}

int simd_remove_idx(int *arr, int start, int end) {
    bool seen[MAX_SIZE] = {false}; // 记录出现过的偶数
    int j = start; // 新数组的索引

    for (int i = start; i < end; i++) {
        int num = arr[i];
        num = num & ~1; // 转换为下一个比当前奇数小的偶数

        // 如果是偶数且未出现过
        if (!seen[num]) {
            arr[j++] = num; // 保存到新数组
            seen[num] = true; // 标记为已出现
        }
    }

    return j; // 返回新数组的大小
}