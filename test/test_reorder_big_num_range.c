//
// Created by mainf on 2024/11/22.
//
#include <stdio.h>
#include <stdlib.h>

#define MAX_SIZE 1000  // 数组最大长度

// 判断一个数是否为大数的函数
int is_large_number(int num, int threshold) {
    return num > threshold;
}

int main() {
    int n; // 数组元素个数
    int arr[MAX_SIZE];
    int threshold;

    printf("请输入数组的元素个数（最多 %d 个）：", MAX_SIZE);
    scanf("%d", &n);
    if (n > MAX_SIZE) {
        printf("元素个数超过最大限制。\n");
        return 1;
    }

    printf("请输入数组的元素：\n");
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }

    printf("请输入判定大数的阈值：");
    scanf("%d", &threshold);

    // 开始处理数组，划分区间
    int i = 0;
    while (i < n) {
        if (!is_large_number(arr[i], threshold)) {
            // 小数区间的开始
            int start = i;
            while (i < n && !is_large_number(arr[i], threshold)) {
                i++;
            }
            printf("小数区间：位置 %d 到 %d\n", start, i - 1);
        } else {
            // 大数区间的开始
            int start = i;
            int count_large = 0; // 大数计数
            int total = 0;       // 总元素计数
            int j = i;
            while (j < n) {
                total++;
                if (is_large_number(arr[j], threshold)) {
                    count_large++;
                }
                double proportion = (double)count_large / total;
                if (proportion > 0.5) {
                    // 继续扩展大数区间
                    j++;
                    i = j;
                } else {
                    // 比例不满足要求，结束当前区间
                    break;
                }
            }
            printf("大数区间：位置 %d 到 %d\n", start, i - 1);
        }
    }

    return 0;
}
