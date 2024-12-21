#include <iostream>
#include <vector>
#include <cmath>
#include <random>

double sum_circle(const std::vector<std::vector<double>>& image, int x, int y, double r) {
    double total_sum = 0.0;
    int height = image.size();
    int width = image[0].size();

    // 遍历图像的每个像素
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            // 计算当前像素到圆心的距离
            double distance = std::sqrt((j - x) * (j - x) + (i - y) * (i - y));
            // 如果在圆内，则累加该像素的值
            if (distance <= r) {
                total_sum += image[i][j];
            }
        }
    }

    return total_sum;
}

int main() {
    int height = 1000;
    int width = 1000;

    // 创建一个1000x1000的随机图像
    std::vector<std::vector<double>> image(height, std::vector<double>(width));
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            image[i][j] = 1; // 使用分布生成随机值
        }
    }

    // 圆心和半径
    int x = 500;
    int y = 500;
    double r = 1;

    // 计算圆内的总和
    double result = sum_circle(image, x, y, r);
    std::cout << "圆内的总和: " << result << std::endl;

    return 0;
}