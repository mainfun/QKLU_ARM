#include <gtest/gtest.h>
#include <cstdio>

void printMessage() {
    printf("Hello, World!\n");
}

TEST(PrintMessageTest, OutputIsCorrect) {
    // 创建一个足够大的缓冲区
    char buffer[256];

    // 使用 fmemopen 创建一个内存流
    FILE* mem_stream = fmemopen(buffer, sizeof(buffer), "w+");
    ASSERT_NE(mem_stream, nullptr); // 确保流创建成功

    // 保存原始的 stdout
    FILE* original_stdout = stdout;

    // 重定向 stdout 到内存流
    stdout = mem_stream;

    // 调用要测试的函数
    printMessage();

    // 恢复原始 stdout
    fflush(stdout);
    stdout = original_stdout;

    // 读取输出内容
    fflush(mem_stream); // 刷新流
    fseek(mem_stream, 0, SEEK_SET); // 移动到流的开始
    fgets(buffer, sizeof(buffer), mem_stream); // 读取输出内容

    // 关闭内存流
    fclose(mem_stream);

    // 测试输出内容
    EXPECT_STREQ(buffer, "Hello, World!\n");
}