//
// Created by mainf on 2024/9/15.
//

#include "plot.h"

void csr2image(const SparseMatrix *matrix, const char *filename, int img_width, int img_height) {
    if (matrix == NULL)
        LOG_ERROR("matrix is NULL");
    if (img_width > matrix->num_row) {
        LOG_WARN("img width so big.");
    }
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *outfile = fopen(filename, "wb");
    if (!outfile) {
        LOG_ERROR("Can't open %s for writing\n", filename);
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);

    cinfo.image_width = img_width;
    cinfo.image_height = img_height;
    cinfo.input_components = 1;
    cinfo.in_color_space = JCS_GRAYSCALE;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 100, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    INDEX_TYPE block_width = matrix->num_col / img_width;
    INDEX_TYPE block_height = matrix->num_row / img_height;

    LOG_DEBUG("block_size is %lld x %lld", block_width, block_height);

    unsigned char *image_buffer = calloc(img_width * img_height, sizeof(unsigned char));
    if (!image_buffer) {
        LOG_ERROR("Memory allocation failed");
    }

    // Calculate pixel values
    for (int r = 0; r < img_height; r++) {
        for (int c = 0; c < img_width; c++) {
            INDEX_TYPE pixel_index = r * img_width + c;
            INDEX_TYPE start_row = r * block_height;
            INDEX_TYPE end_row = (r + 1) * block_height;
            INDEX_TYPE start_col = c * block_width;
            INDEX_TYPE end_col = (c + 1) * block_width;
            INDEX_TYPE non_zero_count = 0;

            for (INDEX_TYPE row = start_row; row < end_row; row++) {
                for (INDEX_TYPE idx = matrix->row_pointers[row]; idx < matrix->row_pointers[row + 1]; ++idx) {
                    if (matrix->col_indices[idx] >= start_col && matrix->col_indices[idx] < end_col) {
                        non_zero_count++;
                    }
                }
            }
            double density = (double) non_zero_count / (double) (block_width * block_height);
            double adjusted_density = density > 0 ? 1 : 0;
            image_buffer[pixel_index] = (unsigned char) (255 * (1 - adjusted_density)); // More non-zeros -> darker
        }
    }

    JSAMPROW row_pointer[1];
    while (cinfo.next_scanline < cinfo.image_height) {
        row_pointer[0] = &image_buffer[cinfo.next_scanline * img_width];
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    free(image_buffer);
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
}

void csr2_RGB_image(const SparseMatrix *matrix, const char *filename, int img_width, int img_height) {
    if (matrix == NULL)
        LOG_ERROR("matrix is NULL");
    if (img_width > matrix->num_row) {
        LOG_WARN("img width so big.");
    }
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *outfile = fopen(filename, "wb");
    if (!outfile) {
        LOG_ERROR("Can't open %s for writing\n", filename);
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);

    cinfo.image_width = img_width;
    cinfo.image_height = img_height;
    cinfo.input_components = 3; // 设置为3，表示RGB图像
    cinfo.in_color_space = JCS_RGB; // 使用RGB颜色模式

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 100, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    INDEX_TYPE block_width = matrix->num_col / img_width;
    INDEX_TYPE block_height = matrix->num_row / img_height;

    LOG_DEBUG("block_size is %lld x %lld", block_width, block_height);

    unsigned char *image_buffer = calloc(img_width * img_height * 3, sizeof(unsigned char)); // 3表示RGB三个分量
    if (!image_buffer) {
        LOG_ERROR("Memory allocation failed");
    }

    // Calculate pixel values
    for (int r = 0; r < img_height; r++) {
        for (int c = 0; c < img_width; c++) {
            INDEX_TYPE pixel_index = (r * img_width + c) * 3; // 每个像素占3个字节
            INDEX_TYPE start_row = r * block_height;
            INDEX_TYPE end_row = (r + 1) * block_height;
            INDEX_TYPE start_col = c * block_width;
            INDEX_TYPE end_col = (c + 1) * block_width;
            INDEX_TYPE non_zero_count = 0;

            for (INDEX_TYPE row = start_row; row < end_row; row++) {
                for (INDEX_TYPE idx = matrix->row_pointers[row]; idx < matrix->row_pointers[row + 1]; ++idx) {
                    if (matrix->col_indices[idx] >= start_col && matrix->col_indices[idx] < end_col) {
                        non_zero_count++;
                    }
                }
            }

            double density = (double) non_zero_count / (double) (block_width * block_height);
            // density 的值在[0, 1]之间
            // 当密度较高时，使用红色。较低时，使用白色（或者其他颜色）
            unsigned char red = 255; // 红色分量随着密度增加
            unsigned char green = 255; // 始终保持绿色为0
            unsigned char blue = 255; // 始终保持蓝色为0
            if (density > 0.2) {
                red = 255;
                green = 0;
                blue = 0;
            } else if (density > 0) {
                red = 0;
                green = 0;
                blue = 0;
            }

            // 设置当前像素的 RGB 值
            image_buffer[pixel_index] = red;
            image_buffer[pixel_index + 1] = green;
            image_buffer[pixel_index + 2] = blue;
        }
    }

    JSAMPROW row_pointer[1];
    while (cinfo.next_scanline < cinfo.image_height) {
        row_pointer[0] = &image_buffer[cinfo.next_scanline * img_width * 3]; // 3个字节表示RGB
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    free(image_buffer);
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
}

///带坐标系
void csr2_RGB_image_coor(const SparseMatrix *matrix, const char *filename, int img_width, int img_height) {
    if (matrix == NULL)
        LOG_ERROR("matrix is NULL");
    if (img_width > matrix->num_row) {
        LOG_WARN("img width so big.");
    }
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *outfile = fopen(filename, "wb");
    if (!outfile) {
        LOG_ERROR("Can't open %s for writing\n", filename);
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);

    cinfo.image_width = img_width;
    cinfo.image_height = img_height;
    cinfo.input_components = 3; // 设置为3，表示RGB图像
    cinfo.in_color_space = JCS_RGB; // 使用RGB颜色模式

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 100, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    INDEX_TYPE block_width = matrix->num_col / img_width;
    INDEX_TYPE block_height = matrix->num_row / img_height;

    LOG_DEBUG("block_size is %lld x %lld", block_width, block_height);

    unsigned char *image_buffer = calloc(img_width * img_height * 3, sizeof(unsigned char)); // 3表示RGB三个分量
    if (!image_buffer) {
        LOG_ERROR("Memory allocation failed");
    }

    // 计算矩阵的像素值
    for (int r = 0; r < img_height; r++) {
        for (int c = 0; c < img_width; c++) {
            INDEX_TYPE pixel_index = (r * img_width + c) * 3; // 每个像素占3个字节
            INDEX_TYPE start_row = r * block_height;
            INDEX_TYPE end_row = (r + 1) * block_height;
            INDEX_TYPE start_col = c * block_width;
            INDEX_TYPE end_col = (c + 1) * block_width;
            INDEX_TYPE non_zero_count = 0;

            for (INDEX_TYPE row = start_row; row < end_row; row++) {
                for (INDEX_TYPE idx = matrix->row_pointers[row]; idx < matrix->row_pointers[row + 1]; ++idx) {
                    if (matrix->col_indices[idx] >= start_col && matrix->col_indices[idx] < end_col) {
                        non_zero_count++;
                    }
                }
            }

            double density = (double) non_zero_count / (double) (block_width * block_height);
            // density 的值在[0, 1]之间
            // 当密度较高时，使用红色。较低时，使用白色（或者其他颜色）
            unsigned char red = 255; // 红色分量随着密度增加
            unsigned char green = 255; // 始终保持绿色为0
            unsigned char blue = 255; // 始终保持蓝色为0
            if (density > 0.2) {
                red = 255;
                green = 0;
                blue = 0;
            } else if (density > 0) {
                red = 0;
                green = 0;
                blue = 0;
            }

            // 设置当前像素的 RGB 值
            image_buffer[pixel_index] = red;
            image_buffer[pixel_index + 1] = green;
            image_buffer[pixel_index + 2] = blue;
        }
    }

    // 1. 绘制黑色边框：四周的黑色边框
    for (int r = 0; r < img_height; r++) {
        for (int c = 0; c < img_width; c++) {
            // 边框在四周，设置为黑色
            if (r == 0 || r == img_height - 1 || c == 0 || c == img_width - 1) {
                int pixel_index = (r * img_width + c) * 3;
                image_buffer[pixel_index] = 0; // 红色
                image_buffer[pixel_index + 1] = 0; // 绿色
                image_buffer[pixel_index + 2] = 0; // 蓝色
            }
        }
    }

    // 2. 绘制坐标系（x轴和y轴）：并添加刻度
    int axis_thickness = 3; // 坐标轴线条的粗细
    int tick_size = 20; // 刻度线的长度
    double max_tick_value = 10.0 / (matrix->num_row); // 刻度的范围

    // 2.1 绘制横坐标轴 (x轴)，从左到右
    for (int c = 1; c < img_width - 1; c++) {
        for (int i = 0; i < axis_thickness; i++) { // 增加厚度
            int pixel_index = ((0 + i) * img_width + c) * 3; // x轴位于第一行和上下各1行
            image_buffer[pixel_index] = 0; // 红色分量
            image_buffer[pixel_index + 1] = 255; // 绿色分量
            image_buffer[pixel_index + 2] = 0; // 蓝色分量（红色）
        }
    }

    // 2.2 绘制纵坐标轴 (y轴)，从上到下
    for (int r = 1; r < img_height - 1; r++) {
        for (int i = 0; i < axis_thickness; i++) { // 增加厚度
            int pixel_index = (r * img_width + i) * 3; // y轴位于第一列和左右各1列
            image_buffer[pixel_index] = 0; // 红色分量
            image_buffer[pixel_index + 1] = 255; // 绿色分量
            image_buffer[pixel_index + 2] = 0; // 蓝色分量（红色）
        }
    }

    // 2.3 绘制刻度线：在x轴和y轴外侧绘制10个刻度
    for (int i = 1; i <= 10; i++) {
        // x轴刻度线：根据比例绘制
        int x_tick_pos = (i * (img_width - 2)) / 10 + 1; // 刻度位置，1为左边框开始，除去边框
        for (int t = 0; t < tick_size; t++) { // 加粗刻度线
            int pixel_index = ((0 + t) * img_width + x_tick_pos) * 3; // 竖着绘制刻度
            image_buffer[pixel_index] = 0; // 红色分量
            image_buffer[pixel_index + 1] = 0; // 绿色分量
            image_buffer[pixel_index + 2] = 0; // 蓝色分量
        }

        // y轴刻度线
        int y_tick_pos = (i * (img_height - 2)) / 10 + 1; // 刻度位置，1为上边框开始，除去边框
        for (int t = 0; t < tick_size; t++) {
            int pixel_index = (y_tick_pos * img_width + (0 + t)) * 3; // 横着绘制刻度
            image_buffer[pixel_index] = 0; // 红色分量
            image_buffer[pixel_index + 1] = 0; // 绿色分量
            image_buffer[pixel_index + 2] = 0; // 蓝色分量
        }
    }

    JSAMPROW row_pointer[1];
    while (cinfo.next_scanline < cinfo.image_height) {
        row_pointer[0] = &image_buffer[cinfo.next_scanline * img_width * 3]; // 3个字节表示RGB
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    free(image_buffer);
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
}


void csr2image_block(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai,
                     const char *filename,
                     INDEX_TYPE row_start, INDEX_TYPE row_end,
                     INDEX_TYPE col_start, INDEX_TYPE col_end) {
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *outfile = fopen(filename, "wb");
    if (!outfile) {
        LOG_ERROR("Can't open %s for writing\n", filename);
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);

    long img_width = row_end - row_start;
    long img_height = col_end - col_start;
    cinfo.image_width = img_width;
    cinfo.image_height = img_height;
    cinfo.input_components = 1;
    cinfo.in_color_space = JCS_GRAYSCALE;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 100, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    unsigned char *image_buffer = calloc(img_width * img_height, sizeof(unsigned char));
    if (!image_buffer) {
        LOG_ERROR("Memory allocation failed");
    }

    // Calculate pixel values
    // int pixel_index = 0;
    // int found;
    // for (INDEX_TYPE i = row_start; i < row_end; i++) {
    //     for (INDEX_TYPE j = col_start; j < col_end; j++) {
    //         found = 0;
    //         for (INDEX_TYPE p = Ap[i]; p < Ap[i + 1]; p++) {
    //             INDEX_TYPE col_index = Ai[p];
    //             if (col_index == j) {
    //                 found = 1;
    //                 image_buffer[pixel_index++] = (unsigned char) 0;
    //                 break;
    //             }
    //         }
    //         if (!found) {
    //             image_buffer[pixel_index++] = (unsigned char) 255;
    //         }
    //     }
    // }
    memset(image_buffer, 255, sizeof(unsigned char) * img_width * img_height);
    for (INDEX_TYPE i = row_start; i < row_end; i++) {
        for (INDEX_TYPE p = Ap[i]; p < Ap[i + 1]; p++) {
            INDEX_TYPE col_index = Ai[p];
            if (col_start <= col_index && col_index < col_end) {
                image_buffer[(i-row_start) * img_width + col_index-col_start] = (unsigned char) 0;
            }
        }
    }

    JSAMPROW row_pointer[1];
    while (cinfo.next_scanline < cinfo.image_height) {
        row_pointer[0] = &image_buffer[cinfo.next_scanline * img_width];
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    free(image_buffer);
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
}

void csr2image_block_cut(const INDEX_TYPE *Ap, const INDEX_TYPE *Ai,
                         const char *filename,
                         INDEX_TYPE row_start, INDEX_TYPE row_end,
                         INDEX_TYPE col_start, INDEX_TYPE col_end,
                         const INDEX_TYPE *cut_points, INDEX_TYPE cut_points_count) {
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *outfile = fopen(filename, "wb");
    if (!outfile) {
        fprintf(stderr, "Can't open %s for writing\n", filename);
        return;
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);

    INDEX_TYPE img_width = row_end - row_start;
    INDEX_TYPE img_height = col_end - col_start;
    cinfo.image_width = img_width;
    cinfo.image_height = img_height;
    cinfo.input_components = 3; // RGB
    cinfo.in_color_space = JCS_RGB;

    // Set the compression mode to "lossless"
    cinfo.raw_data_in = TRUE; // Enable raw data for lossless compression
    cinfo.optimize_coding = TRUE; // Optimize coding for lossless compression

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, 100, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    // Allocate memory for the RGB image
    unsigned char *image_buffer = calloc(img_width * img_height * 3, sizeof(unsigned char));
    if (!image_buffer) {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }

    // Initialize the image with white pixels
    for (int i = 0; i < img_width * img_height; i++) {
        image_buffer[i * 3 + 0] = 255; // Red channel
        image_buffer[i * 3 + 1] = 255; // Green channel
        image_buffer[i * 3 + 2] = 255; // Blue channel
    }

    // Calculate pixel values and add red lines for cut_points
    INDEX_TYPE pixel_index = 0;
    int found;
    for (INDEX_TYPE i = row_start; i < row_end; i++) {
        for (INDEX_TYPE j = col_start; j < col_end; j++) {
            found = 0;
            for (INDEX_TYPE p = Ap[i]; p < Ap[i + 1]; p++) {
                INDEX_TYPE col_index = Ai[p];
                if (col_index == j) {
                    found = 1;
                    image_buffer[pixel_index * 3 + 0] = 0; // Set red to 0
                    image_buffer[pixel_index * 3 + 1] = 0; // Set green to 0
                    image_buffer[pixel_index * 3 + 2] = 0; // Set blue to 0
                    break;
                }
            }
            if (!found) {
                image_buffer[pixel_index * 3 + 0] = 255; // Red channel
                image_buffer[pixel_index * 3 + 1] = 255; // Green channel
                image_buffer[pixel_index * 3 + 2] = 255; // Blue channel
            }
            pixel_index++;
        }
    }

    // Draw red lines at cut_points
    for (int i = 0; i < cut_points_count; i++) {
        INDEX_TYPE cut_point = cut_points[i];

        // Draw red line horizontally (at cut_point-th row)
        if (cut_point >= row_start && cut_point < row_end) {
            for (INDEX_TYPE j = col_start; j < col_end; j++) {
                INDEX_TYPE pixel_idx = (cut_point - row_start) * img_width + (j - col_start);
                image_buffer[pixel_idx * 3 + 0] = 255; // Red
                image_buffer[pixel_idx * 3 + 1] = 0; // Green
                image_buffer[pixel_idx * 3 + 2] = 0; // Blue
            }
        }

        // Draw red line vertically (at cut_point-th column)
        if (cut_point >= col_start && cut_point < col_end) {
            for (INDEX_TYPE i = row_start; i < row_end; i++) {
                INDEX_TYPE pixel_idx = (i - row_start) * img_width + (cut_point - col_start);
                image_buffer[pixel_idx * 3 + 0] = 0; // Red
                image_buffer[pixel_idx * 3 + 1] = 0; // Green
                image_buffer[pixel_idx * 3 + 2] = 255; // Blue
            }
        }
    }

    // Write image to the file
    JSAMPROW row_pointer[1];
    while (cinfo.next_scanline < cinfo.image_height) {
        row_pointer[0] = &image_buffer[cinfo.next_scanline * img_width * 3];
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    // Clean up
    free(image_buffer);
    jpeg_finish_compress(&cinfo);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
}
