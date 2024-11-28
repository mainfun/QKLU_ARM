//
// Created by mainf on 2024/9/15.
//

#include "plot.h"

void csr2image(const SparseMatrix *matrix, const char *filename, int img_width, int img_height) {
    if (matrix == NULL)LOG_ERROR("matrix is NULL");
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
            image_buffer[pixel_index] = (unsigned char) (255 * (1 - adjusted_density));  // More non-zeros -> darker
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
    if (matrix == NULL) LOG_ERROR("matrix is NULL");
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
    cinfo.input_components = 3;  // 设置为3，表示RGB图像
    cinfo.in_color_space = JCS_RGB;  // 使用RGB颜色模式

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
            INDEX_TYPE pixel_index = (r * img_width + c) * 3;  // 每个像素占3个字节
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
            unsigned char green = 255;   // 始终保持绿色为0
            unsigned char blue = 255;    // 始终保持蓝色为0
            if(density>0.2) {
                red = 255;
                green = 0;
                blue = 0;
            }else if(density>0){
                red=0;
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

    int img_width = row_end - row_start;
    int img_height = col_end - col_start;
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
    int pixel_index = 0;
    int found;
    for (INDEX_TYPE i = row_start; i < row_end; i++) {
        for (INDEX_TYPE j = col_start; j < col_end; j++) {
            found = 0;
            for (INDEX_TYPE p = Ap[i]; p < Ap[i + 1]; p++) {
                INDEX_TYPE col_index = Ai[p];
                if (col_index == j) {
                    found = 1;
                    image_buffer[pixel_index++] = (unsigned char) 0;
                    break;
                }
            }
            if (!found) {
                image_buffer[pixel_index++] = (unsigned char) 255;
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