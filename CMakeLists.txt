cmake_minimum_required(VERSION 3.16)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
project(QKLU LANGUAGES C CXX)
#enable_language(CUDA)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "${CMAKE_SOURCE_DIR}/lib")

include_directories(src)
include_directories(include)

set(CMAKE_BUILD_TYPE Debug) # 设置为 Debug 模式
# 启用 AddressSanitizer
if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -g")
    set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")
endif ()
add_compile_options(-g -Wall -fopenmp -O3)

add_executable(QKLU main.c
        src/base/matrix.c
        src/base/plot.c
        src/base/file.c
        src/base/vector_util.cpp
        src/base/sort.c
        src/base/identify_chip.c
        src/base/forest.c
        src/base/malloc.c
        src/pivoting.c
        src/reorder.c
        src/symbolic_analysis.cpp
        src/preprocess.c
        #        src/numerical.c
        src/factor_v00.c
        src/check.c
        src/solving.c
        src/layer_matrix.c
        src/numerical.c
)
#--------------------------依赖库------------------------------
set(CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/jpeg-turbo/3.0.2/")
set(CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/suite-sparse/7.7.0/")
#set(CMAKE_PREFIX_PATH "/opt/intel/oneapi/mkl/latest/")
set(CMAKE_PREFIX_PATH "/opt/homebrew/opt/openblas/")

#JPEG
find_package(JPEG REQUIRED)
if (JPEG_FOUND)
    include_directories(${JPEG_INCLUDE_DIRS})
    target_link_libraries(QKLU ${JPEG_LIBRARIES})
else ()
    message(FATAL_ERROR "JPEG library not found")
endif ()

#AMD
find_library(AMD_LIBRARIES NAMES amd PATHS ${SuiteSparse_LIBRARIES})
if (AMD_LIBRARIES)
    include_directories("/usr/include/suitesparse/")
    include_directories("/opt/homebrew/Cellar/suite-sparse/7.7.0/include/suitesparse")
    target_link_libraries(QKLU ${AMD_LIBRARIES})
else ()
    message(WARNING "AMD library not found")
endif ()

#find_package(BLAS REQUIRED)
#find_package(MKL CONFIG REQUIRED)
#message(STATUS "${MKL_IMPORTED_TARGETS}")
#target_link_libraries(QKLU MKL::MKL)

# 查找Google Test库
#enable_testing()
#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS})

#OpenBLAS
include_directories(/opt/homebrew/opt/openblas/include)
link_directories(/opt/homebrew/opt/openblas/lib)
find_package(OpenBLAS REQUIRED)
#target_link_libraries(QKLU OpenBLAS)

#----复制库----

#--------------------------依赖库end------------------------------

# 自动发现测试
add_executable(gtest
        src/base/matrix.c
        src/base/plot.c
        src/base/file.c
        src/base/vector_util.cpp
        src/base/sort.c
        src/base/identify_chip.c
        src/base/forest.c
        src/base/malloc.c
        src/pivoting.c
        src/reorder.c
        src/symbolic_analysis.cpp
        src/preprocess.c
        #        src/numerical.c
        src/factor_v02.c
        src/check.c
#                        test/test_spgemm.c
        #        test/test_etree.cpp
        #        test/test_load_matrix.cpp
        #        test/test_mc64.cpp
        src/solving.c
        src/layer_matrix.c
#        test/test_gessm.c
        test/test_L2mtriax.c
        src/numerical.c
        #        test/test_lower_solving.c
        #        test/test_upper_solving.cpp
        #        test/mkl_lu_sp.cpp
        #                test/test_DTSTRF.c
        #        test/test_dense_lu.c
        #                test/test_lu_sp.c
        #        src/factor_cuda_v00.cu
)
target_link_libraries(gtest ${AMD_LIBRARIES})
target_link_libraries(gtest ${JPEG_LIBRARIES})
#target_link_libraries(gtest GTest::gtest GTest::gtest_main)
target_link_libraries(gtest OpenBLAS)
find_package(OpenMP REQUIRED)
target_link_libraries(QKLU OpenMP::OpenMP_CXX OpenMP::OpenMP_C)
target_link_libraries(gtest OpenMP::OpenMP_CXX OpenMP::OpenMP_C)
#target_link_libraries(gtest MKL::MKL)
#include(GoogleTest)
#gtest_discover_tests(gtest)

add_executable(test_remove_idx test/test_remove_idx.c)