cmake_minimum_required(VERSION 3.27)
project(M_DCB_C)

# 设置OpenBLAS目录
set(OpenBLAS_DIR "lib\\OpenBLAS-0.3.26-x64")

# 包含OpenBLAS头文件
include_directories(${OpenBLAS_DIR}/include)

# 链接OpenBLAS库
link_directories(${OpenBLAS_DIR}/lib)

#include_directories("D:\\projects\\M_DCB_C\\lib\\eigen-3.4.0")

set(CMAKE_CXX_STANDARD 17)

add_executable(M_DCB_C main.cpp
        ReadFile.cpp
        ReadFile.h)

target_link_libraries(M_DCB_C libopenblas.a)