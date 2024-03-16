cmake_minimum_required(VERSION 3.27)
project(M_DCB_C)

# 设置OpenBLAS目录
set(OpenBLAS_DIR "lib/OpenBLAS")
set(SuiteSparse_DIR "/usr/local")

# 包含OpenBLAS头文件
include_directories(${OpenBLAS_DIR}/include)
include_directories(${SuiteSparse_DIR}/include/suitesparse)

# 链接OpenBLAS库
link_directories(${OpenBLAS_DIR}/lib)
link_directories(${SuiteSparse_DIR}/lib)

set(CMAKE_CXX_STANDARD 17)

add_executable(M_DCB_C main.cpp
)

target_link_libraries(M_DCB_C PRIVATE
        spqr
        ccolamd
        cholmod
        suitesparseconfig
        amd
        camd
        openblas
)







