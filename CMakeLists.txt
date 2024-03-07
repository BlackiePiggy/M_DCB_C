cmake_minimum_required(VERSION 3.27)
project(M_DCB_C)

include_directories("D:\\projects\\M_DCB_C\\lib\\eigen-3.4.0")

set(CMAKE_CXX_STANDARD 17)

add_executable(M_DCB_C main.cpp
        ReadFile.cpp
        ReadFile.h)
