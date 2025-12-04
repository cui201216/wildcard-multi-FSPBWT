cmake_minimum_required(VERSION 3.28)
project(wild-multi-FSPBWT)

set(CMAKE_CXX_STANDARD 14)


# 设置静态链接标志
set(CMAKE_EXE_LINKER_FLAGS "-static")

# 添加可执行文件
add_executable(wild-multiFSPBWT main.cpp tools.h)
