cmake_minimum_required(VERSION 3.28)
project(wild-multi-FSPBWT)

set(CMAKE_CXX_STANDARD 17)


# 设置静态链接标志
set(CMAKE_EXE_LINKER_FLAGS "-static")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# 添加可执行文件
add_executable(wild-multiFSPBWT main.cpp tools.h
)
# 2. 测试程序 (只包含 main_test.cpp)
# 这里的名字不能和上面一样，我们叫它 wild-multiFSPBWT_test
add_executable(wild-multiFSPBWT_test main_test.cpp)
