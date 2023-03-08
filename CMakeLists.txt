cmake_minimum_required(VERSION 3.24)
project(lab3_math C)

set(CMAKE_C_STANDARD 17)

add_executable(lab3_math main.c rkf45.c cmathmsg.c)
target_link_libraries(lab3_math m)
