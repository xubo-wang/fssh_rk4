cmake_minimum_required(VERSION 3.5)
project(fssh_rk4)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    fssh_rk4/fssh_rk4.cc
    fssh_rk4/fssh_rk4.h
    fssh_rk4/main.cpp)

add_executable(fssh_rk4 ${SOURCE_FILES})

INCLUDE_DIRECTORIES("D:/eigen")