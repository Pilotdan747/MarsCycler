cmake_minimum_required(VERSION 3.14)
project(C__)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER /usr/local/opt/llvm/bin/clang++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

add_executable(MarsCycler main.cpp helperFuncs.cpp helperFuncs.h)
