cmake_minimum_required(VERSION 3.11)

project(zero
    LANGUAGES CXX
)

add_library(zero INTERFACE)

target_sources(zero INTERFACE "functional_init_functions.hpp")
target_compile_options(zero INTERFACE -O3 -march=native -ffast-math -std=c++0x -Wall -Wextra)

