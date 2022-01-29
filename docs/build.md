---
layout: default
title: Build
nav_order: 1
---

# Dependencies

In order to build and run Quill, the following dependencies are required:
* C++ compiler with C++11 support;
* CMake (version 3.10 or higher);
* MPI implementation;
* Python3 interpreter (*optional*: NumPy and Matplotlib for data analysis).

Quill is developed and tested on Linux with the use of *g++* and *clang++* compilers and *OpenMPI*.

# Build

To build Quill, run CMake in the `build` folder, e.g.
```
cd build
cmake ..
make
```
Or use your preferred way of building with CMake.

Building in any other folder is not recommended, as scripts rely on executables being present in the `build` folder.

Quill can be compiled without QED support by setting the CMake option QUILL_ENABLE_QED to OFF (ON by default),
```
cmake .. -DQUILL_ENABLE_QED=OFF
```