cmake_minimum_required (VERSION 3.10)

project (Quill C CXX)

if (CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
    message(STATUS "No CMAKE_BUILD_TYPE specified, defaulting to ${CMAKE_BUILD_TYPE}")
else()
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")

find_package(MPI REQUIRED COMPONENTS C CXX)

add_executable(quill src/advance.cpp src/balancing.cpp src/containers.cpp src/deposition.cpp src/fout.cpp src/init.cpp 
               src/low-level.cpp src/main.cpp src/maxwell.cpp src/mc.cpp src/plist.cpp src/pusher.cpp
               src/parameters.cpp)
target_include_directories(quill PRIVATE "include")
target_link_libraries(quill PUBLIC MPI::MPI_CXX)
target_compile_options(quill PRIVATE -Wall -Wextra)

option(QUILL_ENABLE_QED "Compile Quill with QED support" ON)

if (QUILL_ENABLE_QED)
    message(STATUS "QUILL_ENABLE_QED is set to ON, Quill will be compiled with QED support")
else()
    message(STATUS "QUILL_ENABLE_QED is set to OFF, Quill will be compiled without QED support")
    target_compile_definitions(quill PRIVATE QUILL_NOQED)
endif()

add_library(chameleon SHARED lib/chameleon/chameleon-core.cpp)