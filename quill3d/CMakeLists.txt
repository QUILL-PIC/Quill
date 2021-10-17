cmake_minimum_required(VERSION 3.11)

project(quill)

find_package(MPI REQUIRED)

include(FetchContent)
FetchContent_Declare(
    toml11
    GIT_REPOSITORY "https://github.com/ToruNiina/toml11.git"
)
FetchContent_MakeAvailable(toml11)

set(SOURCE_EXE
    advance.cpp deposition.cpp fout.cpp init.cpp low-level.cpp main.cpp
    mc.cpp plist.cpp pusher.cpp maxwell.cpp containers.cpp balancing.cpp
)

set(FLAGS_EXE -MMD -MP -MF -O3 -march=native -ffast-math -std=c++0x -Wall -Wextra)

add_executable(quill ${SOURCE_EXE})

target_include_directories(quill PRIVATE ${MPI_CXX_INCLUDE_PATH})

target_compile_options(quill PRIVATE ${MPI_CXX_COMPILE_FLAGS} ${FLAGS_EXE})

target_link_libraries(quill toml11 ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
