cmake_minimum_required(VERSION 3.11)

project(quill)

# quill deps #

find_package(MPI REQUIRED)

include(FetchContent)
FetchContent_Declare(
    toml11
    GIT_REPOSITORY "https://github.com/ToruNiina/toml11.git"
)
FetchContent_MakeAvailable(toml11)

set(DIR_FUNCTIONAL_INIT ${CMAKE_CURRENT_SOURCE_DIR}/functional_init/zero/)
include(${DIR_FUNCTIONAL_INIT}/CMakeLists.txt)
include_directories(${DIR_FUNCTIONAL_INIT})

# quill #

set(SOURCE_QUILL
    advance.cpp deposition.cpp fout.cpp init.cpp low-level.cpp main.cpp
    mc.cpp plist.cpp pusher.cpp maxwell.cpp containers.cpp balancing.cpp
)
set(FLAGS_QUILL
    -MMD -MP -MF -O3 -march=native -ffast-math -std=c++0x -Wall -Wextra
)

add_executable(quill ${SOURCE_QUILL})
target_include_directories(quill PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_compile_options(quill PRIVATE ${MPI_CXX_COMPILE_FLAGS} ${FLAGS_QUILL})
target_link_libraries(quill toml11 ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})

# chameleon #

set(SOURCE_CHAMELEON "chameleon/chameleon-core.cpp")
set(FLAGS_CHAMELEON -O3 -std=c++11)

add_library(chameleon SHARED ${SOURCE_CHAMELEON})
target_compile_options(chameleon PRIVATE ${FLAGS_CHAMELEON})
set_target_properties(chameleon PROPERTIES LIBRARY_OUTPUT_DIRECTORY "chameleon/")
