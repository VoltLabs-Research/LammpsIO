# Builds the reader core as a plain C++ static library, for consumers that link C++
# directly (CoreToolkit, and through it the analysis plugins).
#
# The Node addon is built by node-gyp from binding.gyp instead, and deliberately not from
# here: it needs Node's headers and produces a .node module. Both builds compile the same
# sources — everything except the N-API layer, which only node-gyp sees.

cmake_minimum_required(VERSION 3.16)
project(lammpsio LANGUAGES CXX VERSION 2.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(LAMMPSIO_BUILD_TESTS "Build the C++ smoke test" OFF)

add_library(lammpsio STATIC
    src/reader_registry.cpp
    src/formats/lammps_dump_text.cpp
    src/formats/lammps_dump_binary.cpp
    src/formats/lammps_dump_yaml.cpp
    src/formats/lammps_data.cpp
    src/formats/extxyz.cpp
)

# One include root for both trees, so `#include <lammpsio/...>` resolves the same way in
# the source checkout and against the installed package.
target_include_directories(lammpsio PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# The text reader parses in parallel with std::thread.
find_package(Threads REQUIRED)
target_link_libraries(lammpsio PUBLIC Threads::Threads)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # No -ffast-math: the hot path is byte scanning and fast_float, which it does not help,
    # and it makes behaviour around non-finite values unpredictable.
    target_compile_options(lammpsio PRIVATE -O3)
endif()

set_target_properties(lammpsio PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS lammpsio ARCHIVE DESTINATION lib)

install(DIRECTORY include/ DESTINATION include)

if(LAMMPSIO_BUILD_TESTS)
    add_executable(lammpsio_smoke test/cpp/smoke.cpp)
    target_link_libraries(lammpsio_smoke PRIVATE lammpsio)
endif()
