cmake_minimum_required(VERSION 3.12)

# 定义项目名称
project(cuda_nblist LANGUAGES CUDA CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(CUDA_ARCHITECTURES_DEFAULT "61")
set(CUDA_ARCHITECTURES ${CUDA_ARCHITECTURES_DEFAULT} CACHE STRING "CUDA architectures")
# compile pylib
set(Pylib ${PROJECT_NAME}_py)

add_library(${Pylib} SHARED)

find_package(Python3 COMPONENTS Development)
target_include_directories(${Pylib} PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external/pybind11-2.11.1/include
    ${Python3_INCLUDE_DIRS}
)
target_sources(${Pylib} PRIVATE
    src/box.cpp
    src/cuda_nblist_py.cpp
)
target_sources(${Pylib} PRIVATE
    src/cuda_nblist.cu
)
set(POSITION_INDEPENDENT_CODE ON)
set_target_properties(${Pylib} PROPERTIES
    PREFIX ""
    CUDA_ARCHITECTURES ${CUDA_ARCHITECTURES}
)
target_compile_options(${Pylib} PRIVATE -g -O3)

# compile lib
add_subdirectory(src)

# compile executable
add_subdirectory(tests)
