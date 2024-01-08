
# compile lib
add_library(${PROJECT_NAME} SHARED)
target_sources(${PROJECT_NAME} PRIVATE
    box.cpp
)
target_sources(${PROJECT_NAME} PRIVATE
    cuda_nblist.cu
)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    #${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    PREFIX ""
    CUDA_ARCHITECTURES ${CUDA_ARCHITECTURES}
)
target_compile_options(${PROJECT_NAME} PRIVATE -g -O3)
