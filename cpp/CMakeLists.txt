cmake_minimum_required(VERSION 3.15)

project(dpnblist LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    # compile target as shared lib
    set(POSITION_INDEPENDENT_CODE ON)
    option(BUILD_SHARED_LIBS "Build shared libraries" ON)
    add_subdirectory(src)
    target_compile_options(dpnblist PRIVATE -g)
    execute_process(COMMAND python3-config --extension-suffix PYTHON_MODULE_EXTENSION)
    set_target_properties(dpnblist PROPERTIES 
        PREFIX ""
        SUFFIX ${PYTHON_MODULE_EXTENSION}
    )
    #add_subdirectory(tests)
    # enable_testing()
endif()