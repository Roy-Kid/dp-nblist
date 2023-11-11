if(BUILD_SHARED_LIBS)
    add_library(dpnblist SHARED)
else()
    add_library(dpnblist STATIC)
endif()

find_package(pybind11 REQUIRED)
target_include_directories(dpnblist PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${pybind11_INCLUDE_DIRS}
)
target_sources(dpnblist PRIVATE
    refCellList.cpp
    box.cpp
)
