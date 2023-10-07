if(BUILD_SHARED_LIBS)
    add_library(dpnblist SHARED)
else()
    add_library(dpnblist STATIC)
endif()

target_include_directories(dpnblist PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)
target_sources(dpnblist PRIVATE
    refCellList.cpp
    box.cpp
)
