set(TEST_EXE ${PROJECT_NAME}_test)

add_executable(${TEST_EXE}
    main.cpp
    test_vec3.cpp
    test_mat3.cpp
    #test_approx.cpp
    test_box.cpp
    test_cuda_nblist.cpp
)

target_include_directories(${TEST_EXE} PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external/doctest-2.4.11
)
target_link_libraries(${TEST_EXE} PRIVATE ${PROJECT_NAME})

#set_target_properties(cuda_nblist_test PROPERTIES 
#    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
#)
