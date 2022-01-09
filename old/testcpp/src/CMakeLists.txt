##Copyright (C) 2020, Advanced Micro Devices, Inc.##

file(GLOB TESTCPP_SRC_FILES *.cc)

list(REMOVE_ITEM TESTCPP_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/libflame_test.cc)

foreach(CPP_SRCFILE ${TESTCPP_SRC_FILES})

    get_filename_component(CPP_SRCFILE_NODIR ${CPP_SRCFILE} NAME)
    string(REGEX MATCH "^(.*)\\.[^.]*$" INVAR ${CPP_SRCFILE_NODIR})

    add_executable(${CMAKE_MATCH_1} ${CPP_SRCFILE_NODIR})

    target_include_directories(${CMAKE_MATCH_1} PRIVATE
                        ${CMAKE_SOURCE_DIR}/src/aocl_dtl
                        ${CMAKE_SOURCE_DIR}/src/base/flamec/include
                        ${CMAKE_SOURCE_DIR}/src/testcpp/src
                        ${CMAKE_SOURCE_DIR}/src/src_cpp
                        ${CMAKE_SOURCE_DIR}/src/lapacke/LAPACKE/include
                        ${CMAKE_SOURCE_DIR}/src/base/flamec/blis/include
                        ${CMAKE_SOURCE_DIR}/src/base/flamec/control
                        ${CMAKE_SOURCE_DIR}/src/base/flamec/include
                        ${CMAKE_CURRENT_SOURCE_DIR}
                        )

    target_link_libraries(${CMAKE_MATCH_1} debug lapacke "${PROJECT_NAME}")
    target_link_libraries(${CMAKE_MATCH_1} optimized lapacke "${PROJECT_NAME}")
    target_compile_definitions(${CMAKE_MATCH_1} PRIVATE -DHAVE_CONFIG_H -DBLIS1_FROM_LIBFLAME)

    if (ENABLE_MULTITHREADING)
        if (ENABLE_OPENMP)
            target_compile_definitions(${CMAKE_MATCH_1} PRIVATE -DFLA_MULTITHREADING_MODEL=1)
        endif ()
        if (ENABLE_PTHREADS)
            target_compile_definitions(${CMAKE_MATCH_1} PRIVATE -DFLA_MULTITHREADING_MODEL=2)
        endif ()
    endif ()

endforeach ()





