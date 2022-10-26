###############################################################################
# Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################
cmake_minimum_required(VERSION 3.0)
project(test_libFLAME_main)

option (ENABLE_ILP64 "Enable ILP64" OFF)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

set(OpenMP_libomp_LIBRARY "C:/Program Files/LLVM/lib/libomp.lib" CACHE STRING "openmp library path")
set(TEST_BLIS_LIBRARY_PATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" CACHE STRING "AOCL-BLIS library path name")
set(TEST_LIBFLAME_LIBRARY_PATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" CACHE STRING "AOCL-LIBFLAME library path name")
set(TEST_BLIS_LIBNAME "AOCL-LibBlis-Win-MT-dll.lib" CACHE STRING "AOCL-BLIS library name")
set(TEST_LIBFLAME_LIBNAME "AOCL-LibFlame-Win-MT-dll.lib" CACHE STRING "AOCL-LIBFLAME library name")

add_definitions(-D_CRT_SECURE_NO_WARNINGS)

if(ENABLE_ILP64)
add_definitions(-DFLA_ENABLE_ILP64)
endif()

link_directories("${TEST_BLIS_LIBRARY_PATH}" "${TEST_LIBFLAME_LIBRARY_PATH}")

add_executable(test_libFLAME_main "")
message("${CMAKE_BUILD_TYPE}")
target_link_libraries(test_libFLAME_main debug "${OpenMP_libomp_LIBRARY}" "${TEST_BLIS_LIBNAME}" "${TEST_LIBFLAME_LIBNAME}")
target_link_libraries(test_libFLAME_main optimized "${OpenMP_libomp_LIBRARY}" "${TEST_BLIS_LIBNAME}" "${TEST_LIBFLAME_LIBNAME}")

target_include_directories(test_libFLAME_main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(test_libFLAME_main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/validate_src)
add_subdirectory(src)
add_subdirectory(validate_src)