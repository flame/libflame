###############################################################################
# Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################
cmake_minimum_required(VERSION 3.0)

project(test_dgetrf)

if(WIN32)
    set(OpenMP_libomp_LIBRARY "C:/Program Files/LLVM/lib/libomp.lib" CACHE STRING "openmp library path")
    set(TEST_FLAGS "${OpenMP_libomp_LIBRARY}")
elseif(UNIX)
    set(OpenMP_libomp_LIBRARY "-fopenmp" CACHE STRING "openmp library path")
    set(Test_C_Flags "-g -O2 -Wall -Wno-comment" CACHE STRING "C Flags")
    set(Test_LD_Flags "-lm -lpthread" CACHE STRING "LD FLags")
    set(TEST_FLAGS "${OpenMP_libomp_LIBRARY} ${Test_C_Flags} ${Test_LD_Flags}")
endif()

set(TEST_BLAS_LIBRARY_PATH "${EXT_BLAS_LIBRARY_DEPENDENCY_PATH}" CACHE STRING "BLAS library path name")
set(TEST_LAPACK_LIBRARY_PATH "${EXT_LAPACK_LIBRARY_PATH}" CACHE STRING "LAPACK library path name")
set(TEST_BLAS_LIBNAME "${EXT_BLAS_LIBNAME}" CACHE STRING "BLAS library name")
set(TEST_LAPACK_LIBNAME "${EXT_LAPACK_LIBNAME}" CACHE STRING "lapack library name")
set(TEST_FLAME_HEADER_PATH "${EXT_FLAME_HEADER_PATH}" CACHE STRING "Flame header path")

include_directories("${EXT_FLAME_HEADER_PATH}")
link_directories("${EXT_BLAS_LIBRARY_DEPENDENCY_PATH}" "${EXT_LAPACK_LIBRARY_PATH}")

add_executable(${PROJECT_NAME} "")

if(WIN32)
    target_link_libraries(${PROJECT_NAME} debug "${TEST_FLAGS}" "${EXT_LAPACK_LIBNAME}" "${EXT_BLAS_LIBNAME}")
    target_link_libraries(${PROJECT_NAME} optimized "${TEST_FLAGS}" "${EXT_LAPACK_LIBNAME}" "${EXT_BLAS_LIBNAME}")
elseif(UNIX)
    target_link_libraries(${PROJECT_NAME} debug "${TEST_FLAGS}" "${EXT_LAPACK_LIBNAME}" "-lstdc++" "${EXT_BLAS_LIBNAME}")
    target_link_libraries(${PROJECT_NAME} optimized "${TEST_FLAGS}" "${EXT_LAPACK_LIBNAME}" "-lstdc++" "${EXT_BLAS_LIBNAME}")
endif()

target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/sample_dgetrf.c)