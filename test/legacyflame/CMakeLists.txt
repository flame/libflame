###############################################################################
# Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################

add_executable(test_libFLAME "")
target_link_libraries(test_libFLAME debug "${PROJECT_NAME}.lib" "${CMAKE_EXT_BLIS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLIS_LIBNAME}")
target_link_libraries(test_libFLAME optimized "${PROJECT_NAME}.lib")
target_include_directories(test_libFLAME PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
add_subdirectory(src)
