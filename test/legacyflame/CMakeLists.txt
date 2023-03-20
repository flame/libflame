###############################################################################
# Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################
add_executable(test_libFLAME "")

if(UNIX)
    set(Legacy_Flame_FLAGS "-lm")
endif()
add_definitions(-DBLIS1_FROM_LIBFLAME)
target_link_libraries(test_libFLAME debug "${CMAKE_PROJECT_NAME}" "${Legacy_Flame_FLAGS}" blas)

target_link_libraries(test_libFLAME optimized "${Legacy_Flame_FLAGS}" "${CMAKE_PROJECT_NAME}" blas)
target_include_directories(test_libFLAME PRIVATE "${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}" blas)
target_include_directories(test_libFLAME PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
# to move executable to the legacy flame folder LINUX
if(UNIX)
    set_target_properties(test_libFLAME PROPERTIES OUTPUT_NAME test_libflame.x)
    set_target_properties(test_libFLAME PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

add_subdirectory(src)
