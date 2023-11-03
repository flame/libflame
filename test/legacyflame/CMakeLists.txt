###############################################################################
# Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################
add_executable(test_libFLAME "")

if(UNIX)
    set(Legacy_Flame_FLAGS "-lm -lpthread")
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

# Added Legacy Test to ctest
if(WIN32)
    set(LEGACY_FLAME_EXE_LOCATION ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})
else()
    set(LEGACY_FLAME_EXE_LOCATION ${CMAKE_CURRENT_SOURCE_DIR})
endif()

add_test(NAME legacy_flame_test
        COMMAND $<TARGET_FILE:test_libFLAME>
        WORKING_DIRECTORY ${LEGACY_FLAME_EXE_LOCATION})
