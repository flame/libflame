include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(test_libFLAME "")
target_link_libraries(test_libFLAME debug "${PROJECT_NAME}.lib" "${CMAKE_EXT_BLIS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLIS_LIBNAME}")
target_link_libraries(test_libFLAME optimized "${PROJECT_NAME}.lib")
add_subdirectory(src)