add_executable(test_libFLAME_aocl "")
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
if(WIN32)
    set(LDFLAGS "")
    set_target_properties(test_libFLAME_aocl PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
elseif(UNIX)
    set(LDFLAGS "-lm -lpthread")
    set_target_properties(test_libFLAME_aocl PROPERTIES OUTPUT_NAME test_libflame_aocl_progress.x)
    set_target_properties(test_libFLAME_aocl PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

target_link_libraries(test_libFLAME_aocl debug  "${LDFLAGS}" "${CMAKE_PROJECT_NAME}" blas)
target_link_libraries(test_libFLAME_aocl optimized "${LDFLAGS}" "${CMAKE_PROJECT_NAME}" blas)

if(BLAS_HEADER_PATH)
    target_include_directories(test_libFLAME_aocl PRIVATE ${BLAS_HEADER_PATH})
else(BLAS_HEADER_PATH)
    message(STATUS "--> BLAS header file path not provided")
endif()

target_include_directories(test_libFLAME_aocl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_subdirectory(src)
