if((NOT EXT_LAPACK_LIBRARY_PATH) OR (NOT EXT_LAPACK_LIBNAME))
		set(EXT_LAPACK_LIBRARY_PATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
        set(EXT_LAPACK_LIBNAME "libflame.a")
endif()
enable_testing()
if(CMAKE_C_COMPILER MATCHES "clang$")
    # Run the netlib test using AOCC if clang compiler is detected 
    set(NETLIB_BASH_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/netlib-test/run-netlib-test-aocc.sh)
else()
    # Run the netlib test using GCC if gcc compiler is detected 
    set(NETLIB_BASH_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/netlib-test/run-netlib-test.sh)
endif()

add_test(netlib-test bash ${NETLIB_BASH_SCRIPT} BLAS_LIB_PATH=${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH} BLAS_LIB=${EXT_BLAS_LIBNAME}
LAPACK_LIB_PATH=${EXT_LAPACK_LIBRARY_PATH} LAPACK_LIB=${EXT_LAPACK_LIBNAME} ILP64=${ENABLE_ILP64} GCOV=${ENABLE_GCOV})

set_tests_properties(netlib-test PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/netlib-test)