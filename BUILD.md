# STEPS TO COMPILE LIBFLAME USING CMAKE
## 1. Generating Libflame library  
    Create a build directory e.g. build1 
        mkdir build1;
        cd build1;

    use the following command to configure project
        GCC:
        -----------
        With LP64
            cmake ../ -DENABLE_AMD_FLAGS=ON -DCMAKE_INSTALL_PREFIX=<path>
        With ILP64
            cmake ../ -DENABLE_ILP64=ON -DENABLE_AMD_FLAGS=ON -DCMAKE_INSTALL_PREFIX=<path>
        AOCC:
        -----------
        export CC=clang
        export CXX=clang++
        export FC=flang
        export FLIBS="-lflang"

        With LP64
            cmake ../ -DENABLE_AMD_AOCC_FLAGS=ON -DCMAKE_INSTALL_PREFIX=<path>
        With ILP64
            cmake ../ -DENABLE_ILP64=ON -DENABLE_AMD_AOCC_FLAGS=ON -DCMAKE_INSTALL_PREFIX=<path>
    
    Shared library is turned on by default. To generate Static library provide additional option
        -DBUILD_SHARED_LIBS=OFF

    compile library using following command. This will generate libflame.a/libflame.so library in the bin directory
        cmake --build . -j OR make -j

## 2. Building Test
    In order to build tests an an additional flag can be set to ON
        -DBUILD_TEST=ON -DCMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH=/path/to/blas/library -DEXT_BLAS_LIBNAME=blas_lib_name
        -DBLAS_HEADER_PATH="<path to BLIS header file blis.h>"
    
    This will enable aocl progress feature tests, main testsuite and legacy libflame test suites. It will generate test_libFLAME_aocl , test_lapack.x ,test_libFLAME executables in the respective directories.
    Note: Building tests require path to an external blas library. Refer to Readme in respective test suite directory for more details
    Recomended to use blis sharedlib with libflame sharedlib

## 3. ENABLE TRACE and LOGS
    User may also enable trace and logs by passing
        -DENABLE_AOCL_DTL=[OPTION]
    along with setting the value of Macros AOCL_DTL_TRACE_ENABLE and AOCL_DTL_LOG_ENABLE to 1 in file libflame/src/aocl_dtl/aocldtlcf.h 
    e.g.
        cmake ../ -DENABLE_ILP64=OFF -DENABLE_AMD_FLAGS=ON -DBUILD_TEST=ON -DENABLE_AOCL_DTL=[DTL_OPTION] -DCMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH=path/to/blas/lib -DEXT_BLAS_LIBNAME=<BLAS_lib_name> -DCMAKE_INSTALL_PREFIX=<path> -DBLAS_HEADER_PATH="<path to BLIS header file blis.h>"
    DTL_OPTION
    1. "ALL" to ENABLE TRACE and LOG
	2. "TRACE" to ENABLE TRACE
	3. "LOG" to ENABLE LOGS
	4. "OFF" to Disable trace and log, if -DENABLE_AOCL_DTL is not passed with the cmake command, DTL is turned off

## 4. using an external Lapack library to run tests
    In order to run tests on an external lapack library an additional option 
    -DEXT_LAPACK_LIBRARY_PATH="path/to/external/lapack/library" and -DEXT_LAPACK_LIBNAME="NAME_OF_THE_LAPACK_LIB" can be passed. 
    if the above options are left blank libflame library will be used

    


