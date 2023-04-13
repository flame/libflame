# STEPS TO COMPILE LIBFLAME USING CMAKE
## 1. Generating Libflame library  
    Create a new build directory e.g. build1 
        mkdir build1;
        cd build1;

    Use the following command to configure project
        GCC:
        -----------
        With LP64
            cmake ../ -DENABLE_AMD_FLAGS=ON -DCMAKE_INSTALL_PREFIX=<path> 
        With ILP64
            cmake ../ -DENABLE_ILP64=ON -DENABLE_AMD_FLAGS=ON -DCMAKE_INSTALL_PREFIX=<path>
     
	Note: Use -DCMAKE_C_COMPILER flag to set the compiler
            -DCMAKE_C_COMPILER=gcc OR 
            export CC=gcc

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

    compile library using following command. This will generate libflame.a/libflame.so library in the lib directory
        cmake --build . -j OR make -j

Linking with AOCL Utilities library
------------------------------------
libflame requires AOCL Utilities library "libaoclutils" for certain functions including CPU architecture detection at runtime. The libflame CMake build system, by default, automatically links with libaoclutils library by downloading the source of libaoclutils from AMD GitHub, compiling it and linking/merging with libflame library. However, user can provide an external path for libaoclutils binary and header files via separate flags, 'LIBAOCLUTILS_LIBRARY_PATH' and 'LIBAOCLUTILS_INCLUDE_PATH' respectively. In this scenario, CMake will use the user provided library and does not download libaoclutils source. Following is a sample command for the same
 
$ cmake ../ -DENABLE_AMD_FLAGS=ON -DCMAKE_INSTALL_PREFIX=<path> -DLIBAOCLUTILS_LIBRARY_PATH=<path/to/libaoclutils/library> -DLIBAOCLUTILS_INCLUDE_PATH=<path/to/libaoclutils/header/files>


## 2. Building main Test and AOCL_FLA_PROGRESS Test Suite
    In order to build tests an an additional flag can be set to ON
        -DBUILD_TEST=ON -DCMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH=/path/to/blas/library -DEXT_BLAS_LIBNAME=blas_lib_name
        -DBLAS_HEADER_PATH="<path to BLIS header file blis.h>"
    
    This will enable aocl progress feature tests, main test suite. It will generate test_libFLAME_aocl , test_lapack.x executables in the respective directories.
    Note: Building tests require path to an external blas library. Refer to Readme in respective test suite directory for more details
    Recomended to use blis sharedlib with libflame sharedlib

## 3 Building Legacy test 
    To build Legacy test suite use 
     -DBUILD_LEGACY_TEST=ON -DCMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH=/path/to/blas/library -DEXT_BLAS_LIBNAME=blas_lib_name
    -DBLAS_HEADER_PATH="<path to BLIS header file blis.h>" 
    Note: On Windows, to build and run legacy test suite, a separate macro flag is enabled during libflame library build because of certain constraints in legacy test suite.

## 4. ENABLE TRACE and LOGS
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

## 5. Using an external Lapack library to run tests
    In order to run tests on an external lapack library an additional option 
    -DEXT_LAPACK_LIBRARY_PATH="path/to/external/lapack/library" and -DEXT_LAPACK_LIBNAME="NAME_OF_THE_LAPACK_LIB" can be passed. 
    if the above options are left blank libflame library will be used

## 6. Linking with an external openmp library
    In Order to link with an external openmp library user can pass 
        -DEXT_OPENMP_PATH=<openmp lib path> -DEXT_OPENMP_LIB=<openmp lib name>
    Note: In order to use openmp from the system -DEXT_OPENMP_PATH is to be left blank

    


