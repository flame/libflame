# AOCL-LAPACK Samples
This directory contains sample source files showing usage of AOCL-LAPACK library functions. 
Use the provided cmake script file "CMakeLists.txt" to compile and run the programs. Same cmake
script can be used for both Linux and Windows platforms.

## Build steps on Linux and windows

1. Move to installed examples directory
        ```
        cd test/examples
        ```
2. Configure the build system to compile sample applications
        ```
        cmake . -DEXT_BLAS_LIBRARY_DEPENDENCY_PATH=/path/to/blas/library
        -DEXT_LAPACK_LIBRARY_PATH=/path/to/libflame/library
	-DAOCLUTILS_LIBRARY_PATH=/path/to/aoclutils/library
        -DEXT_BLAS_LIBNAME=blas_lib_name -DEXT_LAPACK_LIBNAME=libflame_lib_name
        -DEXT_FLAME_HEADER_PATH=/path/to/flame/header/file
        ```
        eg:
        ```
        cmake . -DEXT_BLAS_LIBRARY_DEPENDENCY_PATH=../../lib/ -DEXT_LAPACK_LIBRARY_PATH=../../lib 
        -DAOCLUTILS_LIBRARY_PATH=/home/usr/aoclutils-install/lib 
	-DEXT_BLAS_LIBNAME=libblis-mt.a -DEXT_LAPACK_LIBNAME=libflame.a
        -DEXT_FLAME_HEADER_PATH=../../include/
        ```

3. Compile the sample application
        ```
        For Linux
                cmake --build . or make
        For Windows
                cmake --build .
        ```

4. Run the application
        ```
        For Linux
                ./test_dgetrf.x
        For Windows
                cd Debug
                test_dgetrf.exe
        ```
