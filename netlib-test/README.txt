
===============================================================================
                    Run netlib lapack tests
===============================================================================

1. Use the script "run-netlib-test.sh" to run netlib lapack tests

2. Usage :
   
   $ sh run-netlib-test.sh BLAS_LIB_PATH=<blas library path> LAPACK_TEST_DIR=<lapack library path> 
                           [BLAS_LIB=<blas library] [LAPACK_LIB=<lapack library>] [ILP64=<0/1>] 
                           [LAPACK_TEST_DIR=<path to netlib lapack test directory>]
   
   [] indicates optional argument

   Example: $ sh run-netlib-test.sh BLAS_LIB_PATH="/home/user/blis/lib" LAPACK_LIB_PATH="/home/user/libflame/lib" BLAS_LIB="libblis-mt.a"                                    LAPACK_LIB="libflame.a"   
     
        BLAS_LIB : blas library to use. Default=libblist-mt.a
        LAPACK_LIB : lapac library to use. Default=libflame.a
        ILP64 : LP64 or ILP64 mode. Default=0(Use LP64)
        BLAS_LIB_PATH : path of blas library chosen in BLAS_LIB
        LAPACK_LIB_PATH : path to lapack library chosen in LAPACK_LIB
        LAPACK_TEST_DIR : netlib lapack test directory name. Default=lapack-3.10.0
