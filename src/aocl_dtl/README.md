###############################################################################
Guidelines to enable logging and tracing in libFLAME library
###############################################################################

Following are the steps to enable Trace and Log. 

1.  Open header file "libflame/src/aocl_dtl/aocldtlcf.h"
      i.  Enable Trace by making the following macro as 1 :
            #define AOCL_DTL_TRACE_ENABLE       1
      ii.  Enable Log by making the following macro as 1 :
            #define AOCL_DTL_LOG_ENABLE         1
    
2.  After Step 1, build and install libFLAME. AOCL DTL library, "libaocldtl.a", "libaocldtl.so" along with libflame gets installed in the specified            directory. 

3.  To know the libFLAME APIs used by an application and/or inputs for those APIs, please link DTL library along with libFLAME library in your application.    Example: include '-lflame -laocldtl ' in the application build command

4.  Here are are some examples of usage of DTL in libFLAME
    libFLAME Test Suite : For enabling DTL in libFLAME test suite, provide aocldtl library along with blis library in LIBBLAS path as follows
        From "libflame" directory:
            make check LIBBLAS="-fopenmp <AOCLDTL LIB PATH/libaocldtl.a> <BLIS LIB PATH/libblis-mt.a>"
        From "libflame/test" directory: 
            cd test
            make LIBBLAS="-fopenmp <AOCLDTL LIB PATH/libaocldtl.a> <BLIS LIB PATH/libblis-mt.a>"
            ./test_libflame.x
    
5.  After executing the application, linked to DTL library, log and trace files are generated in libflame/test directory. 
    For Example: "P12466_T12466_aocldtl_trace.txt" and "P12466_T12466_aocldtl_log.txt". 
    
6.  Steps to run netlib tests and link the DTL library along with libflame library:
      i.  bash run-netlib-test.sh BLAS_LIB_PATH="<blas library path>" LAPACK_LIB_PATH="<lapack library path>" DTL_LIB_PATH="<dtl library path>" DTL=1
    
      ii. Generates log and trace files in "libflame/netlib-test/libflame_netlib/TESTING" folder as "P5653_T5653_aocldtl_log.txt",                                   "P5653_T5653_aocldtl_trace.txt".  
             
          
