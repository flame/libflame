###############################################################################
Libflame C++ interface usage guidelines
###############################################################################

Follow below steps to build test cases using C++ interface:
1. The header files for C++ interface are present in src/src_cpp/ path.

2. The test cases using C++ interface are present in test/gtestcpp/ path.
Example: libflame_hbev.cc, libflame_hbevx.cc

3. Open test/gtestcpp/Makefile and change the below paths where BLIS & Libflame built files are present:
By default, they are set within project directory:
BLIS path:
LIBBLAS_PATH   := ../../blis/lib/
LIBBLAS        := $(LIBBLAS_PATH)/libblis-mt.a

Note: libblis-mt.a is built file if multi threading is enabled while build.
libblis.a is built file if multi threading is disabled while build.

Libflame path:
LIBFLAMELIB    := ../../lib/x86_64-unknown-linux-gnu/libflame.a

GoogleTest path:
GTESTDIR       := ../../googletest/

4. Also mention which test case file to build in Makefile:
By default, it has hbev and hbevx test cases:
libflame: libflame_hbev.x \
libflame_hbevx.x

Note: For example, in order to build hbev test case then add libflame_hbev.x in above changes.

5. Also change the LAPACK and BLAS built .so file paths (as needed) in test/gtestcpp/src/main.h:
char NETLIB_LAPACK_LIB[60] = "../../lapack/lib/liblapack.so";
char NETLIB_BLAS_LIB[60] = "../../blis/lib/libblis-mt.so";

6. test/gtestcpp/config/EIG_PARAMS.dat has data input for test cases.

7. Test case threshold is defined in test/gtestcpp/src/main.h which can be modified as needed:
Example: SYM_EIGEN_THRESHOLD

8. Cmd to run test case:
Goto project path and /test/gtestcpp/ path where Makefile is present
Run "make" cmd to build the test cases mentioned in Makefile

Run "make clean" cmd to clean built object/executable files of test cases.