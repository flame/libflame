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

INCLUDE_DIR path:
INCLUDE_DIR     := ../../src/aocl_dtl -I../../windows/include -I../../config/x86_64-unknown-linux-gnu

Change the include file path where FLAME.h, FLA_f2c.h, blis1.h files are present.
In Windows build, these files are present in libflame/windows/include/
In Linux build, these files are present in libflame/include/x86_64-unknown-linux-gnu/

4. Also mention which test case file to build in Makefile:
By default, it has hbev and hbevx test cases:
libflame: libflame_hbev.x \
libflame_hbevx.x

Note: For example, in order to build hbev test case then add libflame_hbev.x in above changes.

5. Also change the LAPACK and BLAS built .so file paths (as needed) in test/gtestcpp/src/main.h:
char NETLIB_LAPACK_LIB[60] = "../../lapack/lib/liblapack.so";
char NETLIB_BLAS_LIB[60] = "../../blis/lib/libblis-mt.so";

6. test/gtestcpp/config/EIG_PARAMS.dat etc has data input for test cases.

7. Test case threshold is defined in test/gtestcpp/src/main.h which can be modified as needed:
Example: SYM_EIGEN_THRESHOLD

8. Cmd to run test case:
Goto project path and /test/gtestcpp/ path where Makefile is present
Run "make" cmd to build the test cases mentioned in Makefile

Run "make clean" cmd to clean built object/executable files of test cases.

Downloading and building other dependent projects:
--------------------------------------------------
1. GoogleTest:
--------------
Steps to clone and build GoogleTest code:
a. git clone https://github.com/google/googletest.git -b release-1.10.0
b. cd googletest        # Main directory of the cloned repository.
c. mkdir build          # Create a directory to hold the build output.
d. cd build
e. cmake ..             # Generate native build scripts for GoogleTest.
or cmake .. -DBUILD_GMOCK=OFF
f. make

If facing some problem in git clone:
fatal: unable to access 'https://github.com/google/googletest.git/': SSL certificate problem: certificate is not yet valid
try this cmd: git config --global http.sslVerify false
And follow above steps again.

2. For BLIS, Libflame, LAPACK project build procedure, please refer to ReadMe or notes file of respective projects.

Notes:
------
1. Used GoogleTest project with branch 'release-1.10.0'
2. Used Lapacke v3.9.0 for liblapack.so
