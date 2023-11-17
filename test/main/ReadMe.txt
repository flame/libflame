##########################################################################################
libFLAME test suite usage guidelines
##########################################################################################

## Introduction

This wiki explains how to use the test suite included with libFLAME.

The test suite directory has (test/main) the following contents,
   1. config - This folder contains config files to control the input combinations to
      test different set of APIs. See ReadMe under this folder for more information.
   2. input.global.operations - This file controls the list of APIs to be tested.
   3. Makefile - Controls how the test suite executable is compiled and linked.
   4. obj - The object files upon being built are placed in this folder.
   5. src - This folder contains the source code.
   6. test_lapack.x - Test suite executable is created upon successfull build
      completion.


## Compiling

Before running the test suite, we must set BLAS, LAPACK and AOCL-Utils library paths.

BLAS library and header paths has to be set using 'LIBBLAS' and 'BLAS_HEADER_PATH' flags
respectively. AOCL-Utils library path has to be set using LIBAOCLUTILS_LIBRARY_PATH flag.
   $ make BLAS_HEADER_PATH=<path to BLAS API prototypes header file> 
          LIBBLAS=<full path to BLAS library including library file>
          LIBAOCLUTILS_LIBRARY_PATH=<full path to AOCL-Utils library including library file>
 
By default, the make file is programmed to look for libflame.a in `../../lib/
x86_64-unknown-linux-gnu` directory for LAPACK library. However, if the users
wish to link different LAPACK library, they must set the envrionment variable `LIB_PATH`
to the install path and `LIBFLAME` to the LAPACK library name like the example given
below.

   $ export LIBFLAME=lapack.a LIB_PATH=/usr/local

Alternatively, you may set the `make` variable `LIB_PATH` on the command line as you
execute `make`:

   $ make LIBFLAME=lapack.a LIB_PATH=/usr/local 
	  BLAS_HEADER_PATH=<path to BLAS API prototypes header file>
          LIBBLAS=<full path to BLAS library including library file>
          LIBAOCLUTILS_LIBRARY_PATH=<full path to AOCL-Utils library including library file>

When you are ready to compile, simply run `make` from the current directory.
After `make` is complete, an executable named `test_lapack.x` is created.
There are different ways to use the executable to perform different tests as given
below.

1. Config file based tests 

   In this method, input parameters to APIs are taken from config files present in
   'config' folder. The APIs to test are selected from the file 
   'input.general.operations'.

   ## Selecting APIs for testing

      ### `input.general.operations`

      The `input.general.operations` file contains the list of all the APIs
      supported by the test suite. User can enable/disable the testing of a
      particular API by setting/resetting the corresponding API flag. Below is
      a representative example of the default contents of `input.general.operations`.

         1 geqrf QR factorization (0 = disable; 1 = enable, 2 = run APIs with value 2)
         1 gerqf RQ factorisation (0 = disable; 1 = enable, 2 = run APIs with value 2)
         1 getrf LU factorization (0 = disable; 1 = enable, 2 = run APIs with value 2)

   ## Selecting Subgroup of APIs for testing

      ### `input.general.operations`

      The `input.general.operations` file also contains the list of all sub-groups of APIs
      like all LIN,EIG etc. User can enable/disable the testing of a particular sub-group
      of API's by setting/resetting the corresponding group's flag. Below is the content
      in `input.general.operations` corresponding to subgroup testing.

         1   LIN     for testing all LIN API's               (0 = disable; 1 = enable)
         1   EIG     for testing all Eigen API's             (0 = disable; 1 = enable)
         1   SVD     for testing all SVD API's               (0 = disable; 1 = enable)

      Note: If any sub-group is enabled then individual API test will not execute.

   ## Running tests

   Once `input.general.operations` and config files have been tailored to your liking,
   simply run the test suit executable:


      $ ./test_lapack.x

   Input matrix sizes and other parameters can be configured by the user by changing the
   config files. Config files support providing input parameters for four tests. For each
   of the four tests, a range of input dimensions can be specified.

   ## Running test with different config directory.

   This method can be used to test APIs with config files from any directory.
   Name of the directory can be specified through command-line option --config_dir as
     given below.
      $ ./test_lapack.x --config-dir=weekly
   Folder chosen for this option will be 'config/weekly' relative to the test-suite folder.
   The default directory chosen when --config-dir option is not specified is 'config'.

   Note:Directory must be inside 'config' directory.

2. Command line tests

   This method can be used to test a single API with a single set of parameters. To run 
   this mode, name of  the API and corresponding parameters are to be specified as
   command line arguments.

   For example, command-line options for the API GGEVX are:
      ggevx <precisions - sdcz> <balanc> <jobvl> <jobvr> <sense>
            <N> <LDA> <LDB> <LDVL> <LDVR> <LWORK> <repeats>

   Specific instances of calling GGEVX are:
      ./test_lapack.x GGEVX d P N N E 10 10 10 10 10 -1 100
      ./test_lapack.x GGEVX sdcz P N N E 10 10 10 10 10 -1 100
   The first instance tests the double precision GGEVX (DGGEVX) for 10x10 matrices and
   the second runs the same test for all precisions.

   The last parameter 'repeats' is the number of times the API will be called repeatedly.
   This is used to get best performance out of the multiple runs.

   Command-line options for any supported API can be obtained by giving only the API name
   as the only argument.


3. Thread Safety Test

   This test is used to check the thread safety of all supported APIs and can be invoked
   by setting the environment variable FLA_TEST_NUM_THREADS to value greater than 1. Both
   the config based and command-line based tests can be used to run this test as given
   below:
      $ FLA_TEST_NUM_THREADS=4 ./test_lapack.x
   These runs create 4 threads and call APIs from all the threads.


4. Tests from the DTL Logs

   Execute run_tests_from_dtl_logs.py located under test/main/scripts folder with
   following parameters as given below:

      --filename  (required) // The filename or filepath of the DTL log files
      --apiname   (optional) // To execute for Specific API with it's name
      --nrepeats  (optional) // Set number of repeats for the execution of each API

      $ run_tests_from_dtl_logs.py --filename="logs.txt"
      $ run_tests_from_dtl_logs.py --filename="logs.txt" --apiname="ggevx"
      $ run_tests_from_dtl_logs.py --filename="logs.txt" --apiname="ggevx" --nrepeats=4

   To execute the test with thread safety, set the environment variable
   FLA_TEST_NUM_THREADS to a value greater than 1.

   For windows system you need to set the environment variable using the set command and 
   then need to execute the script.

      > set FLA_TEST_NUM_THREADS=4
      > run_tests_from_dtl_logs.py --filename="logs.txt"

   For linux system you can set the environment variable along with the execution of the command

      $ OMP_NUM_THREADS=4 run_tests_from_dtl_logs.py --filename="logs.txt"


5. Unaligned Memory test

   To enable allocate dynamic memory unaligned we need to set below flags while building main testsuite 
       Windows -- FLA_MEM_UNALIGNED is set, unaligned memory is allocated
       Linux   -- MEM_UNALN=1 

## Enabling non-default API naming convention in testsuite:

    LAPACK's default API naming convention is lowercase with underscore. (Ex: getrf_ )
    For enabling UPPERCASE w/, w/o underscore and LOWERCASE w/o underscore API naming
    convention, set API_CALLING_CONVENTION to "upper_","upper","lower" strings respectively
    in test/main/Makefile.
    Testsuite default calling convention is lower_  

NOTE:
   To execute test on windows, its recommended to keep the following in same path/folder:

   1) libflame binary(AOCL-LibFlame-Win-MT-dll or -lib)
   2) test_libFLAME_main.exe
   3) config (default path - test/main/config/)
   4) Specify header file paths of libflame and BLIS to get the LAPACK and BLAS API prototypes respectively
      example: $ cmake -DBLAS_HEADER_PATH="<path to BLIS header file blis.h>"

   To execute the test using libflame shared/dynamic binary, AOCL-LibBlis-Win-MT-dll
   should be in the same path along with the above files.

6. Tests with invalid input parameters using --einfo option:

   Tests to check proper functioning of APIs while sending invalid value for any of the input parameters
   can be done using --einfo option. This option is available only through command-line execution.

   Example:
    ./test_lapack.x GGEVX d P N N E -10 10 10 10 10 -1 100 --einfo=-5

   In the above example, the value of the M has been given -10 which is an invalid value.
   The --einfo parameter states the expected value of 'info' coming out of GGEVX API given the invalid input.
   The test-suite checks the actual value of 'info' against this expected value and reurns PASS if they match
   and FAIL if they don't.

   All parameter related testing commands are compiled in test/main/scripts run_negative_test_cases.py which
   can be used for this purpose.

7. Tests with special inputs using --imatrix option:

   Test the API's by intializing matrix with special input values such as NAN or INFINITY using --imatrix.
   This option is available only through command line execution.

   Example:
    ./test_lapack.x GETRF 10 10 10 1 --imatrix=N
    ./test_lapack.x GETRF 10 10 10 1 --imatrix=I
   
   In the above example passing the value of --imatrix as 'N' will intialize the matrix with NAN values
   and if the value is 'I' then matrix will be intialized with the INFINITY.

8. Tests with -1 for leading dimensions from config files

   When -1 is passed as any of the leading dimensions(lda, ldab, ldu, ldvt, ldz etc) from config files,
   least valid value is assigned to the corresponding leading dimension.

   Example: If lda = -1 passed(from config file) to test_geev API
            then main test-suite sets lda = fla_max(1,n) before calling lapack API geev.

            If lda = -1 is passed through command line, then -1 will be taken as the given lda 
            without any change.
 
9. AOCL_FLA_PROGRESS feature test.

   Enable a macro 'AOCL_FLA_SET_PROGRESS_ENABLE' for aocl progress and build libflame main test suite for sequential/multithread and run the
   executable.

   For testing sequential mode : ./test_lapack.x

   output:
   In AOCL Progress thread  0, at API  DGETRF, progress 8 total threads= 1
   In AOCL Progress thread  0, at API  DGETRF, progress 16 total threads= 1
   In AOCL Progress thread  0, at API  DGETRF, progress 24 total threads= 1
   In AOCL Progress thread  0, at API  DGETRF, progress 32 total threads= 1
   In AOCL Progress thread  0, at API  DGETRF, progress 40 total threads= 1
   In AOCL Progress thread  0, at API  DGETRF, progress 48 total threads= 1
   In AOCL Progress thread  0, at API  DGETRF, progress 56 total threads= 1
   
   For testing multithread mode: FLA_TEST_NUM_THREADS=4 ./test_lapack.x
   
   output:
   In AOCL Progress thread  1, at API  DGETRF, progress 8 total threads= 4
   In AOCL Progress thread  1, at API  DGETRF, progress 16 total threads= 4
   In AOCL Progress thread  2, at API  DGETRF, progress 8 total threads= 4
   In AOCL Progress thread  1, at API  DGETRF, progress 24 total threads= 4
   In AOCL Progress thread  2, at API  DGETRF, progress 16 total threads= 4
   In AOCL Progress thread  1, at API  DGETRF, progress 32 total threads= 4
   In AOCL Progress thread  3, at API  DGETRF, progress 8 total threads= 4
   In AOCL Progress thread  2, at API  DGETRF, progress 24 total threads= 4
   In AOCL Progress thread  0, at API  DGETRF, progress 8 total threads= 4
