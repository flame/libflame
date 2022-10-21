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

Before running the test suite, we must link BLAS and LAPACK library.

By default, the make file is programmed to look for libflame.a in `../../lib/
x86_64-unknown-linux-gnu` directory for LAPACK library. However, if the users
wish to link different LAPACK library, they must set the envrionment variable `LIB_PATH`
to the install path and `LIBFLAME` to the LAPACK library name like the example given
below.


   $ export LIBFLAME=lapack.a LIB_PATH=/usr/local
   $ make


Alternatively, you may set the `make` variable `LIB_PATH` on the command line as you
execute `make`:


   $ make LIBFLAME=lapack.a LIB_PATH=/usr/local


Similarly, user has to provide the path for BLAS library by setting the environment
varaiable 'LIBBLAS'.

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

      The `input.general.operations` file conatins the list of all the APIs
      supported by the test suite. User can enable/disable the testing of a
      particular API by setting/resetting the corresponding API flag. Below is
      a representative example of the default contents of `input.general.operations`.

         1 geqrf QR factorization (0 = disable; 1 = enable, 2 = run APIs with value 2)
         1 gerqf RQ factorisation (0 = disable; 1 = enable, 2 = run APIs with value 2)
         1 getrf LU factorization (0 = disable; 1 = enable, 2 = run APIs with value 2)


   ## Running tests

   Once `input.general.operations` and config files have been tailored to your liking,
   simply run the test suit executable:


      $ ./test_lapack.x

   Input matrix sizes and other parameters can be configured by the user by changing the
   config files.


2. Command line tests

   This method can be used to test a single API with a single set of parameters. To run 
   this mode, name of  the API and corresponding parameters are to be specified as
   command line arguments.

   For example, command-line options for the API GGEVX are:
      ggevx <precisions - sdcz> <balanc> <jobvl> <jobvr> <sense> <N> <LDA> <LDB> <LDVL> <LDVR> <LWORK> <repeats>

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
   by setting the environment variable FLA_TEST_NUM_THREADS to value greter than 1. Both
   the config based and command-line based tests can be used to run this test as given
   below:
      $ FLA_TEST_NUM_THREADS = 4 ./test_lapack.x
      $ FLA_TEST_NUM_THREADS = 4 ./test_lapack.x GGEVX d P N N E 10 10 10 10 10 -1 100
   These runs create 4 threads and call APIs from all the threads.
