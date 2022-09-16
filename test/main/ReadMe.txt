###############################################################################
libFLAME test suite usage guidelines
###############################################################################

## Introduction

This wiki explains how to use the test suite included with libFLAME.

The test suite directory has the following contents,
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
There are two ways to use the executable to perform tests and both are explained below.

1. The first method is to test all the supported APIs by running the executable without
any command-line arguments. The input parameters to the APIs are taken from config
files present under 'config' folder. The APIs to test are selected from the file
'input.general.operations'.


## Selecting APIs for testing

### `input.general.operations`

The `input.general.operations` file conatins the list of all the APIs supported by the
test suite. User can enable/disable the testing of a particular API by setting/resetting
the corresponding API flag. Below is a representative example of the default contents of
`input.general.operations`.


1   geqrf   QR factorization                              (0 = disable; 1 = enable)
1   gerqf   RQ factorisation                              (0 = disable; 1 = enable)
1   gerq2   RQ factorisation with unblocked algorithm     (0 = disable; 1 = enable)
1   potrf   Cholesky factorisation                        (0 = disable; 1 = enable)
1   getrf   LU factorization                              (0 = disable; 1 = enable)



## Running tests

Once `input.general.operations` and config files have been tailored to your liking,
simply run the test suit executable:


$ ./test_lapack.x

Input matrix sizes and other parameters can be configured by the user by changing the
config files.


2. The second method is to run one particular API by giving the name and parameters of
the API as command line arguments.

For example, command-line options for the API GGEVX are:
   ggevx <precisions - sdcz> <balanc> <jobvl> <jobvr> <sense> <N> <LDA> <LDB> <LDVL> <LDVR> <LWORK> <repeats>

Specific instances of calling GGEVX are:
   ./test_lapack.x GGEVX d P N N E 10 10 10 10 10 -1 100
   ./test_lapack.x GGEVX sdcz P N N E 10 10 10 10 10 -1 100
The first instance will test the double precision GGEVX (DGGEVX) for 10x10 matrices and
the second will run the same test for all precisions.

The last parameter 'repeats' is the number of times the API will be called repeatedly.
This is used to get best performance out of the multiple runs.

Command-line options for any supported API can be obtained by giving only the API name
as the only argument.
