###############################################################################
libFLAME test suite usage guidelines
###############################################################################

## Introduction

This wiki explains how to use the test suite included with libFLAME.

The test suite directory has the following contents,
   1. Config - This folder contains config files to control the input combinations to
      test different set of APIs. See ReadMe under this folder for more information.
   2. input.global.operations - This file controls the list of APIs to be tested.
   3. Makefile - Controls how the test suite executable is compiled and linked.
   4. obj - The object files upon being built are placed in this folder.
   5. src - This folder contains the source code.
   6. test_libflame.x - Test suite executable is created upon successfull build
      completion.



## Compiling

Before running the test suite, we must link BLAS and LAPACK library.

By default, the make file is programmed to look for libflame.a in `../../lib/
x86_64-unknown-linux-gnu` directory for LAPACK library. However, If the users
wish to link different LAPACK library, you must set the envrionment variable `LIB_PATH`
to the install path and `LIBFLAME` to the LAPACK library name like the example given
below.

   ```
   $ export LIBFLAME=lapack.a LIB_PATH=/usr/local
   $ make
   ```
   
Alternatively, you may set the `make` variable `LIB_PATH` on the command line as you
execute `make`:

   ```
   $ make LIBFLAME=lapack.a LIB_PATH=/usr/local
   ```

Similarly, user has to provide the path for BLAS library by setting the environment
varaiable 'LIBBLAS'.

When you are ready to compile, simply run `make` from the current directory.

After `make` is complete, an executable named `test_libflame.x` is created.



## Selecting APIs for testing

### `input.general.operations`

The `input.general.operations` file conatins the list of all the APIs supported by the
test suite. User can enable/disable the testing of a particular API by setting/resetting
the corresponding API flag. Below is a representative example of the default contents of
`input.general.operations`.

```
1   geqrf   QR factorization                              (0 = disable; 1 = enable)
1   gerqf   RQ factorisation                              (0 = disable; 1 = enable)
1   gerq2   RQ factorisation with unblocked algorithm     (0 = disable; 1 = enable)
1   potrf   Cholesky factorisation                        (0 = disable; 1 = enable)
1   getrf   LU factorization                              (0 = disable; 1 = enable)
```


## Running tests

Once `input.general.operations` and Config files have been tailored to your liking,
simply run the test suit executable:

```
$ ./test_libflame.x
```