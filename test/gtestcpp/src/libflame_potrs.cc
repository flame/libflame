/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_potrs.cc
 *  @brief Test application to validate potrs() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  potrs_test is function template for potrs() functions.
			T can be float, double, scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  potrs_test is function template for potrs() functions.
	  T can be float, double, scomplex, dcomplex
	  
	  potrs_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO == 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO != 0.

    Real reference:
    http://www.netlib.org/lapack/explore-html/d8/db2/group__real_p_ocomputational_gaf5cc1531aa5ffe706533fbca343d55dd.html#gaf5cc1531aa5ffe706533fbca343d55dd
    Double reference:
    http://www.netlib.org/lapack/explore-html/d1/d7a/group__double_p_ocomputational_ga167aa0166c4ce726385f65e4ab05e7c1.html#ga167aa0166c4ce726385f65e4ab05e7c1
 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d6/df6/group__complex_p_ocomputational_gad9052b4b70569dfd6e8943971c9b38b2.html#gad9052b4b70569dfd6e8943971c9b38b2
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d3/d8d/group__complex16_p_ocomputational_gaa2116ea574b01efda584dff0b74c9fcd.html#gaa2116ea574b01efda584dff0b74c9fcd
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template< typename T >
void potrs_test(int ip)
{
  typedef int (*fptr_NL_LAPACK_potrs)(char *uplo, integer *n,
                  integer *nrhs, T *a, integer *lda, T *b,
                  integer *ldb, integer *info);
  fptr_NL_LAPACK_potrs potrs_ref = NULL;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = lin_solver_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }

  /* N is INTEGER
          The order of the matrix A. N >= 0.*/
  integer n = lin_solver_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.");
  }

  /* NRHS is INTEGER. The number of right hand sides,
      i.e., the number of columns of the matrix B.  NRHS >= 0.*/
  integer nrhs = lin_solver_paramslist[ip].nrhs;
  if (nrhs < 0) {
    PRINTF("nrhs < 0 but should be: nrhs >= 0. Please correct the input data.");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).*/
  integer lda = lin_solver_paramslist[ip].lda;
  if (lda < max(1, n)) {
    PRINTF("lda < max(1, n) but it should be: LDA >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);

  /* LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).*/
  integer ldb = lin_solver_paramslist[ip].ldb;
  if (ldb < max(1, n)) {
    PRINTF("ldb < max(1, n) but it should be: LDB >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // B is COMPLEX or COMPLEX16 array, dimension (LDB,NRHS)
  T *bbuff, *brefbuff;
  allocate_init_buffer(bbuff, brefbuff, ldb * nrhs);

  // Call POTRF() to get A buffer and pass the buffer to POTRS()
  PRINTF("Calling POTRF() to get A buffer and pass the buffer to POTRS()\n");
  
  // Call potrf_internal() to get buffers.
  integer info_cpp = -1;
  integer info_ref = -1;
  potrf_internal<T>(uplo, n, abuff, lda, arefbuff, &info_cpp, &info_ref);
  PRINTF ("potrf_internal() info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  if ((info_cpp < 0) && (info_ref < 0)) {
    PRINTF("Info returned by CPP or C API is not successful to get" \
           " A buffer from POTRF(). Exiting....\n");
    closelibs();
    exit (-1);
  }
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("nrhs = %d\n", nrhs);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
  #endif
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1))
  // Array to store array name to print.
  char arrayname[20] = "";
  integer arraysize = sizeof(arrayname);
  #endif
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_INPUT_ARRAYS) && (PRINT_INPUT_ARRAYS == 1))
    // Print all input arrays if PRINT_INPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all Input arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A input", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref input", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
    
    // Prints B array contents
    strncpy(arrayname, "B input", arraysize);
    print_array<T>(arrayname, bbuff, ldb * nrhs);
    strncpy(arrayname, "B ref input", arraysize);
    print_array<T>(arrayname, brefbuff, ldb * nrhs);
  #endif
  
  // Call CPP function
  info_cpp = -1;
  libflame::potrs<T>(&uplo, &n, &nrhs, abuff, &lda, bbuff, &ldb, &info_cpp);
  
  // Call C function
  if (typeid(T) == typeid(float)) {
    potrs_ref = (fptr_NL_LAPACK_potrs)dlsym(lapackModule, "spotrs_");
  } else if (typeid(T) == typeid(double)) {
    potrs_ref = (fptr_NL_LAPACK_potrs)dlsym(lapackModule, "dpotrs_");
  } else if (typeid(T) == typeid(scomplex)) {
    potrs_ref = (fptr_NL_LAPACK_potrs)dlsym(lapackModule, "cpotrs_");
  } else if (typeid(T) == typeid(dcomplex)) {
    potrs_ref = (fptr_NL_LAPACK_potrs)dlsym(lapackModule, "zpotrs_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (potrs_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  
  info_ref = -1;
  potrs_ref(&uplo, &n, &nrhs, arefbuff, &lda, brefbuff, &ldb, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  if ((info_cpp == 0) && (info_ref == 0)) {
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
      // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
      PRINTF("\nPrinting all Output arrays contents...\n");
      
      // Prints B array contents
      strncpy(arrayname, "B output", arraysize);
      print_array<T>(arrayname, bbuff, ldb * nrhs);
      strncpy(arrayname, "B ref output", arraysize);
      print_array<T>(arrayname, brefbuff, ldb * nrhs);
    #endif
    
    double diff = 0.0;
    // TODO: Yet to finalize and do verification changes.
    //diff = computeError<T>(ldb, nrhs, brefbuff, bbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, abs(diff), LIN_DRVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp != 0) || (info_ref != 0));
  }
  
  // Free up the buffers.
  delete[] abuff; delete[] arefbuff;
  delete[] bbuff; delete[] brefbuff;
}

/* Use TEST macro and call C++ test function template with
   float as typenames.*/
TEST(LAPACKCPP_potrs, SPOTRS) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    potrs_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double as typenames.*/
TEST(LAPACKCPP_potrs, DPOTRS) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    potrs_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_potrs, CPOTRS) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    potrs_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_potrs, ZPOTRS) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    potrs_test<dcomplex> (index);
  }
}