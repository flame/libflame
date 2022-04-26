/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_potrf.cc
 *  @brief Test application to validate potrf() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  potrf_test is function template for potrf() functions.
			T can be scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  potrf_test is function template for potrf() functions.
	  T can be float, double scomplex, dcomplex
	  
	  potrf_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

    Real Reference:
    http://www.netlib.org/lapack/explore-html/d8/db2/group__real_p_ocomputational_gaaf31db7ab15b4f4ba527a3d31a15a58e.html#gaaf31db7ab15b4f4ba527a3d31a15a58e
    Double Reference:
    http://www.netlib.org/lapack/explore-html/d1/d7a/group__double_p_ocomputational_ga2f55f604a6003d03b5cd4a0adcfb74d6.html#ga2f55f604a6003d03b5cd4a0adcfb74d6
 	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d8/d6c/group__variants_p_ocomputational_ga4e85f48dbd837ccbbf76aa077f33de19.html#ga4e85f48dbd837ccbbf76aa077f33de19
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d3/d8d/group__complex16_p_ocomputational_ga93e22b682170873efb50df5a79c5e4eb.html#ga93e22b682170873efb50df5a79c5e4eb
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void potrf_test(int ip)
{
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
   /* UPLO is CHARACTER*1
          Specifies whether the upper or lower triangular part of the
          Hermitian matrix A is stored, and how B has been factorized.
          = 'U':  Upper triangular
          = 'L':  Lower triangular*/
  char uplo = lin_solver_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = lin_solver_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
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
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
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
  #endif
  
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
  
  // Calculate the differences of buffers.
  if ((info_cpp >= 0) && (info_ref >= 0)) {
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
      // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
      PRINTF("\nPrinting all Output arrays contents...\n");
      
      // Prints A array contents
      strncpy(arrayname, "A output", arraysize);
      print_array<T>(arrayname, abuff, lda * n);
      strncpy(arrayname, "A ref output", arraysize);
      print_array<T>(arrayname, arefbuff, lda * n);
    #endif
    
    double diff = computeError<T>(lda, n, abuff, arefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
}

/* Use TEST macro and call C++ test function template with
   float as typenames.*/
TEST(LAPACKCPP_potrf, SPOTRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    potrf_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double as typenames.*/
TEST(LAPACKCPP_potrf, DPOTRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    potrf_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_potrf, CPOTRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    potrf_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_potrf, ZPOTRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    potrf_test<dcomplex> (index);
  }
}