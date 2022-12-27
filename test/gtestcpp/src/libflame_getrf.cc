/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_getrf.cc
 *  @brief Test application to validate getrf() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  getrf_test is function template for getrf() functions.
			T can be scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  getrf_test is function template for getrf() functions.
	  T can be scomplex, dcomplex
	  
	  getrf_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

    Real reference:
    http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga8d99c11b94db3d5eac75cac46a0f2e17.html#ga8d99c11b94db3d5eac75cac46a0f2e17
    Double reference:
    https://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html#ga0019443faea08275ca60a734d0593e60
 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d4/d7e/group__complex_g_ecomputational_gaed8e85049ecfb314d259bfdb3908a60d.html#gaed8e85049ecfb314d259bfdb3908a60d
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/df/dc5/group__variants_g_ecomputational_ga5b625680e6251feb29e386193914981c.html#ga5b625680e6251feb29e386193914981c
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void getrf_test(int ip)
{
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* M is INTEGER
          The number of rows of the matrix A.  M >= 0.*/
  integer m = lin_solver_paramslist[ip].m;
  if (m < 0) {
    PRINTF("m should be >=0. Please correct the input data.\n");
  }

  /* N is INTEGER
          The number of columns of the matrix A.  N >= 0.*/
  integer n = lin_solver_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= fla_max(1,M).*/
  integer lda = lin_solver_paramslist[ip].lda_getrf;
  if (lda < fla_max(1, m)) {
    PRINTF("lda < fla_max(1, m) but it should be: LDA >= fla_max(1,M). Please " \
           "correct the input data.\n");
  }
  
  // A is FLOAT/DOUBLE/COMPLEX/COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // IPIV is INTEGER array, dimension (fla_min(M,N))
  integer *ipivbuff, *ipivrefbuff;
  allocate_init_buffer(ipivbuff, ipivrefbuff, fla_min(m, n));
  
  // Call getrf_internal() to get buffers.
  integer info_cpp = -1;
  integer info_ref = -1;
  getrf_internal<T>(m, n, abuff, lda, ipivbuff,
              arefbuff, ipivrefbuff, &info_cpp, &info_ref);
  PRINTF ("getrf_internal info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  if ((info_cpp >= 0) && (info_ref >= 0)) {
    double diff = computeError<T>(lda, n, abuff, arefbuff);
    diff += computeError<integer>(1, fla_min(m, n), ipivbuff, ipivrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] ipivbuff; delete[] ipivrefbuff;
}

/* Use TEST macro and call C++ test function template with
   float as typenames.*/
TEST(LAPACKCPP_getrf, SGETRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getrf_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double as typenames.*/
TEST(LAPACKCPP_getrf, DGETRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getrf_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_getrf, CGETRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getrf_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_getrf, ZGETRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getrf_test<dcomplex> (index);
  }
}