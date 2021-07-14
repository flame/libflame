/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_heequb.cc
 *  @brief Test application to validate heequb() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  heequb_test is function template for heequb() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  heequb_test is function template for heequb() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  heequb_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output.
	  
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d4/d74/group__complex_h_ecomputational_ga986174490b3d9eb0d10502d96883e153.html#ga986174490b3d9eb0d10502d96883e153
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d3/d80/group__complex16_h_ecomputational_ga0626d54efa3610a40077cf6685df73f1.html#ga0626d54efa3610a40077cf6685df73f1
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
		      Used to pass Index of Eigen Parameters array present in config file.

 * @return DOUBLE
          Returns Differences value after comparing output of C and CPP based
          library APIs.
 * */
template< typename T, typename Ta >
double heequb_test(int ip)
{
  typedef int (*Fptr_NL_LAPACKE_heequb)(char* uplo, integer* n, T* a,
                  integer* lda, Ta* s, Ta* scond, Ta* amax, T* work,
                  integer* info);
  Fptr_NL_LAPACKE_heequb HEEQUB;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
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

  // S is REAL or DOUBLE PRECISION array, dimension (N).
  Ta *sbuff, *srefbuff;
  allocate_init_buffer(sbuff, srefbuff, n);

  // scond and amax are REAL or DOUBLE PRECISION.
  Ta scond = 0.0, scondref = 0.0;
  Ta amax = 0.0, amaxref = 0.0;
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (2*N))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, 2*n);
  
  integer info_cpp = -1;
  
  // Call CPP function
  libflame::heequb<T, Ta>(&uplo, &n, abuff, &lda, sbuff, &scond, &amax,
              workbuff, &info_cpp);

  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    HEEQUB = (Fptr_NL_LAPACKE_heequb)dlsym(lapackModule, "cheequb_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HEEQUB = (Fptr_NL_LAPACKE_heequb)dlsym(lapackModule, "zheequb_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (HEEQUB == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }
  
  integer info_ref = -1;
  // Call C function
  HEEQUB(&uplo, &n, arefbuff, &lda, srefbuff, &scondref, &amaxref,
        workrefbuff, &info_ref);
  PRINTF ("info_cpp: %u, info_ref: %u\n", info_cpp, info_ref);
  PRINTF ("scond: %lf, scondref: %lf\n", scond, scondref);
  PRINTF ("amax: %lf, amaxref: %lf\n", amax, amaxref);
  
  // Calculate the differences of buffers.
  double diff = -1;
  if ((info_cpp == 0) && (info_ref == 0)) {
    diff = computeError<Ta>(1, 1, &amaxref, &amax);
    diff = computeError<Ta>(1, 1, &scondref, &scond);
	  diff += computeError<Ta>(1, n, srefbuff, sbuff);
    diff += computeError<T>(2, n, workbuff, workrefbuff);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] sbuff; delete[] srefbuff;
  delete[] workbuff; delete[] workrefbuff;
  
  // Return the difference.
  return abs(diff);
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_heequb, CHEEQUB) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    diff = heequb_test<scomplex, float> (index);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
    PRINTF("diff: %lf\n", diff);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_heequb, ZHEEQUB) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    diff = heequb_test<dcomplex, double> (index);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
    PRINTF("diff: %lf\n", diff);
  }
}