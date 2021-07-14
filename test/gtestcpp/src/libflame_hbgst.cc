/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hbgst.cc
 *  @brief Test application to validate hbgst() using CPP template interface
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hbgst_test is function template for hbgst() functions.
            T can be scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  hbgst_test is function template for hbevx() functions.
	  T can be scomplex, dcomplex
	  
	  hbgst_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output.
	  
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d3/db9/group__complex_o_t_h_e_rcomputational_ga808bf06bc4d353a18ab94f5eaf7c67f0.html#ga808bf06bc4d353a18ab94f5eaf7c67f0
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d0/da6/group__complex16_o_t_h_e_rcomputational_ga4c139408320128b94a42695614ae2646.html#ga4c139408320128b94a42695614ae2646
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return DOUBLE
          Returns Differences value after comparing output of C and CPP based
          library APIs.
 * */
template< typename T, typename Ta >
double hbgst_test(int ip)
{
  typedef integer (*Fptr_NL_LAPACKE_hbgst)(char* vect, char* uplo, integer* n,
                      integer* ka, integer* kb, T* ab, integer* ldab,  T* bb,
                      integer* ldbb, T* x, integer* ldx, T* work, Ta* rwork,
                      integer* info);
  Fptr_NL_LAPACKE_hbgst HBGST;

  // Initialise random number generators with timestamp.
  srand (time(NULL));
  
  /* N is INTEGER
          The order of the matrices A and B.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.\n");
  }
  
  /* VECT is CHARACTER*1
          = 'N':  do not form the transformation matrix X;
          = 'V':  form X.*/
  char vect = eig_paramslist[ip].vect_hbgst;
  if ((vect != 'N') && (vect != 'V')) {
    PRINTF("vect should be N or V. Please correct the input data.\n");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }

  /* KA is INTEGER
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.*/
  integer ka = 0;
  if (uplo == 'U') {
	  ka = eig_paramslist[ip].sda;
  } else if (uplo == 'L') {
	  ka = eig_paramslist[ip].subda;
  }
  if (ka < 0) {
    PRINTF("ka < 0 but it should be: KA >= 0. Please correct the input" \
           " data.\n");
  }
  
  /* KB is INTEGER
          The number of superdiagonals of the matrix B if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KA >= KB >= 0.*/
  integer kb = 0;
  if (uplo == 'U') {
	  kb = eig_paramslist[ip].sdb;
  } else if (uplo == 'L') {
	  kb = eig_paramslist[ip].subdb;
  }
  if (kb < 0) {
    PRINTF("kb < 0 but it should be: KB >= 0. Please correct the input" \
           " data.\n");
  }
  
  /* LDAB is INTEGER
          The leading dimension of the array AB.  LDAB >= KA+1.*/
  integer ldab = eig_paramslist[ip].ldab_hbgv;
  if (ldab < (ka + 1)) {
    PRINTF("ldab < (ka + 1) but should be: LDAB >= (KA + 1). Please correct" \
           " the input data.\n");
  }
  
  // AB is COMPLEX or COMPLEX*16 array, dimension (LDAB,N)
  T *abbuff, *abrefbuff;
  allocate_init_buffer(abbuff, abrefbuff, ldab * n);
  
  /* LDBB is INTEGER
          The leading dimension of the array BB.  LDBB >= KB+1.*/
  integer ldbb = eig_paramslist[ip].ldbb;
  if (ldbb < (kb + 1)) {
    PRINTF("ldbb < (kb + 1) but it should be: LDBB >= (KB + 1). Please " \
           "correct the input data.\n");
  }
  
  /* BB is COMPLEX or COMPLEX*16 array, dimension (LDBB,N) */
  T *bbbuff, *bbrefbuff;
  allocate_init_buffer(bbbuff, bbrefbuff, ldbb * n);
  
  /* LDX is INTEGER
          The leading dimension of the array X.
          LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise.*/
  integer ldx = eig_paramslist[ip].ldx;
  if (ldx < max(1, n)) {
    PRINTF("ldx < max(1, n) but it should be: LDX >= max(1,N). Please " \
           "correct the input data.\n");
  }
  if ((vect == 'V') && (ldx < 1)) {
    PRINTF("When Vect = V, ldx < 1 but it should be: LDX >=1. Please " \
           "correct the input data.\n");
  }
  
  /* X is COMPLEX or COMPLEX*16 array, dimension (LDX,N)
          If VECT = 'V', the n-by-n matrix X.
          If VECT = 'N', the array X is not referenced.*/
  T *xbuff, *xrefbuff; // output buffer
  allocate_init_buffer(xbuff, xrefbuff, ldx * n);
  
  // WORK is COMPLEX array, dimension (N)
  T *workbuff, *workrefbuff;
  allocate_init_buffer(workbuff, workrefbuff, n);
  
  // RWORK is REAL array, dimension (N)
  Ta *rworkbuff, *rworkrefbuff;
  allocate_init_buffer(rworkbuff, rworkrefbuff, n);
  
  // Call CPP function
  integer info_cpp = -1;
  libflame::hbgst<T, Ta>(&vect, &uplo, &n, &ka, &kb, abbuff, &ldab,
              bbbuff, &ldbb, xbuff, &ldx, workbuff, rworkbuff, &info_cpp);

  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    HBGST = (Fptr_NL_LAPACKE_hbgst)dlsym(lapackModule, "chbgst_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HBGST = (Fptr_NL_LAPACKE_hbgst)dlsym(lapackModule, "zhbgst_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HBGST == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }
  integer info_ref = -1;
  HBGST(&vect, &uplo, &n, &ka, &kb, abrefbuff, &ldab, bbrefbuff, &ldbb,
    xrefbuff, &ldx, workrefbuff, rworkrefbuff, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  double diff = 0.0;
  if ((info_cpp == 0) && (info_ref == 0)) {
    diff =  computeError<T>(ldab, n, abrefbuff, abbuff);
    if (vect == 'V') {
      diff +=  computeError<T>(ldx, n, xbuff, xrefbuff);
    }
    diff +=  computeError<T>(1, n, workbuff, workrefbuff);
    diff +=  computeError<Ta>(1, n, rworkbuff, rworkrefbuff);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
  }
  
  // Free up the buffers
  delete[] abbuff; delete[] abrefbuff;
  delete[] bbbuff; delete[] bbrefbuff;
  delete[] xbuff; delete[] xrefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;

  // Return the difference.
  return abs(diff);
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_hbgst, CHBGST) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    diff = hbgst_test<scomplex, float> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
	  PRINTF("diff: %lf\n", diff);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_hbgst, ZHBGST) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    diff = hbgst_test<dcomplex, double> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
	  PRINTF("diff: %lf\n", diff);
  }
}
