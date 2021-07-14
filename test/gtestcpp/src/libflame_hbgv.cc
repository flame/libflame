/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hbgv.cc
 *  @brief Test application to validate hbgv() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hbgv_test is function template for hbgst() functions.
            T can be scomplex, dcomplex
            Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hbgv_test is function template for hbevx() functions.
	  T can be scomplex, dcomplex
    Ta can be float, double.
	  
	  hbgv_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output.
    
    Complex reference:
	  http://www.netlib.org/lapack/explore-html/d9/d98/group__complex_o_t_h_e_reigen_gae30c26efa0a7b94048c00cad17532044.html#gae30c26efa0a7b94048c00cad17532044
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/db/d61/group__complex16_o_t_h_e_reigen_ga76a8cfc758f8dc17ac37f6eed2ef18a4.html#ga76a8cfc758f8dc17ac37f6eed2ef18a4
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return DOUBLE
          Returns Differences value after comparing output of C and CPP based
          library APIs.
 * */
template< typename T, typename Ta >
double hbgv_test(int ip)
{
  typedef integer (*Fptr_NL_LAPACKE_hbgv)(char* jobz, char* uplo, integer* n,
                      integer* ka, integer* kb, T* ab, integer* ldab, T* bb,
                      integer* ldbb, Ta* w, T* z, integer* ldz, T* work,
                      Ta* rwork, integer* info);
  Fptr_NL_LAPACKE_hbgv HBGV = NULL;
  
  // Initialise random number generators with timestamp.
  srand (time(NULL));
  
  /* N is INTEGER
          The order of the matrices A and B.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.\n");
  }
  
  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.*/
  char jobz = eig_paramslist[ip].jobz;
  if ((jobz != 'N') && (jobz != 'V')) {
    PRINTF("jobz should be N or V. Please correct the input data.\n");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangles of A and B are stored;
          = 'L':  Lower triangles of A and B are stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("jobz should be N or V. Please correct the input data.\n");
  }

  /* KA is INTEGER
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'. KA >= 0.*/
  integer ka = 0;
  if (uplo == 'U') {
	  ka = eig_paramslist[ip].sda;
  } else if (uplo == 'L') {
	  ka = eig_paramslist[ip].subda;
  }
  if (ka < 0) {
    PRINTF("ka < 0 but should be: KA >= 0. Please correct the input data.\n");
  }
  
  /* KB is INTEGER
          The number of superdiagonals of the matrix B if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'. KB >= 0.*/
  integer kb = 0;
  if (uplo == 'U') {
	  kb = eig_paramslist[ip].sdb;
  } else if (uplo == 'L') {
	  kb = eig_paramslist[ip].subdb;
  }
  if (kb < 0) {
    PRINTF("kb < 0 but should be: KB >= 0. Please correct the input data.\n");
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
    PRINTF("ldbb < (kb + 1) but should be: LDBB >= (KB + 1). Please correct" \
           " the input data.\n");
  }
  
  /* BB is COMPLEX or COMPLEX*16 array, dimension (LDBB,N) */
  T *bbbuff, *bbrefbuff;
  allocate_init_buffer(bbbuff, bbrefbuff, ldbb * n);
  
  Ta *wbuff, *wrefbuff;
  // W is REAL or DOUBLE PRECISION array, dimension (N)
  allocate_init_buffer(wbuff, wrefbuff, n);
  
  /* LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= N.*/
  integer ldz = eig_paramslist[ip].ldz_hbgv;
  if (ldz < 1) {
    PRINTF("ldz < 1 but should be: LDZ >= 1. Please correct the input data.\n");
  }
  if ((jobz == 'V') && (ldz < n)) {
    PRINTF("If jobz = V and ldz < n but should be: LDZ >= N. Please correct" \
           " the input data.\n");
  }
  
  /* Z is COMPLEX or COMPLEX*16 array, dimension (LDZ, N)
          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
          eigenvectors, with the i-th column of Z holding the
          eigenvector associated with W(i). The eigenvectors are
          normalized so that Z**H*B*Z = I.
          If JOBZ = 'N', then Z is not referenced.*/
  T *zbuff, *zrefbuff; // output buffer
  allocate_init_buffer(zbuff, zrefbuff, ldz * n);
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (N)
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, n);
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (LRWORK)
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, 3 * n);
  
  integer info_cpp = -1;
  // Call CPP function
  libflame::hbgv<T, Ta>(&jobz, &uplo, &n, &ka, &kb, abbuff, &ldab, bbbuff,
              &ldbb, wbuff, zbuff, &ldz, workbuff, rworkbuff, &info_cpp);

  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    HBGV = (Fptr_NL_LAPACKE_hbgv)dlsym(lapackModule, "chbgv_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HBGV = (Fptr_NL_LAPACKE_hbgv)dlsym(lapackModule, "zhbgv_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HBGV == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  integer info_ref = -1;
  
  HBGV(&jobz, &uplo, &n, &ka, &kb, abrefbuff, &ldab, bbrefbuff, &ldbb,
        wrefbuff, zrefbuff, &ldz, workrefbuff, rworkrefbuff, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  double diff = 0.0;
  if ((info_cpp == 0) && (info_ref == 0)) {
    diff =  computeError<T>(ldab, n, abrefbuff, abbuff);
    diff += computeError<T>(ldbb, n, bbrefbuff, bbbuff);
    if (jobz == 'V') {
      diff +=  computeError<T>(ldz, n, zbuff, zrefbuff);
    }
    diff +=  computeError<Ta>(1, n, wbuff, wrefbuff);
    diff +=  computeError<T>(1, n, workbuff, workrefbuff);
    diff +=  computeError<Ta>(3, n, rworkbuff, rworkrefbuff);
  } else {
    PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  // Free up the buffers.
  delete[] abbuff; delete[] abrefbuff;
  delete[] bbbuff; delete[] bbrefbuff;
  delete[] wbuff; delete[] wrefbuff;
  delete[] zbuff; delete[] zrefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
  
  // Return the difference.
  return abs(diff);
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_hbgv, CHBGV) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    diff = hbgv_test<scomplex, float> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
	  PRINTF("diff: %lf\n", diff);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_hbgv, ZHBGV) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    diff = hbgv_test<dcomplex, double> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
	  PRINTF("diff: %lf\n", diff);
  }
}