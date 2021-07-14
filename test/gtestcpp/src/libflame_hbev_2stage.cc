/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hbev_2stage.cc
 *  @brief Test application to validate hbev_2stage() using CPP template
           interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief hbev_2stage_test is function template for hbev_2stage() functions.
			T can be scomplex, dcomplex.
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hbev_2stage_test is function template for hbev_2stage() functions.
	  T can be scomplex, dcomplex.
	  Ta can be float, double.
	  
	  hbev_2stage_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output.
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d9/d98/group__complex_o_t_h_e_reigen_ga4ef30f4426bc3e5e88d1c833b53aeadc.html#ga4ef30f4426bc3e5e88d1c833b53aeadc
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/db/d61/group__complex16_o_t_h_e_reigen_gaf637994a7cb287906efc0254d7d58f69.html#gaf637994a7cb287906efc0254d7d58f69
    \endverbatim

 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return DOUBLE
          Returns differences value after comparing output of C and CPP based
          library APIs.
 * */
template<typename T, typename Ta>
double hbev_2stage_test(int ip)
{
  typedef integer (*Fptr_NL_LAPACK_hbev_2stage)(char *jobz, char *uplo,
                      integer *n, integer *kd, T *ab, integer *ldab, Ta *w,
                      T *z, integer *ldz, T* work, integer *lwork, Ta* rwork,
                      integer *info);
  Fptr_NL_LAPACK_hbev_2stage HBEV_2STAGE = NULL;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.
                  Not available in this release.*/
  char jobz = eig_paramslist[ip].jobz_2stage;
  if (jobz != 'N') {
    PRINTF("jobz should be N as V is not supported by hbev_2stage API." \
           "Please correct the input data.\n");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("jobz should be N or V. Please correct the input data.");
  }
  
  /* N is INTEGER
          The order of the matrix A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.");
  }
  
  /* KD is INTEGER
	    The number of superdiagonals of the matrix A if UPLO = 'U',
	    or the number of subdiagonals if UPLO = 'L'.  KD >= 0.*/
  integer kd = 0;
  if (uplo == 'U') {
    kd = eig_paramslist[ip].sda;
  } else if (uplo == 'L') {
    kd = eig_paramslist[ip].subda;
  }
  if (kd < 0) {
    PRINTF("kd is 0 but should be: KD >= 0. Please correct the input data.");
  }
  
  /* LDAB is INTEGER
          The leading dimension of the array AB.  LDAB >= KD + 1.*/
  integer ldab = eig_paramslist[ip].ldab;
  if (ldab < (kd+1)) {
    PRINTF("ldab < (kd+1) but it should be: LDAB >= KD + 1. Please correct" \
          " the input data.\n");
  }
  
  // AB is COMPLEX or COMPLEX*16 array, dimension (LDAB, N)
  T *abbuff = NULL, *abrefbuff = NULL;
  allocate_init_buffer(abbuff, abrefbuff, ldab*n);
  
  /* W is REAL or DOUBLE PRECISION array, dimension (N)
          If INFO = 0, the eigenvalues in ascending order.*/
  Ta *wbuff = NULL, *wrefbuff = NULL;
  allocate_init_buffer(wbuff, wrefbuff, n);
  
  /* LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= max(1,N).*/
  integer ldz = eig_paramslist[ip].ldz;
  if (ldz < 1) {
    PRINTF("ldz < 1 but it should be: ldz >= 1. Please correct the input" \
          " data.\n");
  }
  if ((jobz == 'V') && (ldz < max(1,n))) {
    PRINTF("When jobz is V, ldz < max(1,n) but it should be: ldz >= max(1,n)" \
          ". Please correct the input data.\n");
  }
  
  // Z is COMPLEX or COMPLEX*16 array, dimension (LDZ, N)
  T *zbuff = NULL, *zrefbuff = NULL;
  allocate_init_buffer(zbuff, zrefbuff, ldz*n);
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension max(1,3*N-2)
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, max(1, 3*n-2));
  
  //  LWORK is INTEGER. The length of the array WORK.
  integer lwork = eig_paramslist[ip].lwork;
  integer lwork_size = lwork;
  integer info_cpp = -1;
  
  /* If lwork is -1, then call this API with default work variable.
     In return, work variable will be updated with lwork needed.*/
  if (lwork == -1) {
    PRINTF("lwork is -1, so call hbev_2stage() to get the WORK array size.\n");
    T worksize = {0};
    
    // Call CPP function
    libflame::hbev_2stage<T, Ta>(&jobz, &uplo, &n, &kd, abbuff, &ldab, wbuff,
                  zbuff, &ldz, (T *)&worksize, &lwork, rworkbuff, &info_cpp);
    PRINTF("info_cpp: %d, worksize: %f\n", info_cpp, worksize.real);
    if (info_cpp == 0) {
      lwork_size = worksize.real;
    }
  }
  
  // WORK is COMPLEX or COMPLEX*16 array, dimension (LWORK)
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, lwork_size);
  
  info_cpp = -1;
  // Call CPP function
  libflame::hbev_2stage<T, Ta>(&jobz, &uplo, &n, &kd, abbuff, &ldab, wbuff,
                  zbuff, &ldz, workbuff, &lwork, rworkbuff, &info_cpp);

  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    HBEV_2STAGE = (Fptr_NL_LAPACK_hbev_2stage)dlsym(lapackModule, 
                      "chbev_2stage_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HBEV_2STAGE = (Fptr_NL_LAPACK_hbev_2stage)dlsym(lapackModule, 
                      "zhbev_2stage_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HBEV_2STAGE == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  
  integer info_ref = -1;
  
  // Call C function - NetLib Lapack API
  HBEV_2STAGE(&jobz, &uplo, &n, &kd, abrefbuff, &ldab, wrefbuff, 
      zrefbuff, &ldz, workrefbuff, &lwork, rworkrefbuff, &info_ref);

  PRINTF("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  double diff = 0.0;
  if ((info_cpp == 0) && (info_ref == 0)) {
    diff = computeError<T>(ldab, n, abrefbuff, abbuff);
    if (jobz == 'V') {
      diff += computeError<T>(ldz, n, zrefbuff, zbuff);
    }
    diff += computeError<Ta>(1, n, wbuff, wrefbuff);
    if (lwork_size != 0) {
      diff += computeError<T>(1, lwork_size, workbuff, workrefbuff);
    }
    diff += computeError<Ta>(1, max(1, 3*n-2), rworkbuff, rworkrefbuff);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
           " differences.\n");
  }
  
  // Free up the buffers
  delete[] abbuff; delete[] abrefbuff;
  delete[] zbuff; delete[] zrefbuff;
  delete[] wbuff; delete[] wrefbuff;
  if (lwork_size != 0) {
    delete[] workbuff; delete[] workrefbuff; 
  }
  delete[] rworkbuff; delete[] rworkrefbuff;
  
  // Return the difference.
  return abs(diff);
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_hbev_2stage, CHBEV_2STAGE) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    diff = hbev_2stage_test<scomplex, float> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
    PRINTF("diff: %lf\n", diff);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_hbev_2stage, ZHBEV_2STAGE) {
 double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    diff = hbev_2stage_test<dcomplex, double> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
    PRINTF("diff: %lf\n", diff);
  }
}