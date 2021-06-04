/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hbev.cc
 *  libflame_hbev.cc Test application to validate CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hbev_test is function template for hbev() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hbev_test is function template for hbev() functions.
	  T can be lapack_complex_float, lapack_complex_double
	  Ta can be float, double.
	  
	  hbev_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output.
    
    Complex reference:
    http://www.netlib.org/lapack/explore-html/d9/d98/group__complex_o_t_h_e_reigen_ga0f6d76a9363386f7fe3d13b8e6a19229.html#ga0f6d76a9363386f7fe3d13b8e6a19229
	  Double Complex reference:
    http://www.netlib.org/lapack/explore-html/db/d61/group__complex16_o_t_h_e_reigen_ga72184c03c8976891c11e42f3463c2d38.html#ga72184c03c8976891c11e42f3463c2d38
    \endverbatim
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return DOUBLE
          Returns Differences value after comparing output of C and CPP based
          library APIs.
 * */
template<typename T, typename Ta>
double hbev_test(int ip)
{
  typedef void (*Fptr_NL_LAPACK_hbev)(char *jobz, char *uplo, int *n,
                    int *kd, T *ab, int *ldab, Ta *w, T *z,
                    int *ldz, T* work, Ta* rwork, int *info);
  Fptr_NL_LAPACK_hbev HBEV = NULL;

  // Initialise random number generators with timestamp
  srand (time(NULL));
  
  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.*/
  char jobz = eig_paramslist[ip].jobz;
  
  if ((jobz != 'N') && (jobz != 'V')) {
    printf("jobz should be N or V. Please correct the input data.");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = eig_paramslist[ip].uplo;
  
  if ((uplo != 'U') && (uplo != 'L')) {
    printf("jobz should be N or V. Please correct the input data.");
  }
  
  /* N is INTEGER
          The order of the matrix A.  N >= 0.*/
  int n = eig_paramslist[ip].n;
  
  if (n < 0) {
    printf("n < 0 but should be: n >= 0. Please correct the input data.");
  }
  
  /* KD is INTEGER
	    The number of superdiagonals of the matrix A if UPLO = 'U',
	    or the number of subdiagonals if UPLO = 'L'.  KD >= 0.*/
  int kd = 0;
  if (uplo == 'U') {
    kd = eig_paramslist[ip].sda;
  } else if (uplo == 'L') {
    kd = eig_paramslist[ip].subda;
  }
  
  if (kd < 0) {
    printf("kd is 0 but should be: KD >= 0. Please correct the input data.");
  }
  
  /* LDAB is INTEGER
          The leading dimension of the array AB.  LDAB >= KD + 1.*/
  int ldab = eig_paramslist[ip].ldab;
  
  if (ldab < (kd+1)) {
    printf("ldab < (kd+1) but it should be: LDAB >= KD + 1. Please correct" \
          " the input data.\n");
  }
  
  // AB is COMPLEX or COMPLEX*16 array, dimension (LDAB, N)
  T *abbuff = NULL, *abrefbuff = NULL;
  
  allocate_init_buffer(abbuff, abrefbuff, ldab * n);
  
  /* W is REAL or DOUBLE PRECISION array, dimension (N)
          If INFO = 0, the eigenvalues in ascending order.*/
  Ta *wbuff = NULL, *wrefbuff = NULL;
  allocate_init_buffer(wbuff, wrefbuff, n);
  
  /* LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= max(1,N).*/
  int ldz = eig_paramslist[ip].ldz;
  
  if (ldz < 1) {
    printf("ldz < 1 but it should be: ldz >= 1. Please correct the input" \
          " data.\n");
  }
  
  if ((jobz == 'V') && (ldz < max(1,n))) {
    printf("When jobz is V, ldz < max(1,n) but it should be: ldz >= max(1,n)." \
          "Please correct the input data.\n");
  }
  
  // Z is COMPLEX or COMPLEX*16 array, dimension (LDZ, N)
  T *zbuff = NULL, *zrefbuff = NULL;
  allocate_init_buffer(zbuff, zrefbuff, ldz * n);
  
  // WORK is COMPLEX or COMPLEX*16 array, dimension (N)
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, n, 0);
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (max(1,3*N-2))
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, max(1, 3*n-2), 0);

  int info_cpp = -1, info_ref = -1;
  
  // Call CPP function
  libflame::hbev<T, Ta>(&jobz, &uplo, &n, &kd, abbuff, &ldab, wbuff,
							zbuff, &ldz, workbuff, rworkbuff, &info_cpp);

  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    HBEV = (Fptr_NL_LAPACK_hbev)dlsym(lapackModule, "chbev_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HBEV = (Fptr_NL_LAPACK_hbev)dlsym(lapackModule, "zhbev_");
  } else {
	  printf("Invalid typename is passed to hbev_test function template.\n");
  }
  
  if (HBEV == NULL) {
    printf("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  
  // Call C function - NetLib Lapack API
  HBEV(&jobz, &uplo, &n, &kd, abrefbuff, &ldab, wrefbuff, 
      zrefbuff, &ldz, workrefbuff, rworkrefbuff, &info_ref);

  // Calculate the differences of buffers.
  double diff = -1;
  if ((info_cpp == 0) && (info_ref == 0)) {
    diff = computeError<T>(ldab, n, abrefbuff, abbuff);
    diff += computeError<T>(ldz, n, zrefbuff, zbuff);
    diff += computeError<Ta>(1, n, wbuff, wrefbuff);
    diff += computeError<T>(1, n, workbuff, workrefbuff);
    diff += computeError<Ta>(1, max(1, 3*n-2), rworkbuff, rworkrefbuff);
  } else {
    printf("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
  }
  
  // Free up the buffers
  delete[] abbuff; delete[] abrefbuff;
  delete[] zbuff; delete[] zrefbuff;
  delete[] wbuff; delete[] wrefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
  
  // Return the difference.
  return diff;
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_hbev, CHBEV) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  diff = hbev_test<scomplex, float> (0);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_hbev, ZHBEV) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    diff = hbev_test<dcomplex, double> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  }
}