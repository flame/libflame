/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hbevd.cc
 *  @brief Test application to validate hbevd() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hbevd_test is function template for hbevd() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hbevd_test is function template for hbevd() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  hbevd_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output.
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d9/d98/group__complex_o_t_h_e_reigen_ga90a780f962a7affa9b0ee155985e80ab.html#ga90a780f962a7affa9b0ee155985e80ab
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/db/d61/group__complex16_o_t_h_e_reigen_ga3854721b34946230321131ed9e3c08e3.html#ga3854721b34946230321131ed9e3c08e3
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return DOUBLE
          Returns Differences value after comparing output of C and CPP based
          library APIs.
 * */
template<typename T, typename Ta>
double hbevd_test(int ip)
{
  typedef integer (*Fptr_NL_LAPACKE_hbevd)(char* jobz, char* uplo, integer* n,
                      integer* kd, T* ab, integer* ldab, Ta* w, T* z,
                      integer* ldz, T* work, integer* lwork, Ta* rwork,
                      integer* lrwork, integer* iwork, integer* liwork,
                      integer* info);
  Fptr_NL_LAPACKE_hbevd HBEVD = NULL;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));

  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.*/
  char jobz = eig_paramslist[ip].jobz;
  if ((jobz != 'N') && (jobz != 'V')) {
    PRINTF("jobz should be N or V. Please correct the input data.\n");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("jobz should be N or V. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrix A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.\n");
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
    PRINTF("kd is 0 but should be: KD >= 0. Please correct the input data.\n");
  }
  
  /*LDAB is INTEGER
          The leading dimension of the array AB.  LDAB >= KD + 1.*/
  integer ldab = eig_paramslist[ip].ldab;
  if (ldab < (kd+1)) {
    PRINTF("ldab < (kd+1) but it should be: LDAB >= KD + 1. Please correct" \
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
  
  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = eig_paramslist[ip].lwork_hbevd;
  integer lwork_size = 0;
  
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  
  /* LRWORK is INTEGER.
     The length of the array RWORK.*/
  integer lrwork = eig_paramslist[ip].lrwork_hbevd;
  integer lrwork_size = 0;
  if ((lrwork < -1) || (lrwork == 0)) {
    PRINTF("lrwork is 0 or less than -1 and array cannot be allocated with" \
           "this size. Please change the input data.\n");
  }
  
  /* LIWORK is INTEGER.
     The length of the array IWORK.*/
  integer liwork = eig_paramslist[ip].liwork_hbevd;
  integer liwork_size = 0;
  if ((liwork < -1) || (liwork == 0)) {
    PRINTF("liwork is 0 or less than -1 and array cannot be allocated with" \
           "this size. Please change the input data.\n");
  }
  
  integer info_cpp = -1;
  
  /* If lwork/lrwork/liwork is -1, then call this API with default work
     variables. In return, work variables will be updated with array sizes
     needed.*/
  if ((lwork == -1) || (lrwork == -1) || (liwork == -1)) {
    PRINTF("lwork/lrwork/liwork is -1, so call hbevd to get the array" \
           " sizes.\n");
    T worksize = {0};
    Ta rworksize = 0.0;
    integer iworksize = 0;
    
    // Call CPP function
    libflame::hbevd<T, Ta>(&jobz, &uplo, &n, &kd, abbuff, &ldab, wbuff, zbuff,
                  &ldz, &worksize, &lwork, &rworksize, &lrwork, &iworksize,
                  &liwork, &info_cpp);
    PRINTF("info_cpp: %d, worksize: %f, rworksize: %f, iworksize: %d\n", \
           info_cpp, worksize.real, rworksize, iworksize);
    if (info_cpp == 0) {
      if (lwork == -1) {
        lwork_size = worksize.real;
      }
      if (lrwork == -1) {
        lrwork_size = rworksize;
      }
      if (liwork == -1) {
        liwork_size = iworksize;
      }
    }
  }
  
  // Check if lwork = -1 for intializing random buffer size.
  if (lwork != -1) {
    if (n <= 1) {
      if (lwork < 1) {
        PRINTF("lwork must be atleast 1, Please correct the input data.\n");
      }
    } else {
      if (jobz == 'N') {
        if (lwork < n) {
          PRINTF("lwork must be atleast n, Please correct the input data.\n");
        }
      } else if (jobz == 'V') {
        if (lwork < 2*n*2) {
          PRINTF("lwork must be atleast 2*n*2. Please correct the input" \
                 " data.\n");
        }
      }
    }
    lwork_size = lwork;
  }

  // WORK is COMPLEX or COMPLEX*16  array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, max(1, lwork_size));
  
  // Check if lrwork = -1 for intializing random buffer size.
  if (lrwork != -1) {
    if (n <= 1) {
      if (lrwork < 1) {
        PRINTF("lrwork must be atleast 1, Please correct the input data.\n");
      }
    } else {
      if (jobz == 'N') {
        if (lrwork < n) {
          PRINTF("lrwork must be atleast n, Please correct the input data.\n");
        }
      } else if (jobz == 'V') {
        if (lrwork < (1 + 5*n + 2*n*2)) {
          PRINTF("lrwork must be atleast 1 + 5*n + 2*n*2, Please correct the" \
                  " input data.\n");
        }
      }
    }
    lrwork_size = lrwork;
  }
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (LRWORK)
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, lrwork_size);
  
  // Check if liwork = -1 for intializing random buffer size.
  if (liwork != -1) {
    if ((jobz == 'N') && (n <= 1)) {
      if (liwork < 1) {
        PRINTF("liwork must be atleast 1, Please correct the input data.\n");
      }
    } else if ((jobz == 'V') && (n > 1)) {
      if (liwork < (3+5*n)) {
        PRINTF("liwork must be atleast 3 + 5*n, Please correct the" \
                " input data.\n");
      }
    }
    liwork_size = liwork;
  }
  
  // IWORK is INTEGER array, dimension (MAX(1,LIWORK))
  integer *iworkbuff = NULL, *iworkrefbuff = NULL;
  allocate_init_buffer(iworkbuff, iworkrefbuff, max(1, liwork_size));

  // Call CPP function
  integer info_ref = -1;
  libflame::hbevd<T, Ta>(&jobz, &uplo, &n, &kd, abbuff, &ldab, wbuff, zbuff,
                  &ldz, workbuff, &lwork_size, rworkbuff, &lrwork_size,
                  iworkbuff, &liwork_size, &info_cpp);
  
  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    HBEVD = (Fptr_NL_LAPACKE_hbevd)dlsym(lapackModule, "chbevd_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HBEVD = (Fptr_NL_LAPACKE_hbevd)dlsym(lapackModule, "zhbevd_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HBEVD == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  
  HBEVD(&jobz, &uplo, &n, &kd, abrefbuff, &ldab, wrefbuff, zrefbuff, &ldz,
        workrefbuff, &lwork_size, rworkrefbuff, &lrwork_size, iworkrefbuff,
        &liwork_size, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  double diff = 0.0;
  if ((info_cpp == 0) && (info_ref == 0)) {
    diff = computeError<T>(ldab, n, abrefbuff, abbuff);
    if (jobz == 'V') {
      diff += computeError<T>(ldz, n, zrefbuff, zbuff);
    }
    diff += computeError<Ta>(1, n, wbuff, wrefbuff);
    diff += computeError<T>(1, max(1, lwork_size), workbuff, workrefbuff);
    diff += computeError<Ta>(1, lrwork_size, rworkbuff, rworkrefbuff);
    diff += computeError<integer>(1, max(1, liwork_size), iworkbuff,
                                  iworkrefbuff);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
  }
  
  // Free up the buffers
  delete[] abbuff; delete[] abrefbuff;
  delete[] zbuff; delete[] zrefbuff;
  delete[] wbuff; delete[] wrefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
  delete[] iworkbuff; delete[] iworkrefbuff;
  
  // Return the difference.
  return abs(diff);
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_hbevd, CHBEVD) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    diff = hbevd_test<scomplex, float> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
    PRINTF("diff: %lf\n", diff);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_hbevd, ZHBEVD) {
 double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    diff = hbevd_test<dcomplex, double> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
    PRINTF("diff: %lf\n", diff);
  }
}