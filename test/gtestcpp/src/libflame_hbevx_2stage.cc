/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hbevx.cc
 *  @brief Test application to validate hbevx_2stage() using CPP template
           interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hbevx_2stage_test is function template for hbevx_2stage() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hbevx_2stage_test is function template for hbevx_2stage() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  hbevx_2stage_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output.
	  
    Complex reference:
	  http://www.netlib.org/lapack/explore-html/d9/d98/group__complex_o_t_h_e_reigen_gad8820f8a084e192bc6e3713b766a6cdd.html#gad8820f8a084e192bc6e3713b766a6cdd
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/db/d61/group__complex16_o_t_h_e_reigen_ga42094c90148d11dd8eb7627231144c09.html#ga42094c90148d11dd8eb7627231144c09

    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return DOUBLE
          Returns differences value after comparing output of C and CPP based
          library APIs.
 * */
template<typename T, typename Ta>
double hbevx_2stage_test(int ip)
{
  typedef integer (*Fptr_NL_LAPACKE_hbevx_2stage)(char* jobz, char* range,
                  char* uplo, integer* n, integer* kd, T* ab, integer* ldab,
                  T* q, integer* ldq, Ta* vl, Ta* vu, integer* il, integer* iu,
                  Ta* abstol, integer* m, Ta* w, T* z, integer* ldz, T* work,
                  integer* lwork, Ta* rwork, integer* iwork, integer* ifail,
                  integer* info);
  Fptr_NL_LAPACKE_hbevx_2stage HBEVX_2STAGE = NULL;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.
                  Not available in this release.*/
  char jobz = eig_paramslist[ip].jobz_2stage;
  if (jobz != 'N') {
    PRINTF("jobz should be N as V is not supported by hbevx_2stage API." \
           "Please correct the input data.");
  }
  
  /* RANGE is CHARACTER*1
          = 'A': all eigenvalues will be found;
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found;
          = 'I': the IL-th through IU-th eigenvalues will be found.*/
  char range = eig_paramslist[ip].range;
  
  if ((range != 'A') && (range != 'V') && (range != 'I')) {
    PRINTF("range should be A or V or I. Please correct the input data.");
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
    PRINTF("kd < 0 but should be: KD >= 0. Please correct the input data.");
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
  allocate_init_buffer(abbuff, abrefbuff, ldab * n);
  
  /* LDQ is INTEGER
          The leading dimension of the array Q.  If JOBZ = 'V', then
          LDQ >= max(1,N).*/
  integer ldq = eig_paramslist[ip].ldq;
  
  if ((jobz == 'V') && (ldq < max(1,n))) {
    PRINTF("When jobz is V, ldz < max(1,n) but it should be: ldz >= max(1,n)" \
          ". Please correct the input data.\n");
  }
  
  /* Q is COMPLEX or COMPLEX*16 array, dimension (LDQ, N)
        If JOBZ = 'V', the N-by-N unitary matrix used in the
                        reduction to tridiagonal form.
        If JOBZ = 'N', the array Q is not referenced.*/
  T *qbuff = NULL, *qrefbuff = NULL;
  allocate_init_buffer(qbuff, qrefbuff, ldq * n, 0);
  
  /*  VL is REAL or DOUBLE PRECISION
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'.*/
  Ta vl = (Ta)eig_paramslist[ip].vl;
  
  /* VU is REAL or DOUBLE PRECISION
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'.*/
  Ta vu = (Ta)eig_paramslist[ip].vu;
  
  if ((range == 'V' )&& (vl > vu)) {
    PRINTF("vl > vu but it should be: vl < vu. Please correct the input data.");
  }
  
  /* IL is INTEGER
          If RANGE='I', the index of the
          smallest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'.*/
  integer il = eig_paramslist[ip].il;
  
  /* IU is INTEGER
          If RANGE='I', the index of the
          largest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'.*/
  integer iu = eig_paramslist[ip].iu;
  
  if (range == 'I' ) {
    if (n > 0) {
      if ((1 > il) || (il > iu) || (iu > n)) {
        PRINTF("if n > 0, then 1 <= il <= iu <= n is not satisfied." \
                " Please correct the input data.");
      }
    } else if (n = 0) {
      if ((il != 1) && (iu != 0)) {
        PRINTF ("If n = 0, then il = 1 and iu = 0. Please correct the " \
                "input data.");
      }
    }
  }
  
  // ABSTOL is REAL or DOUBLE PRECISION
  Ta abstol = eig_paramslist[ip].abstol;
  
  /* M is INTEGER
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.*/
  integer m = 0, mref = 0;
  
  if (range == 'A') {
    m = mref = n;
  } else if (range == 'I') {
    m = mref = iu-il+1;
  }
  integer mtemp = m; // For comparing Z buffer.
  
  /* W is REAL or DOUBLE PRECISION array, dimension (N)
          The first M elements contain the selected eigenvalues in
          ascending order.
     Allocate & initialize the array.*/
  Ta *wbuff = NULL, *wrefbuff = NULL;
  allocate_init_buffer(wbuff, wrefbuff, n, 0);
  
  /*LDZ is INTEGER
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
  
  // Z is COMPLEX or COMPLEX*16 array, dimension (LDZ, max(1,M))
  T *zbuff = NULL, *zrefbuff = NULL;
  allocate_init_buffer(zbuff, zrefbuff, ldz * (max(1, m)), 0);
  
  // LWORK is INTEGER. The length of the array WORK. LWORK >= 1, when N <= 1;
  integer lwork = eig_paramslist[ip].lwork;
  integer lwork_size = lwork;
  
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           "this size. Please change the input data.\n");
  }
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (7*N)
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, 7 * n, 0);
  
  // IWORK is INTEGER array, dimension (5*N)
  integer *iworkbuff = NULL, *iworkrefbuff = NULL;
  allocate_init_buffer(iworkbuff, iworkrefbuff, 5 * n, 0);
  
  /* IFAIL is INTEGER array, dimension (N). Allocate and initialize
     the output buffer.*/
  integer *ifail = NULL, *ifailref = NULL;
  allocate_init_buffer(ifail, ifailref, n, 0);
  
  integer info_cpp = -1;
  /* If lwork is -1, then call this API with default work variables.
     In return, work variables will be updated with array sizes needed.*/
  if (lwork == -1) {
    PRINTF("lwork is -1, so call hbevx_2stage to get the array sizes.\n");
    T worksize = {0};
    
    libflame::hbevx_2stage<T, Ta>(&jobz, &range, &uplo, &n, &kd, abbuff, &ldab,
                qbuff, &ldq, &vl, &vu, &il, &iu, &abstol, &m, wbuff, zbuff,
                &ldz, &worksize, &lwork, rworkbuff, iworkbuff, ifail,
                &info_cpp);
    PRINTF("info_cpp: %d, worksize: %f\n", info_cpp, worksize.real);
    if (info_cpp == 0) {
      lwork_size = worksize.real;
    }
  }
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (LWORK)
  T *workbuff = NULL, *workrefbuff = NULL;
  if (lwork_size > 0) {
    allocate_init_buffer(workbuff, workrefbuff, lwork_size);
  } else {
    PRINTF("Incorrect LWORK size given/returned to allocate WORK buffer.\n");
  }
  
  // Call CPP function
  info_cpp = -1;
  libflame::hbevx_2stage<T, Ta>(&jobz, &range, &uplo, &n, &kd, abbuff, &ldab,
                qbuff, &ldq, &vl, &vu, &il, &iu, &abstol, &m, wbuff, zbuff,
                &ldz, workbuff, &lwork_size, rworkbuff, iworkbuff, ifail,
                &info_cpp);

  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    HBEVX_2STAGE = (Fptr_NL_LAPACKE_hbevx_2stage)dlsym(lapackModule,
                              "chbevx_2stage_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HBEVX_2STAGE = (Fptr_NL_LAPACKE_hbevx_2stage)dlsym(lapackModule,
                              "zhbevx_2stage_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HBEVX_2STAGE == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  integer info_ref = -1;
  
  HBEVX_2STAGE(&jobz, &range, &uplo, &n, &kd, abrefbuff, &ldab, qrefbuff,
        &ldq, &vl, &vu, &il, &iu, &abstol, &mref, wrefbuff, zrefbuff, &ldz,
        workrefbuff, &lwork_size, rworkrefbuff, iworkrefbuff, ifailref,
        &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  PRINTF ("m: %d, mref: %d\n", m, mref);
  
  // Calculate the differences of buffers.
  double diff = 0.0;
  if ((info_cpp == 0) && (info_ref == 0)) {
    diff = computeError<T>(ldab, n, abrefbuff, abbuff);
    if (jobz == 'V') {
      diff += computeError<T>(ldq, n, qrefbuff, qbuff);
      diff += computeError<integer>(1, n, ifailref, ifail);
      diff +=  computeError<T>(ldz, max(1, mtemp), zrefbuff, zbuff);
                  // Using mtemp, because m will be modified by lapacke func.
    }
    diff += computeError<integer>(1, 1, &mref, &m);
    diff +=  computeError<Ta>(1, n, wrefbuff, wbuff);
    diff += computeError<T>(1, lwork_size, workbuff, workrefbuff);
    diff += computeError<Ta>(7, n, rworkbuff, rworkrefbuff);
    diff += computeError<integer>(5, n, iworkbuff, iworkrefbuff);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
  }
  
  // Free up the buffers
  delete[] abbuff; delete[] abrefbuff;
  delete[] qbuff; delete[] qrefbuff;
  delete[] ifail; delete[] ifailref;
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
TEST(LAPACKCPP_hbevx_2stage, CHBEVX_2STAGE) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
	  diff = hbevx_2stage_test<scomplex, float> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
    PRINTF("diff: %lf\n", diff);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_hbevx_2stage, ZHBEVX_2STAGE) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    diff = hbevx_2stage_test<dcomplex, double> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
    PRINTF("diff: %lf\n", diff);
  }
}