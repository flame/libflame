/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hbgvx.cc
 *  @brief Test application to validate hbgvx() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hbgvx_test is function template for hbgvx() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hbgvx_test is function template for hbgvx() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  hbgvx_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.
	  
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d9/d98/group__complex_o_t_h_e_reigen_gadc31cb38fcdd9a79731c5ca924c1756c.html#gadc31cb38fcdd9a79731c5ca924c1756c
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/db/d61/group__complex16_o_t_h_e_reigen_gae902a3301300bbb1a46459cb990192ad.html#gae902a3301300bbb1a46459cb990192ad
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template<typename T, typename Ta>
void hbgvx_test(int ip)
{
  typedef integer (*fptr_NL_LAPACK_hbgvx)(char* jobz, char* range, char* uplo,
                      integer* n, integer* ka, integer* kb, T* ab,
                      integer* ldab, T* bb, integer* ldbb, T* q, integer* ldq,
                      Ta* vl, Ta* vu, integer* il, integer* iu, Ta* abstol,
                      integer* m, Ta* w, T* z, integer* ldz, T* work,
                      Ta* rwork, integer* iwork, integer* ifail,
                      integer* info);
  fptr_NL_LAPACK_hbgvx hbgvx_ref;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.*/
  char jobz = eig_paramslist[ip].jobz;
  if ((jobz != 'N') && (jobz != 'V')) {
    PRINTF("jobz should be N or V. Please correct the input data.\n");
  }
  
  /* RANGE is CHARACTER*1
          = 'A': all eigenvalues will be found;
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found;
          = 'I': the IL-th through IU-th eigenvalues will be found.*/
  char range = eig_paramslist[ip].range;
  if ((range != 'A') && (range != 'V') && (range != 'I')) {
    PRINTF("range should be A or V or I. Please correct the input data.\n");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A and B are stored;
          = 'L':  Lower triangle of A and B are stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrix A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.");
  }
  
  /* KA is INTEGER
	    The number of superdiagonals of the matrix A if UPLO = 'U',
	    or the number of subdiagonals if UPLO = 'L'.  KD >= 0.*/
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
	    or the number of subdiagonals if UPLO = 'L'.  KD >= 0.*/
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
          The leading dimension of the array AB.  LDAB >= KA + 1.*/
  integer ldab = eig_paramslist[ip].ldab;
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
  
  // BB is COMPLEX or COMPLEX*16 array, dimension (LDBB,N)
  T *bbbuff, *bbrefbuff;
  allocate_init_buffer(bbbuff, bbrefbuff, ldbb * n);
  
  /* LDQ is INTEGER
	  The leading dimension of the array Q.  If JOBZ = 'N',
	  LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).*/
  integer ldq = eig_paramslist[ip].ldq_hbgvx;
  if ((jobz == 'N') && (ldq < 1)) {
    PRINTF("If jobz= N and ldq is < 1, but should be LDQ >= 1. Please " \
           "correct the input data.\n");
  }
  if ((jobz == 'V') && (ldq < max(1, n))) {
    PRINTF("If jobz= N and ldq is < max(1, n), but should be " \
           " LDQ >= max(1, n). Please correct the input data.\n");
  }
  
  T *qbuff = NULL, *qrefbuff = NULL;
  // Q is COMPLEX or COMPLEX*16 array, dimension (LDQ, N)
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
  if ((range == 'V') && (vl >= vu)) {
    PRINTF("If range = V and vl > vu, but should be VL <= VU. Please " \
           "correct the input data.\n");
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
  
  if (range == 'I') {
    if (n > 0) {
      if ((1 > il) || (il > iu) || (iu > n)) {
        PRINTF("If n > 0, then 1 <= il <= iu <= n is not satisfied. Please " \
               "correct the input data.");
      }
    } else if (n = 0) {
      if ((il != 1) && (iu != 0)) {
        PRINTF("If n = 0, then il = 1 and iu = 0. Please correct the " \
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
    m = mref = iu - il + 1;
  }
  
  /* W is REAL or DOUBLE PRECISION array, dimension (N)
	    If INFO = 0, the eigenvalues in ascending order.
     Allocate & initialize the array.*/
  Ta *wbuff, *wrefbuff;
  allocate_init_buffer(wbuff, wrefbuff, n, 0);
  
  /*LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= N.*/
  integer ldz = eig_paramslist[ip].ldz;
  if (ldz < 1) {
    PRINTF("ldz < 1 but it should be: ldz >= 1. Please correct the input" \
          " data.\n");
  }
  if ((jobz == 'V') && (ldz < n)) {
    PRINTF("When jobz is V, ldz < n but it should be: ldz >= n." \
          "Please correct the input data.\n");
  }
  
  // Z is COMPLEX or COMPLEX*16 array, dimension (LDZ, N))
  T *zbuff, *zrefbuff;
  allocate_init_buffer(zbuff, zrefbuff, ldz * n, 0);
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (N)
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, n, 0);
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (7*N)
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, 7 * n, 0);
  
  // IWORK is INTEGER array, dimension (5*N)
  integer *iworkbuff = NULL, *iworkrefbuff = NULL;
  allocate_init_buffer(iworkbuff, iworkrefbuff, 5 * n, 0);
  
  /* IFAIL is INTEGER array, dimension (N). Allocate and initialize the output
     buffer.*/
  integer *ifail = NULL, *ifailref = NULL;
  allocate_init_buffer(ifail, ifailref, n, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("jobz = %c\n", jobz);
    PRINTF("range = %c\n", range);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("ka = %d\n", ka);
    PRINTF("kb = %d\n", kb);
    PRINTF("ldab = %d\n", ldab);
    PRINTF("Size of AB array (ldab*n) = %d\n", ldab * n);
    PRINTF("ldbb = %d\n", ldbb);
    PRINTF("Size of BB array (ldbb*n) = %d\n", ldbb * n);
    PRINTF("ldq = %d\n", ldq);
    PRINTF("Size of Q array (ldq*n) = %d\n", ldq * n);
    PRINTF("vl = %f\n", vl);
    PRINTF("vu = %f\n", vu);
    PRINTF("il = %d\n", il);
    PRINTF("iu = %d\n", iu);
    PRINTF("abstol = %f\n", abstol);
    PRINTF("m = %d\n", m);
    PRINTF("Size of W array (n) = %d\n", n);
    PRINTF("ldz = %d\n", ldz);
    PRINTF("Size of Z array (ldz*n) = %d\n", ldz * n);
    PRINTF("Size of WORK array (n)) = %d\n", n);
    PRINTF("Size of RWORK array (7*n) = %d\n", 7 * n);
    PRINTF("Size of IWORK array (5*n) = %d\n", 5 * n);
    PRINTF("Size of IFAIL array (n) = %d\n", n);
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
    
    // Prints AB array contents
    strncpy(arrayname, "AB input", arraysize);
    print_array<T>(arrayname, abbuff, ldab * n);
    strncpy(arrayname, "AB ref input", arraysize);
    print_array<T>(arrayname, abrefbuff, ldab * n);
    
    // Prints BB array contents
    strncpy(arrayname, "BB input", arraysize);
    print_array<T>(arrayname, bbbuff, ldbb * n);
    strncpy(arrayname, "BB ref input", arraysize);
    print_array<T>(arrayname, bbrefbuff, ldbb * n);
    
    // Prints Q array contents
    strncpy(arrayname, "Q input", arraysize);
    print_array<T>(arrayname, qbuff, ldq * n);
    strncpy(arrayname, "Q ref input", arraysize);
    print_array<T>(arrayname, qrefbuff, ldq * n);
    
    // Prints W array contents
    strncpy(arrayname, "W input", arraysize);
    print_array<Ta>(arrayname, wbuff, n);
    strncpy(arrayname, "W ref input", arraysize);
    print_array<Ta>(arrayname, wrefbuff, n);
    
    // Prints Z array contents
    strncpy(arrayname, "Z input", arraysize);
    print_array<T>(arrayname, zbuff, ldz * n);
    strncpy(arrayname, "Z ref input", arraysize);
    print_array<T>(arrayname, zrefbuff, ldz * n);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, n);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, n);
    
    // Prints RWORK array contents
    strncpy(arrayname, "RWORK input", arraysize);
    print_array<Ta>(arrayname, rworkbuff, 7 * n);
    strncpy(arrayname, "RWORK ref input", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, 7 * n);
    
    // Prints IWORK array contents
    strncpy(arrayname, "IWORK input", arraysize);
    print_array<integer>(arrayname, iworkbuff, 5 * n);
    strncpy(arrayname, "IWORK ref input", arraysize);
    print_array<integer>(arrayname, iworkrefbuff, 5 * n);
    
    // Prints IFAIL array contents
    strncpy(arrayname, "IFAIL input", arraysize);
    print_array<integer>(arrayname, ifail, n);
    strncpy(arrayname, "IFAIL ref input", arraysize);
    print_array<integer>(arrayname, ifailref, n);
  #endif
  
  integer info_cpp = -1;
  // Call CPP function
  libflame::hbgvx<T, Ta>(&jobz, &range, &uplo, &n, &ka, &kb, abbuff, &ldab,
              bbbuff, &ldbb, qbuff, &ldq, &vl, &vu, &il, &iu, &abstol, &m,
              wbuff, zbuff, &ldz, workbuff, rworkbuff, iworkbuff, ifail,
              &info_cpp);

  // Call C function
  
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    hbgvx_ref = (fptr_NL_LAPACK_hbgvx)dlsym(lapackModule, "chbgvx_");
  } else if (typeid(T) == typeid(dcomplex)) {
    hbgvx_ref = (fptr_NL_LAPACK_hbgvx)dlsym(lapackModule, "zhbgvx_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (NULL == hbgvx_ref) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  integer info_ref = -1;
  hbgvx_ref(&jobz, &range, &uplo, &n, &ka, &kb, abrefbuff, &ldab, bbrefbuff, &ldbb,
    qrefbuff, &ldq, &vl, &vu, &il, &iu, &abstol, &mref, wrefbuff, zrefbuff,
    &ldz, workrefbuff, rworkrefbuff, iworkrefbuff, ifailref, &info_ref);
  PRINTF ("info_cpp: %u, info_ref: %u\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  if ((info_cpp >= 0) && (info_ref >= 0)) {
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
      // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
      PRINTF("\nPrinting all Output arrays contents...\n");
      // Prints AB array contents
      strncpy(arrayname, "AB output", arraysize);
      print_array<T>(arrayname, abbuff, ldab * n);
      strncpy(arrayname, "AB ref output", arraysize);
      print_array<T>(arrayname, abrefbuff, ldab * n);
      
      // Prints BB array contents
      strncpy(arrayname, "BB output", arraysize);
      print_array<T>(arrayname, bbbuff, ldbb * n);
      strncpy(arrayname, "BB ref output", arraysize);
      print_array<T>(arrayname, bbrefbuff, ldbb * n);
      
      // Prints Q array contents
      strncpy(arrayname, "Q output", arraysize);
      print_array<T>(arrayname, qbuff, ldq * n);
      strncpy(arrayname, "Q ref output", arraysize);
      print_array<T>(arrayname, qrefbuff, ldq * n);
      
      // Prints M value after API call.
      PRINTF("m = %d, mref = %d\n", m, mref);
      
      // Prints W array contents
      strncpy(arrayname, "W output", arraysize);
      print_array<Ta>(arrayname, wbuff, n);
      strncpy(arrayname, "W ref output", arraysize);
      print_array<Ta>(arrayname, wrefbuff, n);
      
      // Prints Z array contents
      strncpy(arrayname, "Z output", arraysize);
      print_array<T>(arrayname, zbuff, ldz * n);
      strncpy(arrayname, "Z ref output", arraysize);
      print_array<T>(arrayname, zrefbuff, ldz * n);
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, n);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, n);
      
      // Prints RWORK array contents
      strncpy(arrayname, "RWORK output", arraysize);
      print_array<Ta>(arrayname, rworkbuff, 7 * n);
      strncpy(arrayname, "RWORK ref output", arraysize);
      print_array<Ta>(arrayname, rworkrefbuff, 7 * n);
      
      // Prints IWORK array contents
      strncpy(arrayname, "IWORK output", arraysize);
      print_array<integer>(arrayname, iworkbuff, 5 * n);
      strncpy(arrayname, "IWORK ref output", arraysize);
      print_array<integer>(arrayname, iworkrefbuff, 5 * n);
      
      // Prints IFAIL array contents
      strncpy(arrayname, "IFAIL output", arraysize);
      print_array<integer>(arrayname, ifail, n);
      strncpy(arrayname, "IFAIL ref output", arraysize);
      print_array<integer>(arrayname, ifailref, n);
    #endif
    
    double diff = computeError<T>(ldab, n, abrefbuff, abbuff);
    diff += computeError<T>(ldab, n, bbrefbuff, bbbuff);
    diff += computeError<integer>(1, 1, &mref, &m);
    if (jobz == 'V') {
      diff += computeError<T>(ldq, n, qrefbuff, qbuff);
      diff += computeError<integer>(1, n, ifailref, ifail);
      diff += computeError<T>(ldz, n, zrefbuff, zbuff);
    }
    diff += computeError<Ta>(1, n, wrefbuff, wbuff);
    diff += computeError<T>(1, n, workbuff, workrefbuff);
    diff += computeError<Ta>(7, n, rworkbuff, rworkrefbuff);
    diff += computeError<integer>(5, n, iworkbuff, iworkrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, abs(diff), SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  if ((info_cpp > 0) || (info_ref > 0)) {
    if ((info_cpp <= n) || (info_ref <= n)) {
      PRINTF("then i eigenvectors failed to converge. Their" \
             " indices are stored in array IFAIL.\n");
    }
    else if ((info_cpp > n) || (info_ref > n)) {
      PRINTF("if INFO = N + i, for 1 <= i <= N, then CPBSTF" \
             " returned INFO = i: B is not positive definite." \
             " The factorization of B could not be completed and" \
             " no eigenvalues or eigenvectors were computed.\n");
    }
  }
  
  // Free up the buffers
  delete[] abbuff; delete[] abrefbuff;
  delete[] bbbuff; delete[] bbrefbuff;
  delete[] qbuff; delete[] qrefbuff;
  delete[] ifail; delete[] ifailref;
  delete[] zbuff; delete[] zrefbuff;
  delete[] wbuff; delete[] wrefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
  delete[] iworkbuff; delete[] iworkrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_hbgvx, CHBGVX) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hbgvx_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_hbgvx, ZHBGVX) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hbgvx_test<dcomplex, double> (index);
  }
}