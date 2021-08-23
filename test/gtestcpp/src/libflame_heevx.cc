/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_heevx.cc
 *  @brief Test application to validate heevx() using CPP template
           interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  heevx_test is function template for heevx() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  heevx_test is function template for heevx() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  heevx_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d9/de3/group__complex_h_eeigen_ga9f7c713a0119e777afe726e54feb6ef7.html#ga9f7c713a0119e777afe726e54feb6ef7
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaabef68a9c7b10df7aef8f4fec89fddbe.html#gaabef68a9c7b10df7aef8f4fec89fddbe
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template<typename T, typename Ta>
void heevx_test(int ip)
{
  typedef integer (*fptr_NL_LAPACK_heevx)(char* jobz, char* range, char* uplo,
                  integer* n, T* a, integer* lda, Ta* vl, Ta* vu, integer* il,
                  integer* iu, Ta* abstol, integer* m, Ta* w, T* z,
                  integer* ldz, T* work, integer* lwork, Ta* rwork,
                  integer* iwork, integer* ifail, integer* info);
  fptr_NL_LAPACK_heevx heevx_ref;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.*/
  char jobz = eig_paramslist[ip].jobz;
  if ((jobz != 'N') && (jobz != 'V')) {
    PRINTF("jobz should be N or V. Please correct the input data.");
  }
  
  /* RANGE is CHARACTER*1
          = 'A': all eigenvalues will be found.
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found.
          = 'I': the IL-th through IU-th eigenvalues will be found.*/
  char range = eig_paramslist[ip].range;
  if ((range != 'A') && (range != 'V') && (range != 'I')) {
    PRINTF("range should be A or V or I. Please correct the input data.\n");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrix A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.\n");
  }

  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).*/
  integer lda = eig_paramslist[ip].lda;
  if (lda < max(1, n)) {
    PRINTF("lda < max(1, n) but it should be: lda >= max(1, n). Please " \
           "correct the input data.");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA, N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
    
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
  if ((range == 'V') && (vl > vu)) {
    PRINTF("vl > vu but it should be vl < vu when range is V. Please correct the input data.\n");
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
      if (il > iu) {
        PRINTF("il should be <= iu when range is I and n > 0. Please correct the input data.\n");
      }
      if (il < 1) {
        PRINTF("il should be >= 1 when range is I and n > 0. Please correct the input data.\n");
      }
      if (iu > n) {
        PRINTF("iu should be <= n when range is I and n > 0. Please correct the input data.\n");
      }
    } else if (n == 0) {
      il = 1;
      iu = 0;
    }
  }
  
  // ABSTOL is REAL or DOUBLE PRECISION
  Ta abstol = eig_paramslist[ip].abstol;
  
  /* M is INTEGER
		  The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.*/
  integer m = 0;
  if (range == 'A') {
    m = n;
  } else if (range == 'I') {
    m = iu - il + 1;
  }
  
  /*LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= max(1,N).*/
  integer ldz = eig_paramslist[ip].ldz;
  if (ldz < 1) {
    PRINTF("ldz < 1 but it should be: ldz >= 1. Please correct the input" \
          " data.\n");
  }
  if ((jobz == 'V') && (ldz < max(1, n))) {
    PRINTF("When jobz is V, ldz < max(1, n) but it should be: ldz >= max(1, n)." \
          "Please correct the input data.\n");
  }
  
  // Z is COMPLEX array, dimension (LDZ, max(1,M))
  T *zbuff, *zrefbuff;
  allocate_init_buffer(zbuff, zrefbuff, ldz * (max(1, m)), 0);
  
  /* W is REAL or DOUBLE PRECISION array, dimension (N)
          The first M elements contain the selected eigenvalues in
          ascending order.
     Allocate & initialize the array.*/
  Ta *wbuff, *wrefbuff;
  allocate_init_buffer(wbuff, wrefbuff, n, 0);
  
  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = eig_paramslist[ip].lwork_heevr;
  integer lwork_size = lwork;
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  if ((n <= 1) && (lwork < 1)) {
    PRINTF("lwork should be >= 1 when n <= 1.\n");
  }
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (7*N))
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, 7 * n, 0);

  // IWORK is INTEGER array, dimension (5*N)
  integer *iworkbuff = NULL, *iworkrefbuff = NULL;
  allocate_init_buffer(iworkbuff, iworkrefbuff, 5 * n, 0);
  
  integer info_cpp = -1;
  
  // IFAIL is INTEGER array, dimension (N).
  integer *ifail = NULL, *ifailref = NULL;
  allocate_init_buffer(ifail, ifailref, n, 0);
  
  /* If lwork is -1, then call this API with default work variables.
     In return, work variables will be updated with array sizes needed.*/
  if (lwork == -1) {
    PRINTF("lwork/lrwork/liwork is -1, so call heevx() to get the array" \
           " sizes.\n");
    T worksize = {0};
    
    libflame::heevx<T, Ta>(&jobz, &range, &uplo, &n, abuff, &lda, &vl, &vu,
                &il, &iu, &abstol, &m,  wbuff, zbuff, &ldz, &worksize,
                &lwork_size, rworkbuff, iworkbuff, ifail, &info_cpp);
    if (info_cpp == 0) {
      lwork_size = worksize.real;
    }
  }
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, max(1, lwork_size), 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("jobz = %c\n", jobz);
    PRINTF("range = %c\n", range);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("vl = %f\n", vl);
    PRINTF("vu = %f\n", vu);
    PRINTF("il = %d\n", il);
    PRINTF("iu = %d\n", iu);
    PRINTF("abstol = %f\n", abstol);
    PRINTF("m = %d\n", m);
    PRINTF("Size of W array (n) = %d\n", n);
    PRINTF("ldz = %d\n", ldz);
    PRINTF("Size of Z array (ldz*max(1,m)) = %d\n", (ldz * max(1, m)));
    PRINTF("Size of WORK array (MAX(1, LWORK)) = %d\n", max(1, lwork_size));
    PRINTF("LWORK = %d\n", lwork_size);
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
    
    // Prints A array contents
    strncpy(arrayname, "A input", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref input", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
    
    // Prints W array contents
    strncpy(arrayname, "W input", arraysize);
    print_array<Ta>(arrayname, wbuff, n);
    strncpy(arrayname, "W ref input", arraysize);
    print_array<Ta>(arrayname, wrefbuff, n);
    
    // Prints Z array contents
    strncpy(arrayname, "Z input", arraysize);
    print_array<T>(arrayname, zbuff, (ldz * max(1, m)));
    strncpy(arrayname, "Z ref input", arraysize);
    print_array<T>(arrayname, zrefbuff, (ldz * max(1, m)));
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, max(1, lwork_size));
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
    
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
    
    // Prints IFAIl array contents
    strncpy(arrayname, "IFAIL input", arraysize);
    print_array<integer>(arrayname, ifail, n);
    strncpy(arrayname, "IFAIL ref input", arraysize);
    print_array<integer>(arrayname, ifailref, n);
  #endif
  
  info_cpp = -1;
  // Call CPP function
  libflame::heevx<T, Ta>(&jobz, &range, &uplo, &n, abuff, &lda, &vl, &vu,
                &il, &iu, &abstol, &m,  wbuff, zbuff, &ldz, workbuff,
                &lwork_size, rworkbuff, iworkbuff, ifail, &info_cpp);

  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    heevx_ref = (fptr_NL_LAPACK_heevx)dlsym(lapackModule, "cheevx_");
  } else if (typeid(T) == typeid(dcomplex)) {
    heevx_ref = (fptr_NL_LAPACK_heevx)dlsym(lapackModule, "zheevx_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (heevx_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  
  integer info_ref = -1;
  heevx_ref(&jobz, &range, &uplo, &n, arefbuff, &lda, &vl, &vu, &il, &iu, &abstol,
        &m, wrefbuff, zrefbuff, &ldz, workrefbuff, &lwork_size, rworkrefbuff,
        iworkrefbuff, ifailref, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
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
      
      // Prints W array contents
      strncpy(arrayname, "W output", arraysize);
      print_array<Ta>(arrayname, wbuff, n);
      strncpy(arrayname, "W ref output", arraysize);
      print_array<Ta>(arrayname, wrefbuff, n);
      
      // Prints Z array contents
      strncpy(arrayname, "Z output", arraysize);
      print_array<T>(arrayname, zbuff, (ldz * max(1, m)));
      strncpy(arrayname, "Z ref output", arraysize);
      print_array<T>(arrayname, zrefbuff, (ldz * max(1, m)));
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, max(1, lwork_size));
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
      
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
      
      // Prints IFAIl array contents
      strncpy(arrayname, "IFAIL output", arraysize);
      print_array<integer>(arrayname, ifail, n);
      strncpy(arrayname, "IFAIL ref output", arraysize);
      print_array<integer>(arrayname, ifailref, n);
    #endif
    
    double diff = computeError<T>(lda, n, arefbuff, abuff);
    if (jobz == 'V') {
      diff += computeError<T>(ldz, max(1, m), zrefbuff, zbuff);
      diff += computeError<integer>(1, n, ifailref, ifail);
    }
    diff += computeError<Ta>(1, n, wrefbuff, wbuff);
    diff += computeError<T>(1, max(1, lwork_size), workbuff, workrefbuff);
    diff += computeError<Ta>(1, 7*n, rworkbuff, rworkrefbuff);
    diff += computeError<integer>(1, 5*n, iworkbuff, iworkrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, abs(diff), SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
           " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] ifail; delete[] ifailref;
  delete[] zbuff; delete[] zrefbuff;
  delete[] wbuff; delete[] wrefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
  delete[] iworkbuff; delete[] iworkrefbuff;
}

  /* Use TEST macro and call C++ test function template with
     scomplex and float as typenames.*/
TEST(LAPACKCPP_heevx, CHEEVX) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("Index: %d\n", index);
    heevx_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_heevx, ZHEEVX) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("Index: %d\n", index);
    heevx_test<dcomplex, double> (index);
  }
}