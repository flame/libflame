/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hegv_2stage.cc
 *  @brief Test application to validate hegv_2stage() using CPP template
           interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hegv_2stage_test is function template for hegv_2stage() functions.
			T can be scomplex, dcomplex
      Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hegv_2stage_test is function template for hegv_2stage() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  hegv_2stage_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d9/de3/group__complex_h_eeigen_ga7027121a5599a07ebab2a56100f05f32.html#ga7027121a5599a07ebab2a56100f05f32
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_ga415db46f3fbb92fd4c7af056e6a040d2.html#ga415db46f3fbb92fd4c7af056e6a040d2
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template< typename T, typename Ta >
void hegv_2stage_test(int ip)
{
  typedef integer (*fptr_NL_LAPACK_hegv_2stage)(integer* itype, char* jobz,
                      char* uplo, integer* n, T* a, integer* lda, T* b,
                      integer* ldb, Ta* w, T* work, integer* lwork,
                      Ta* rwork, integer* info);
  fptr_NL_LAPACK_hegv_2stage hegv_2stage_ref = NULL;
  
  srand (SRAND_SEED_VALUE);
  
  /* ITYPE is INTEGER
          Specifies the problem type to be solved:
          = 1:  A*x = (lambda)*B*x
          = 2:  A*B*x = (lambda)*x
          = 3:  B*A*x = (lambda)*x*/
  int itype = eig_paramslist[ip].itype;
  if ((itype < 1) || (itype > 3)) {
    PRINTF("itype should be from 1 to 3. Please correct the input data.\n");
  }
  
  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.
          Not available in this release.*/
  char jobz = eig_paramslist[ip].jobz;
  if ((jobz != 'N') && (jobz != 'V')) {
    PRINTF("jobz should be N or V. Please correct the input data.");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangles of A and B are stored;
          = 'L':  Lower triangles of A and B are stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A and B.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  // LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
  integer lda = eig_paramslist[ip].lda;
  if (lda < max(1, n)) {
    PRINTF("lda < max(1, n) but it should be: LDA >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA, N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // LDB is INTEGER. The leading dimension of the array B.  LDB >= max(1,N).
  integer ldb = eig_paramslist[ip].ldb;
  if (ldb < max(1, n)) {
    PRINTF("ldb < max(1, n) but it should be: LDB >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // B is COMPLEX or COMPLEX*16 array, dimension (LDB, N)
  T *bbuff, *brefbuff;
  allocate_init_buffer(bbuff, brefbuff, ldb * n);
  
  // W is REAL or DOUBLE PRECISION array, dimension (N)
  Ta *wbuff, *wrefbuff;
  allocate_init_buffer(wbuff, wrefbuff, n, 0);
  
  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = eig_paramslist[ip].lwork_heev_2stage;
  integer lwork_size = lwork;
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  if ((n <= 1) && (lwork < 1)) {
    PRINTF("lwork should be >= 1 when n <= 1.\n");
  }
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (max(1, 3*N-2))
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, max(1, 3*n - 2), 0);
  
  integer info_cpp = -1;
  
  /* If lwork is -1, then call this API with default work variables.
     In return, work variables will be updated with array sizes needed.*/
  if (lwork == -1) {
    PRINTF("lwork/lrwork/liwork is -1, so call hegv_2stage() to get the" \
           "array sizes.\n");
    T worksize = {0};
    
    libflame::hegv_2stage<T, Ta>(&itype, &jobz, &uplo, &n, abuff, &lda, bbuff,
                &ldb, wbuff, &worksize, &lwork_size, rworkbuff, &info_cpp);
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
    PRINTF("itype = %d\n", itype);
    PRINTF("jobz = %c\n", jobz);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("ldb = %d\n", ldb);
    PRINTF("Size of B array (ldb*n) = %d\n", ldb * n);
    PRINTF("Size of W array (n) = %d\n", n);
    PRINTF("Size of WORK array (MAX(1, LWORK)) = %d\n", max(1, lwork_size));
    PRINTF("LWORK = %d\n", lwork_size);
    PRINTF("Size of RWORK array (7*n) = %d\n", max(1, 3*n-2));
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
    
    // Prints B array contents
    strncpy(arrayname, "B input", arraysize);
    print_array<T>(arrayname, bbuff, ldb * n);
    strncpy(arrayname, "B ref input", arraysize);
    print_array<T>(arrayname, brefbuff, ldb * n);
    
    // Prints W array contents
    strncpy(arrayname, "W input", arraysize);
    print_array<Ta>(arrayname, wbuff, n);
    strncpy(arrayname, "W ref input", arraysize);
    print_array<Ta>(arrayname, wrefbuff, n);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, max(1, lwork_size));
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
    
    // Prints RWORK array contents
    strncpy(arrayname, "RWORK input", arraysize);
    print_array<Ta>(arrayname, rworkbuff, max(1, 3*n-2));
    strncpy(arrayname, "RWORK ref input", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, max(1, 3*n-2));
  #endif
  
  info_cpp = -1;
  // Call CPP function
  libflame::hegv_2stage(&itype, &jobz, &uplo, &n, abuff, &lda, bbuff, &ldb,
              wbuff, workbuff, &lwork_size, rworkbuff, &info_cpp);

  // Call C function - NetLib Lapack API
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    hegv_2stage_ref = (fptr_NL_LAPACK_hegv_2stage)dlsym(lapackModule, \
                      "chegv_2stage_");
  } else if (typeid(T) == typeid(dcomplex)) {
    hegv_2stage_ref = (fptr_NL_LAPACK_hegv_2stage)dlsym(lapackModule, \
                      "zhegv_2stage_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (hegv_2stage_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  integer info_ref = -1;
  
  hegv_2stage_ref(&itype, &jobz, &uplo, &n, arefbuff, &lda, brefbuff, &ldb,
              wrefbuff, workrefbuff, &lwork_size, rworkrefbuff, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
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
      
      // Prints B array contents
      strncpy(arrayname, "B output", arraysize);
      print_array<T>(arrayname, bbuff, ldb * n);
      strncpy(arrayname, "B ref output", arraysize);
      print_array<T>(arrayname, brefbuff, ldb * n);
      
      // Prints W array contents
      strncpy(arrayname, "W output", arraysize);
      print_array<Ta>(arrayname, wbuff, n);
      strncpy(arrayname, "W ref output", arraysize);
      print_array<Ta>(arrayname, wrefbuff, n);
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, max(1, lwork_size));
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
      
      // Prints RWORK array contents
      strncpy(arrayname, "RWORK output", arraysize);
      print_array<Ta>(arrayname, rworkbuff, max(1, 3*n-2));
      strncpy(arrayname, "RWORK ref output", arraysize);
      print_array<Ta>(arrayname, rworkrefbuff, max(1, 3*n-2));
    #endif
    
    double diff = computeError<T>(lda, n, arefbuff, abuff);
    diff += computeError<T>(ldb, n, brefbuff, bbuff);
    diff += computeError<Ta>(1, n, wbuff, wrefbuff);
    if (lwork_size != 0) {
      diff += computeError<T>(1, lwork_size, workbuff, workrefbuff);
    }
    diff += computeError<Ta>(1, max(1, 3*n-2), rworkbuff, rworkrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, abs(diff), SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
           " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  if ((info_cpp > 0) || (info_ref > 0)) {
    if ((info_cpp <= n) || (info_ref <= n)) {
      PRINTF("if INFO = i, CHEEV failed to converge; " \
             "i off-diagonal elements of an intermediate " \
             "tridiagonal form did not converge to zero;\n");
    }
    else if ((info_cpp > n) || (info_ref > n)) {
      PRINTF("if INFO = N + i, for 1 <= i <= N, then the leading" \
             " minor of order i of B is not positive definite." \
             " The factorization of B could not be completed and" \
             " no eigenvalues or eigenvectors were computed.\n");
    }
  }
  
  //Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] bbuff; delete[] brefbuff;
  delete[] wbuff; delete[] wrefbuff;
  if (lwork_size != 0) {
    delete[] workbuff; delete[] workrefbuff; 
  }
  delete[] rworkbuff; delete[] rworkrefbuff;
}

TEST(LAPACKCPP_hegv_2stage, CHEGV_2STAGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hegv_2stage_test<scomplex, float> (0);
  }
}

TEST(LAPACKCPP_hegv_2stage, ZHEGV_2STAGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hegv_2stage_test<dcomplex, double> (0);
  }
}