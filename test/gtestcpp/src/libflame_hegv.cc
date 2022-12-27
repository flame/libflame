/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hegv.cc
 *  @brief Test application to validate hegv() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hegv_test is function template for hegv() functions.
			T can be scomplex, dcomplex
      Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hegv_test is function template for hegv() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  hegv_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d9/de3/group__complex_h_eeigen_gab2f86fb41df5ae239798c9c3081a2d49.html#gab2f86fb41df5ae239798c9c3081a2d49
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaf7b790b3b89de432a423c9006c1cc1ac.html#gaf7b790b3b89de432a423c9006c1cc1ac
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T, typename Ta>
void hegv_test(int ip)
{
  typedef integer (*fptr_NL_LAPACK_hegv)(integer* itype, char* jobz, char* uplo,
                  integer* n, T* a, integer* lda, T* b, integer* ldb, Ta* w,
                  T* work, integer* lwork, Ta* rwork, integer* info);
  fptr_NL_LAPACK_hegv hegv_ref = NULL;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* ITYPE is INTEGER
          Specifies the problem type to be solved:
          = 1:  A*x = (lambda)*B*x
          = 2:  A*B*x = (lambda)*x
          = 3:  B*A*x = (lambda)*x*/
  integer itype = eig_paramslist[ip].itype;
  if ((itype < 1) || (itype > 3)) {
    PRINTF("itype should be from 1 to 3. Please correct the input data.\n");
  }
  
  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.*/
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
  
  // LDA is INTEGER. The leading dimension of the array A.  LDA >= fla_max(1,N).
  integer lda = eig_paramslist[ip].lda;
  if (lda < fla_max(1, n)) {
    PRINTF("lda < fla_max(1, n) but it should be: LDA >= fla_max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA, N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // LDB is INTEGER. The leading dimension of the array B.  LDB >= fla_max(1,N).
  integer ldb = eig_paramslist[ip].ldb;
  if (ldb < fla_max(1, n)) {
    PRINTF("ldb < fla_max(1, n) but it should be: LDB >= fla_max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // B is COMPLEX or COMPLEX*16 array, dimension (LDB, N)
  T *bbuff, *brefbuff;
  allocate_init_buffer(bbuff, brefbuff, ldb * n);
  
  // W is REAL or DOUBLE PRECISION array, dimension (N)
  Ta *wbuff, *wrefbuff;
  allocate_init_buffer(wbuff, wrefbuff, n, 0);
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (max(1,3*N-2))
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, fla_max(1, 3*n-2), 0);
  
  //  LWORK is INTEGER. The length of the array WORK.
  integer lwork = eig_paramslist[ip].lwork_hegv;
  integer lwork_size = lwork;
  integer info_cpp = -1;
  
  /* If lwork is -1, then call this API with default work variable.
     In return, work variable will be updated with lwork needed.*/
  if (lwork == -1) {
    PRINTF("lwork is -1, so call hegv() to get the WORK array size.\n");
    T worksize = {0};
    
    // Call CPP function
    libflame::hegv<T, Ta>(&itype, &jobz, &uplo, &n, abuff, &lda, bbuff, &ldb,
                wbuff, &worksize, &lwork_size, rworkbuff, &info_cpp);
    PRINTF("info_cpp: %d, worksize: %f\n", info_cpp, worksize.real);
    if (info_cpp == 0) {
      lwork_size = worksize.real;
    }
  } else if (lwork < fla_max(1, 2*n-1)) {
    PRINTF("lwork should be >= fla_max(1,2*n-1). Please correct the input data.\n");
  }
  
  // WORK is COMPLEX or COMPLEX*16 array, dimension (MAX(1, LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, fla_max(1, lwork_size), 0);
  
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
    PRINTF("Size of WORK array (MAX(1, LWORK)) = %d\n", fla_max(1, lwork_size));
    PRINTF("LWORK = %d\n", lwork_size);
    PRINTF("Size of RWORK array (7*n) = %d\n", fla_max(1, 3*n-2));
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
    print_array<T>(arrayname, workbuff, fla_max(1, lwork_size));
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, fla_max(1, lwork_size));
    
    // Prints RWORK array contents
    strncpy(arrayname, "RWORK input", arraysize);
    print_array<Ta>(arrayname, rworkbuff, fla_max(1, 3*n-2));
    strncpy(arrayname, "RWORK ref input", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, fla_max(1, 3*n-2));
  #endif
  
  info_cpp = -1;
  // Call CPP function
  libflame::hegv<T, Ta>(&itype, &jobz, &uplo, &n, abuff, &lda, bbuff, &ldb,
              wbuff, workbuff, &lwork_size, rworkbuff, &info_cpp);

  // Call C function - NetLib Lapack API
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    hegv_ref = (fptr_NL_LAPACK_hegv)dlsym(lapackModule, "chegv_");
  } else if (typeid(T) == typeid(dcomplex)) {
    hegv_ref = (fptr_NL_LAPACK_hegv)dlsym(lapackModule, "zhegv_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (hegv_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  integer info_ref = -1;
  hegv_ref(&itype, &jobz, &uplo, &n, arefbuff, &lda, brefbuff, &ldb, wrefbuff,
      workrefbuff, &lwork_size, rworkrefbuff, &info_ref);
  PRINTF("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
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
      print_array<T>(arrayname, workbuff, fla_max(1, lwork_size));
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, fla_max(1, lwork_size));
      
      // Prints RWORK array contents
      strncpy(arrayname, "RWORK output", arraysize);
      print_array<Ta>(arrayname, rworkbuff, fla_max(1, 3*n-2));
      strncpy(arrayname, "RWORK ref output", arraysize);
      print_array<Ta>(arrayname, rworkrefbuff, fla_max(1, 3*n-2));
    #endif
    
    double diff = computeError<T>(lda, n, arefbuff, abuff);
    diff += computeError<T>(ldb, n, brefbuff, bbuff);
    diff += computeError<Ta>(1, n, wbuff, wrefbuff);
    if (lwork_size != 0) {
      diff += computeError<T>(1, lwork_size, workbuff, workrefbuff);
    }
    diff += computeError<Ta>(1, fla_max(1, 3*n-2), rworkbuff, rworkrefbuff);
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
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] bbuff; delete[] brefbuff;
  delete[] wbuff; delete[] wrefbuff;
  if (lwork_size != 0) {
    delete[] workbuff; delete[] workrefbuff; 
  }
  delete[] rworkbuff; delete[] rworkrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_hegv, CHEGV) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hegv_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_hegv, ZHEGV) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hegv_test<dcomplex, double> (index);
  }
}