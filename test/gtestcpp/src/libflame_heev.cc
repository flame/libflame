/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_heev.cc
 *  @brief Test application to validate heev() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  heev_test is function template for heev() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  heev_test is function template for heev() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  heev_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.
	  
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d9/de3/group__complex_h_eeigen_ga003ee37091d65ee62fd72da1035f06e2.html#ga003ee37091d65ee62fd72da1035f06e2
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaf23fb5b3ae38072ef4890ba43d5cfea2.html#gaf23fb5b3ae38072ef4890ba43d5cfea2
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template<typename T, typename Ta>
void heev_test(int ip)
{
  typedef int (*Fptr_NL_LAPACK_heev)(char* jobz, char* uplo, integer* n, T* a,
                  integer* lda, Ta*  w, T* work, integer* lwork, Ta* rwork,
                  integer* info);
  Fptr_NL_LAPACK_heev HEEV;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.*/
  char jobz = eig_paramslist[ip].jobz;
  if ((jobz != 'N') && (jobz != 'V')) {
    PRINTF("jobz should be N or V. Please correct the input data.");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.");
  }
  
  /* N is INTEGER
          The order of the matrix A. N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).*/
  integer lda = eig_paramslist[ip].lda;
  if (lda < max(1, n)) {
    PRINTF("lda < max(1, n) but it should be: lda >= max(1, n). Please " \
           "correct the input data.");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);

  // W is REAL or DOUBLE PRECISION array, dimension (N).
  Ta *wbuff, *wrefbuff;
  allocate_init_buffer(wbuff, wrefbuff, n, 0);

  // RWORK is REAL or DOUBLE PRECISION array, dimension max(1,3*N-2)
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, max(1, 3*n-2), 0);
  
  //  LWORK is INTEGER. The length of the array WORK.
  integer lwork = eig_paramslist[ip].lwork;
  integer lwork_size = lwork;
  integer info_cpp = -1;
  
  /* If lwork is -1, then call this API with default work variable.
     In return, work variable will be updated with lwork needed.*/
  if (lwork == -1) {
    PRINTF("lwork is -1, so call heev() to get the WORK array size.\n");
    T worksize = {0};
    
    // Call CPP function
    libflame::heev<T, Ta>(&jobz, &uplo, &n, abuff, &lda, wbuff, &worksize,
                &lwork, rworkbuff, &info_cpp);
    PRINTF("info_cpp: %d, worksize: %f\n", info_cpp, worksize.real);
    if (info_cpp == 0) {
      lwork_size = worksize.real;
    }
  }
  
  // WORK is COMPLEX or COMPLEX*16 array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, max(1, lwork_size), 0);
  
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    // Print input values other than arrays.
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("jobz = %c\n", jobz);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", (lda*n));
    PRINTF("Size of W array (n) = %d\n", n);
    PRINTF("lwork = %d\n", lwork_size);
    PRINTF("Size of WORK array (max(1, lwork)) = %d\n", max(1, lwork_size));
    PRINTF("Size of RWORK array (max(1, 3*n-2)) = %d\n", max(1, 3*n-2));
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
  libflame::heev<T, Ta>(&jobz, &uplo, &n, abuff, &lda, wbuff, workbuff,
              &lwork_size, rworkbuff, &info_cpp);

  // Call C function
  
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    HEEV = (Fptr_NL_LAPACK_heev)dlsym(lapackModule, "cheev_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HEEV = (Fptr_NL_LAPACK_heev)dlsym(lapackModule, "zheev_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HEEV == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  integer info_ref = -1;
  HEEV(&jobz, &uplo, &n, arefbuff, &lda, wrefbuff, workrefbuff, &lwork_size,
              rworkrefbuff, &info_ref);
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
	  diff += computeError<Ta>(1, n, wrefbuff, wbuff);
    if (lwork_size != 0) {
      diff += computeError<T>(1, max(1, lwork_size), workbuff, workrefbuff);
    }
    diff += computeError<Ta>(1, max(1, 3*n-2), rworkbuff, rworkrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, abs(diff), SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
           " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] wbuff; delete[] wrefbuff;
  if (lwork_size != 0) {
    delete[] workbuff; delete[] workrefbuff; 
  }
  delete[] rworkbuff; delete[] rworkrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_heev, CHEEV) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    heev_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_heev, ZHEEV) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    heev_test<dcomplex, double> (index);
  }
}