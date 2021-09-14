/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_heevd_2stage.cc
 *  @brief Test application to validate heevd_2stage() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  heevd_2stage_test is function template for heevd_2stage() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  heevd_2stage_test is function template for heevd_2stage() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  heevd_2stage_test() function template calls C and CPP based lbrary APIs
	  with valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.
	  
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d9/de3/group__complex_h_eeigen_ga7d77c89933d557e7131ebb6474e2a3c0.html#ga7d77c89933d557e7131ebb6474e2a3c0
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_ga59996f9a96d264d7eba3edef73e32c43.html#ga59996f9a96d264d7eba3edef73e32c43
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template<typename T, typename Ta>
void heevd_2stage_test(int ip)
{
  typedef int (*fptr_NL_LAPACK_heevd_2stage)(char* jobz, char* uplo, integer* n,
                  T* a, integer* lda, Ta* w, T* work, integer* lwork,
                  Ta* rwork, integer* lrwork, integer* iwork, integer* liwork,
                  integer* info);
  fptr_NL_LAPACK_heevd_2stage heevd_2stage_ref;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
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

  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = eig_paramslist[ip].lwork_heevd_2stage;
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
  
  if ((lwork != -1) && (jobz == 'N') && (n > 1)) {
    PRINTF("If JOBZ = 'N' and N > 1, lwork must be queried.\n");
    lwork = -1;
  }
  
  /* If lwork/lrwork/liwork is -1, then call this API with default work
     variables. In return, work variables will be updated with array sizes
     needed.*/
  if ((lwork == -1) || (lrwork == -1) || (liwork == -1)) {
    PRINTF("lwork/lrwork/liwork is -1, so call heevd_2stage() to get the " \
           "array sizes.\n");
    T worksize = {0};
    Ta rworksize = 0.0;
    integer iworksize = 0;
    
    // Call CPP function
    libflame::heevd_2stage<T, Ta>(&jobz, &uplo, &n, abuff, &lda, wbuff,
                  &worksize, &lwork, &rworksize, &lrwork, &iworksize,
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
    } else if (jobz == 'V') {
      if (lwork < (2*n+n*2)) {
        PRINTF("lwork must be atleast 2*n+n*2. Please correct the input" \
               " data.\n");
      }
    }
    lwork_size = lwork;
  }

  // WORK is COMPLEX or COMPLEX*16  array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, max(1, lwork_size), 0);
  
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
  allocate_init_buffer(rworkbuff, rworkrefbuff, lrwork_size, 0);
  
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
  allocate_init_buffer(iworkbuff, iworkrefbuff, max(1, liwork_size), 0);
  
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
    PRINTF("lrwork = %d\n", lrwork_size);
    PRINTF("Size of RWORK array (lrwork) = %d\n", lrwork_size);
    PRINTF("liwork = %d\n", liwork_size);
    PRINTF("Size of IWORK array (max(1, liwork)) = %d\n", max(1, liwork_size));
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
    print_array<Ta>(arrayname, rworkbuff, lrwork_size);
    strncpy(arrayname, "RWORK ref input", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, lrwork_size);
    
    // Prints IWORK array contents
    strncpy(arrayname, "IWORK input", arraysize);
    print_array<integer>(arrayname, iworkbuff, max(1, liwork_size));
    strncpy(arrayname, "IWORK ref input", arraysize);
    print_array<integer>(arrayname, iworkrefbuff, max(1, liwork_size));
  #endif
  
  info_cpp = -1;
  
  // Call CPP function
  libflame::heevd_2stage<T, Ta>(&jobz, &uplo, &n, abuff, &lda, wbuff,
                  workbuff, &lwork_size, rworkbuff, &lrwork_size, iworkbuff,
                  &liwork_size, &info_cpp);

  // Call C function
  
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    heevd_2stage_ref = (fptr_NL_LAPACK_heevd_2stage)dlsym(lapackModule, \
                        "cheevd_2stage_");
  } else if (typeid(T) == typeid(dcomplex)) {
    heevd_2stage_ref = (fptr_NL_LAPACK_heevd_2stage)dlsym(lapackModule, \
                        "zheevd_2stage_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (heevd_2stage_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  integer info_ref = -1;
  heevd_2stage_ref(&jobz, &uplo, &n, arefbuff, &lda, wrefbuff, workrefbuff,
        &lwork_size, rworkrefbuff, &lrwork_size, iworkrefbuff, &liwork_size,
        &info_ref);
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
      print_array<Ta>(arrayname, rworkbuff, lrwork_size);
      strncpy(arrayname, "RWORK ref output", arraysize);
      print_array<Ta>(arrayname, rworkrefbuff, lrwork_size);
      
      // Prints IWORK array contents
      strncpy(arrayname, "IWORK output", arraysize);
      print_array<integer>(arrayname, iworkbuff, max(1, liwork_size));
      strncpy(arrayname, "IWORK ref output", arraysize);
      print_array<integer>(arrayname, iworkrefbuff, max(1, liwork_size));
    #endif
    
    double diff = computeError<T>(lda, n, arefbuff, abuff);
	  diff += computeError<Ta>(1, n, wrefbuff, wbuff);
    diff += computeError<T>(1, max(1, lwork_size), workbuff, workrefbuff);
    diff += computeError<Ta>(1, lrwork_size, rworkbuff, rworkrefbuff);
    diff += computeError<integer>(1, max(1, liwork_size), iworkbuff,
                                  iworkrefbuff);
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
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
  delete[] iworkbuff; delete[] iworkrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_heevd_2stage, CHEEVD_2STAGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    heevd_2stage_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_heevd_2stage, ZHEEVD_2STAGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    heevd_2stage_test<dcomplex, double> (index);
  }
}