/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_geqp3.cc
 *  @brief Test application to validate geqp3() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  geqp3_test is function template for geqp3() functions.
			T can be float, double
 * @details
 * \b Purpose:
    \verbatim
	  geqp3_test is function template for geqp3() functions.
	  T can be float, double
	  
	  geqp3_test() function template calls C and CPP based library APIs with
	  valid test values, calculate the differences in output if INFO == 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

    Real Reference:
    http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga63f9e3af96fa42609e41bf3d77660bdf.html#ga63f9e3af96fa42609e41bf3d77660bdf
    Double Reference:
    https://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga1b0500f49e03d2771b797c6e88adabbb.html#ga1b0500f49e03d2771b797c6e88adabbb
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void geqp3_test(int ip)
{
  typedef integer (*fptr_NL_LAPACK_geqp3)(integer* m, integer* n, T* a,
                      integer* lda, integer* jpvt, T* tau, T* work,
                      integer* lwork, integer* info);
  fptr_NL_LAPACK_geqp3 geqp3_ref = NULL;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* M is INTEGER
          The number of rows of the matrix A. M >= 0.*/
  integer m = eig_paramslist[ip].m;
  if (m < 0) {
    PRINTF("m < 0 but it should be: m >= 0. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(M,1).*/
  integer lda = eig_paramslist[ip].lda_lange;
  if (lda < max(1, m)) {
    PRINTF("lda < max(1, m) but it should be: LDA >= max(1,M). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // JPVT is INTEGER array, dimension (N)
  integer *jpvtbuff, *jpvtrefbuff;
  allocate_init_buffer(jpvtbuff, jpvtrefbuff, n);
  
  // TAU is REAL or DOUBLE PRECISION array, dimension (min(M,N))
  T *taubuff, *taurefbuff;
  allocate_init_buffer(taubuff, taurefbuff, min(m, n), 0);
  
  /* LWORK is INTEGER
          The dimension of the array WORK. LWORK >= 3*N+1.
          For optimal performance LWORK >= 2*N+( N+1 )*NB, where NB
          is the optimal blocksize.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.*/
  integer lwork = eig_paramslist[ip].lwork_geqp3;
  integer lwork_size = lwork;
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  integer info_cpp = -1;
  
  /* If lwork is -1, then call this API with default work variables.
     In return, work variable will be updated with array sizes needed.*/
  if (lwork == -1) {
    PRINTF("lwork is -1, so call geqp3() to get the array sizes.\n");
    T worksize = {0};
    libflame::geqp3<T>(&m, &n, abuff, &lda, jpvtbuff, taubuff, &worksize,
                &lwork_size, &info_cpp);
    PRINTF("info_cpp: %d\n", info_cpp);
    if (info_cpp == 0) {
      lwork_size = worksize;
    }
  }
  
  // WORK is REAL or DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, max(1, lwork_size), 0);
 
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("m = %d\n", m);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of JPVT array (n) = %d\n", n);
    PRINTF("Size of TAU array (min(M,N)) = %d\n", min(m, n));
    PRINTF("lwork = %d\n", lwork_size);
    PRINTF("Size of WORK array (MAX(1,LWORK)) = %d\n", max(1, lwork_size));
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
    
    // Prints JPVT array contents
    strncpy(arrayname, "JPVT input", arraysize);
    print_array<integer>(arrayname, jpvtbuff, n);
    strncpy(arrayname, "JPVT ref input", arraysize);
    print_array<integer>(arrayname, jpvtrefbuff, n);
    
    // Prints TAU array contents
    strncpy(arrayname, "TAU input", arraysize);
    print_array<T>(arrayname, taubuff, min(m, n));
    strncpy(arrayname, "TAU ref input", arraysize);
    print_array<T>(arrayname, taurefbuff, min(m, n));
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, max(1, lwork_size));
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
  #endif
  
  // Call CPP function
  libflame::geqp3<T>(&m, &n, abuff, &lda, jpvtbuff, taubuff, workbuff,
              &lwork_size, &info_cpp);
  
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(float)) {
    geqp3_ref = (fptr_NL_LAPACK_geqp3)dlsym(lapackModule, "sgeqp3_");
  } else if (typeid(T) == typeid(double)) {
    geqp3_ref = (fptr_NL_LAPACK_geqp3)dlsym(lapackModule, "dgeqp3_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (geqp3_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  integer info_ref = -1;
  geqp3_ref(&m, &n, arefbuff, &lda, jpvtrefbuff, taurefbuff, workrefbuff,
      &lwork_size, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all Output arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A output", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref output", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
    
    // Prints JPVT array contents
    strncpy(arrayname, "JPVT output", arraysize);
    print_array<integer>(arrayname, jpvtbuff, n);
    strncpy(arrayname, "JPVT ref output", arraysize);
    print_array<integer>(arrayname, jpvtrefbuff, n);
    
    // Prints TAU array contents
    strncpy(arrayname, "TAU output", arraysize);
    print_array<T>(arrayname, taubuff, min(m, n));
    strncpy(arrayname, "TAU ref output", arraysize);
    print_array<T>(arrayname, taurefbuff, min(m, n));
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK output", arraysize);
    print_array<T>(arrayname, workbuff, max(1, lwork_size));
    strncpy(arrayname, "WORK ref output", arraysize);
    print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
  #endif
  
  if ((info_cpp == 0) && (info_ref == 0)) {
    double diff = 0.0;
    // TODO: Yet to finalize and do verification changes.
    /*diff = computeError<T>(lda, n, arefbuff, abuff);
    diff += computeError<integer>(1, n, jpvtrefbuff, jpvtbuff);
    diff += computeError<T>(1, min(m, n), taurefbuff, taubuff);
    diff += computeError<T>(1, max(1, lwork_size), workbuff, workrefbuff);*/
    
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] jpvtbuff; delete[] jpvtrefbuff;
  delete[] taubuff; delete[] taurefbuff;
  delete[] workbuff; delete[] workrefbuff;
}

/*! @brief  geqp3_test is function template for geqp3() functions.
			T can be scomplex, dcomplex
      Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  geqp3_test is function template for geqp3() functions.
	  T can be scomplex, dcomplex
    Ta can be float, double.
	  
	  geqp3_test() function template calls C and CPP based library APIs with
	  valid test values, calculate the differences in output if INFO == 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.
    
    Complex reference:
	  http://www.netlib.org/lapack/explore-html/d4/d7e/group__complex_g_ecomputational_ga3947eb2e884bf32f7380f22c501151e9.html#ga3947eb2e884bf32f7380f22c501151e9
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d0/d9e/group__complex16_g_eauxiliary_ga7908bb12a6f02dbfa4d5a92a27c0e9b7.html#ga7908bb12a6f02dbfa4d5a92a27c0e9b7
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T, typename Ta>
void geqp3_test_cmplx(int ip)
{
  typedef integer (*fptr_NL_LAPACK_geqp3)(integer* m, integer* n, T* a,
                      integer* lda, integer* jpvt, T* tau, T* work,
                      integer* lwork, Ta* rwork, integer* info);
  fptr_NL_LAPACK_geqp3 geqp3_ref = NULL;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* M is INTEGER
          The number of rows of the matrix A. M >= 0.*/
  integer m = eig_paramslist[ip].m;
  if (m < 0) {
    PRINTF("m < 0 but it should be: m >= 0. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(M,1).*/
  integer lda = eig_paramslist[ip].lda_lange;
  if (lda < max(1, m)) {
    PRINTF("lda < max(1, m) but it should be: LDA >= max(1,M). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // JPVT is INTEGER array, dimension (N)
  integer *jpvtbuff, *jpvtrefbuff;
  allocate_init_buffer(jpvtbuff, jpvtrefbuff, n);
  
  // TAU is REAL or DOUBLE PRECISION array, dimension (min(M,N))
  T *taubuff, *taurefbuff;
  allocate_init_buffer(taubuff, taurefbuff, min(m, n), 0);
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (2*N)
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, 2*n, 0);
  
  /* LWORK is INTEGER
          The dimension of the array WORK. LWORK >= 3*N+1.
          For optimal performance LWORK >= 2*N+( N+1 )*NB, where NB
          is the optimal blocksize.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.*/
  integer lwork = eig_paramslist[ip].lwork_geqp3;
  integer lwork_size = lwork;
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  integer info_cpp = -1;
  
  /* If lwork is -1, then call this API with default work variables.
     In return, work variable will be updated with array sizes needed.*/
  if (lwork == -1) {
    PRINTF("lwork is -1, so call geqp3() to get the array sizes.\n");
    T worksize = {0};
    libflame::geqp3<T, Ta>(&m, &n, abuff, &lda, jpvtbuff, taubuff, &worksize,
                &lwork_size, rworkbuff, &info_cpp);
    PRINTF("info_cpp: %d\n", info_cpp);
    if (info_cpp == 0) {
        lwork_size = worksize.real;
    }
  }
  
  // WORK is COMPLEX or COMPLEX*16 array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, max(1, lwork_size), 0);

  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("m = %d\n", m);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of JPVT array (n) = %d\n", n);
    PRINTF("Size of TAU array (min(M,N)) = %d\n", min(m, n));
    PRINTF("lwork = %d\n", lwork_size);
    PRINTF("Size of WORK array (MAX(1,LWORK)) = %d\n", max(1, lwork_size));
    PRINTF("Size of RWORK array (2*n) = %d\n", 2*n);
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
    
    // Prints JPVT array contents
    strncpy(arrayname, "JPVT input", arraysize);
    print_array<integer>(arrayname, jpvtbuff, n);
    strncpy(arrayname, "JPVT ref input", arraysize);
    print_array<integer>(arrayname, jpvtrefbuff, n);
    
    // Prints TAU array contents
    strncpy(arrayname, "TAU input", arraysize);
    print_array<T>(arrayname, taubuff, min(m, n));
    strncpy(arrayname, "TAU ref input", arraysize);
    print_array<T>(arrayname, taurefbuff, min(m, n));
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, max(1, lwork_size));
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
    
    // Prints RWORK array contents
    strncpy(arrayname, "RWORK input", arraysize);
    print_array<Ta>(arrayname, rworkbuff, 2*n);
    strncpy(arrayname, "RWORK ref input", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, 2*n);
  #endif
  
  // Call CPP function
  libflame::geqp3<T, Ta>(&m, &n, abuff, &lda, jpvtbuff, taubuff, workbuff,
              &lwork_size, rworkbuff, &info_cpp);
  
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    geqp3_ref = (fptr_NL_LAPACK_geqp3)dlsym(lapackModule, "cgeqp3_");
  } else if (typeid(T) == typeid(dcomplex)) {
    geqp3_ref = (fptr_NL_LAPACK_geqp3)dlsym(lapackModule, "zgeqp3_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (geqp3_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  integer info_ref = -1;
  geqp3_ref(&m, &n, arefbuff, &lda, jpvtrefbuff, taurefbuff, workrefbuff,
      &lwork_size, rworkrefbuff, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all Output arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A output", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref output", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
    
    // Prints JPVT array contents
    strncpy(arrayname, "JPVT output", arraysize);
    print_array<integer>(arrayname, jpvtbuff, n);
    strncpy(arrayname, "JPVT ref output", arraysize);
    print_array<integer>(arrayname, jpvtrefbuff, n);
    
    // Prints TAU array contents
    strncpy(arrayname, "TAU output", arraysize);
    print_array<T>(arrayname, taubuff, min(m, n));
    strncpy(arrayname, "TAU ref output", arraysize);
    print_array<T>(arrayname, taurefbuff, min(m, n));
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK output", arraysize);
    print_array<T>(arrayname, workbuff, max(1, lwork_size));
    strncpy(arrayname, "WORK ref output", arraysize);
    print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
    
    // Prints RWORK array contents
    strncpy(arrayname, "RWORK output", arraysize);
    print_array<Ta>(arrayname, rworkbuff, 2*n);
    strncpy(arrayname, "RWORK ref output", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, 2*n);
  #endif
  
  if ((info_cpp == 0) && (info_ref == 0)) {
    double diff = 0.0;
    // TODO: Yet to finalize and do verification changes.
    /*diff = computeError<T>(lda, n, arefbuff, abuff);
    diff += computeError<integer>(1, n, jpvtrefbuff, jpvtbuff);
    diff += computeError<T>(1, min(m, n), taurefbuff, taubuff);
    diff += computeError<T>(1, max(1, lwork_size), workbuff, workrefbuff);
    diff += computeError<Ta>(1, 2*n, rworkbuff, rworkrefbuff);*/
    
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] jpvtbuff; delete[] jpvtrefbuff;
  delete[] taubuff; delete[] taurefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
}

/* Use TEST macro and call C++ test function template with
   float as typenames.*/
TEST(LAPACKCPP_geqp3, SGEQP3) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    geqp3_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double as typenames.*/
TEST(LAPACKCPP_geqp3, DGEQP3) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    geqp3_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_geqp3, CGEQP3) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    geqp3_test_cmplx<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_geqp3, ZGEQP3) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    geqp3_test_cmplx<dcomplex, double> (index);
  }
}