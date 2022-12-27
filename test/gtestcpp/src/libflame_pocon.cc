/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_pocon.cc
 *  @brief Test application to validate pocon() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  pocon_test is function template for pocon() functions.
			T can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  pocon_test is function template for pocon() functions.
	  T can be float, double.
	  
	  pocon_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output if INFO == 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO != 0.
	  
    Real reference:
	  http://www.netlib.org/lapack/explore-html/d8/db2/group__real_p_ocomputational_gaca094dd6ef3db9ecb580ea731ecb5365.html#gaca094dd6ef3db9ecb580ea731ecb5365
	  Double reference:
	  http://www.netlib.org/lapack/explore-html/d1/d7a/group__double_p_ocomputational_ga42c90b8fcfef1a8f7c87a45e8176d643.html#ga42c90b8fcfef1a8f7c87a45e8176d643
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void pocon_test(int ip) {
  typedef integer (*fptr_NL_LAPACKE_pocon)(char* uplo, integer* n, T* a,
                      integer* lda, T* anorm, T* rcond, T* work,
                      integer* iwork, integer* info);
  fptr_NL_LAPACKE_pocon pocon_ref;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = lin_solver_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = lin_solver_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= fla_max(1,N).*/
  integer lda = lin_solver_paramslist[ip].lda;
  if (lda < fla_max(1, n)) {
    PRINTF("lda < fla_max(1, n) but it should be: LDA >= fla_max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is REAL or DOUBLE PRECISION array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // Call POTRF() to get A buffer and pass the buffer to POCON()
  PRINTF("Calling POTRF() to get A buffer and pass the buffer to POCON()\n");
  
  // Call potrf_internal() to get buffers.
  integer info_cpp = -1;
  integer info_ref = -1;
  potrf_internal<T>(uplo, n, abuff, lda, arefbuff, &info_cpp, &info_ref);
  PRINTF ("potrf_internal() info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  if ((info_cpp < 0) && (info_ref < 0)) {
    PRINTF("Info returned by CPP or C API is not successful to get" \
           " A buffer from POTRF(). Exiting....\n");
    closelibs();
    exit (-1);
  }
  
  // ANORM is REAL or DOUBLE PRECISION. The 1-norm of the original matrix A.
  T anorm =  (T)lin_solver_paramslist[ip].anorm;
  
  // RCOND is REAL or DOUBLE PRECISION.
  T rcond = 0.0, rcondref = 0.0;
  
  // WORK is REAL or DOUBLE PRECISION  array, dimension (3*N))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, 3*n, 0);
  
  // IWORK is INTEGER array, dimension (N)
  integer *iworkbuff = NULL, *iworkrefbuff = NULL;
  allocate_init_buffer(iworkbuff, iworkrefbuff, n, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("anorm = %lf\n", anorm);
    PRINTF("Size of WORK array (3*n) = %d\n", 3*n);
    PRINTF("Size of IWORK array (n) = %d\n", n);
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
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, 3 * n);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, 3 * n);
    
    // Prints IWORK array contents
    strncpy(arrayname, "IWORK input", arraysize);
    print_array<integer>(arrayname, iworkbuff, n);
    strncpy(arrayname, "IWORK ref input", arraysize);
    print_array<integer>(arrayname, iworkrefbuff, n);
  #endif
  
  info_cpp = -1;
  // Call CPP function
  libflame::pocon<T>(&uplo, &n, abuff, &lda, &anorm, &rcond, workbuff,
                    iworkbuff, &info_cpp);
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(float)) {
    pocon_ref = (fptr_NL_LAPACKE_pocon)dlsym(lapackModule, "spocon_");
  } else if (typeid(T) == typeid(double)) {
    pocon_ref = (fptr_NL_LAPACKE_pocon)dlsym(lapackModule, "dpocon_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (pocon_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  info_ref = -1;
  pocon_ref(&uplo, &n, arefbuff, &lda, &anorm, &rcondref, workrefbuff,
        iworkrefbuff, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  PRINTF ("rcond: %lf, rcondref: %lf\n", rcond, rcondref);
  
  // Calculate the differences of buffers.
  if ((info_cpp == 0) && (info_ref == 0)) {
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
      // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
      PRINTF("\nPrinting all Output arrays contents...\n");
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, 3 * n);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, 3 * n);
      
      // Prints IWORK array contents
      strncpy(arrayname, "IWORK output", arraysize);
      print_array<integer>(arrayname, iworkbuff, n);
      strncpy(arrayname, "IWORK ref output", arraysize);
      print_array<integer>(arrayname, iworkrefbuff, n);
    #endif
    double diff = computeError<T>(1, 1, &rcond, &rcondref);
    diff += computeError<T>(3, n, workbuff, workrefbuff);
    diff += computeError<integer>(1, n, iworkbuff, iworkrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp != 0) || (info_ref != 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] iworkbuff; delete[] iworkrefbuff;
}

/*! @brief  pocon_test_cmplx is function template for pocon() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  pocon_test_cmplx is function template for pocon() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  pocon_test_cmplx() function template calls C and CPP based lbrary APIs
	  with valid test values and returns the differences in output if INFO == 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO != 0.
	  
    Complex reference:
	  http://www.netlib.org/lapack/explore-html/d6/df6/group__complex_p_ocomputational_ga2ddc05543f7ed596609cdce0478ca8a3.html#ga2ddc05543f7ed596609cdce0478ca8a3
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d3/d8d/group__complex16_p_ocomputational_gaa3938ab5d7bc02f1d7115794d242b7d0.html#gaa3938ab5d7bc02f1d7115794d242b7d0
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template<typename T, typename Ta>
void pocon_test_cmplx(int ip)
{
  typedef integer (*fptr_NL_LAPACKE_pocon)(char* uplo, integer* n, T* a,
                      integer* lda, Ta* anorm, Ta* rcond, T* work,
                      Ta* rwork, integer* info);
  fptr_NL_LAPACKE_pocon pocon_ref;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = lin_solver_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }

  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = lin_solver_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= fla_max(1,N).*/
  integer lda = lin_solver_paramslist[ip].lda;
  if (lda < fla_max(1, n)) {
    PRINTF("lda < fla_max(1, n) but it should be: LDA >= fla_max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // Call POTRF() to get A buffer and pass the buffer to POCON()
  PRINTF("Calling POTRF() to get A buffer and pass the buffer to POCON()\n");
  
  // Call potrf_internal() to get buffers.
  integer info_cpp = -1;
  integer info_ref = -1;
  potrf_internal<T>(uplo, n, abuff, lda, arefbuff, &info_cpp, &info_ref);
  PRINTF ("potrf_internal() info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  if ((info_cpp < 0) && (info_ref < 0)) {
    PRINTF("Info returned by CPP or C API is not successful to get" \
           " A buffer from POTRF(). Exiting....\n");
    closelibs();
    exit (-1);
  }
  
  // ANORM is REAL or DOUBLE PRECISION. The 1-norm of the original matrix A.
  Ta anorm =  (Ta)lin_solver_paramslist[ip].anorm;
  
  // RCOND is REAL or DOUBLE PRECISION.
  Ta rcond = 0.0, rcondref = 0.0;
  
  // WORK is COMPLEX or COMPLEX*16 array, dimension (2*N)
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, 2*n, 0);
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (N)
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, n, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("anorm = %lf\n", anorm);
    PRINTF("Size of WORK array (2*n) = %d\n", 2*n);
    PRINTF("Size of RWORK array (n) = %d\n", n);
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
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, 2 * n);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, 2 * n);
    
    // Prints RWORK array contents
    strncpy(arrayname, "RWORK input", arraysize);
    print_array<Ta>(arrayname, rworkbuff, n);
    strncpy(arrayname, "RWORK ref input", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, n);
  #endif
  
  info_cpp = -1;
  // Call CPP function
  libflame::pocon<T, Ta>(&uplo, &n, abuff, &lda, &anorm, &rcond, workbuff,
                rworkbuff, &info_cpp);
  
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    pocon_ref = (fptr_NL_LAPACKE_pocon)dlsym(lapackModule, "cpocon_");
  } else if (typeid(T) == typeid(dcomplex)) {
    pocon_ref = (fptr_NL_LAPACKE_pocon)dlsym(lapackModule, "zpocon_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (pocon_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  info_ref = -1;
  pocon_ref(&uplo, &n, arefbuff, &lda, &anorm, &rcondref, workrefbuff,
        rworkrefbuff, &info_ref);
  PRINTF ("info_cpp: %u, info_ref: %u\n", info_cpp, info_ref);
  PRINTF ("rcond: %lf, rcondref: %lf\n", rcond, rcondref);
  
  // Calculate the differences of buffers.
  if ((info_cpp == 0) && (info_ref == 0)) {
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
      // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
      PRINTF("\nPrinting all Output arrays contents...\n");
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, 2 * n);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, 2 * n);
      
      // Prints RWORK array contents
      strncpy(arrayname, "RWORK output", arraysize);
      print_array<Ta>(arrayname, rworkbuff, n);
      strncpy(arrayname, "RWORK ref output", arraysize);
      print_array<Ta>(arrayname, rworkrefbuff, n);
    #endif
    double diff = computeError<Ta>(1, 1, &rcond, &rcondref);
    diff += computeError<T>(2, n, workbuff, workrefbuff);
    diff += computeError<Ta>(1, n, rworkbuff, rworkrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp != 0) || (info_ref != 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
}

/* Use TEST macro and call C++ test function template with
   float as typenames.*/
TEST(LAPACKCPP_pocon, SPOCON) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    pocon_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double as typenames.*/
TEST(LAPACKCPP_pocon, DPOCON) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    pocon_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_pocon, CPOCON) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    pocon_test_cmplx<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_pocon, ZPOCON) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    pocon_test_cmplx<dcomplex, double> (index);
  }
}