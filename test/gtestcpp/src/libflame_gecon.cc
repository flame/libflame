/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_gecon.cc
 *  @brief Test application to validate gecon() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  gecon_test is function template for gecon() functions.
			T can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  gecon_test is function template for gecon() functions.
	  T can be float, double.
	  
	  gecon_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output.
	  
    Real reference:
	  http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga89f21d7700aaccc5fc72ca3316c33463.html#ga89f21d7700aaccc5fc72ca3316c33463
	  Double reference:
	  https://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga188b8d30443d14b1a3f7f8331d87ae60.html#ga188b8d30443d14b1a3f7f8331d87ae60
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void gecon_test(int ip) {
  typedef integer (*Fptr_NL_LAPACKE_gecon)(char* norm, integer* n, T* a,
                      integer* lda, T* anorm, T* rcond, T* work,
                      integer* iwork, integer* info);
  Fptr_NL_LAPACKE_gecon gecon_ref;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
  /* NORM is CHARACTER*1
          Specifies whether the 1-norm condition number or the
          infinity-norm condition number is required:
          = '1' or 'O':  1-norm;
          = 'I':         Infinity-norm.*/
  char norm = lin_solver_paramslist[ip].norm;
  if ((norm != '1') && (norm != 'O') && (norm != 'I')) {
    PRINTF("norm should be 1 or O or I. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = lin_solver_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).*/
  integer lda = lin_solver_paramslist[ip].lda;
  if (lda < max(1, n)) {
    PRINTF("lda < max(1, n) but it should be: LDA >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is REAL or DOUBLE PRECISION array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // ANORM is REAL or DOUBLE PRECISION. The 1-norm of the original matrix A.
  T anorm =  (T)lin_solver_paramslist[ip].anorm;
  
  // RCOND is REAL or DOUBLE PRECISION.
  T rcond = 0.0, rcondref = 0.0;
  
  // WORK is REAL or DOUBLE PRECISION  array, dimension (4*N))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, 4*n, 0);
  
  // IWORK is INTEGER array, dimension (N)
  integer *iworkbuff = NULL, *iworkrefbuff = NULL;
  allocate_init_buffer(iworkbuff, iworkrefbuff, n, 0);
  
  // Call GETRF() to get A and IPIV buffers and pass these buffers to GECON()
  PRINTF("Calling GETRF() to get A, IPIV buffers and pass these buffers " \
         "to GECON()\n");
  
  // IPIV is INTEGER array, dimension (N)
  integer *ipivbuff, *ipivrefbuff;
  allocate_init_buffer(ipivbuff, ipivrefbuff, n, 0);
  
  /* M is INTEGER
          The number of rows of the matrix A.  M >= 0.*/
  // m is assigned with n as GECON() needs square matrix.
  integer m = lin_solver_paramslist[ip].n;
  if (m < 0) {
    PRINTF("m should be >= 0. Please correct the input data.\n");
  }
  
  // Call getrf_internal() to get buffers.
  integer info_cpp = -1;
  integer info_ref = -1;
  getrf_internal<T>(m, n, abuff, lda, ipivbuff,
              arefbuff, ipivrefbuff, &info_cpp, &info_ref);
  PRINTF ("getrf_internal() info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  if ((info_cpp < 0) && (info_ref < 0)) {
    PRINTF("Info returned by CPP or C API is not successful to get" \
           " A, IPIV buffers from GETRF(). Exiting....\n");
    closelibs();
    exit (-1);
  }
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("norm = %c\n", norm);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("anorm = %lf\n", anorm);
    PRINTF("Size of WORK array (4*n) = %d\n", 4*n);
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
    print_array<T>(arrayname, workbuff, 4 * n);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, 4 * n);
    
    // Prints IWORK array contents
    strncpy(arrayname, "IWORK input", arraysize);
    print_array<integer>(arrayname, iworkbuff, n);
    strncpy(arrayname, "IWORK ref input", arraysize);
    print_array<integer>(arrayname, iworkrefbuff, n);
  #endif
  
  info_cpp = -1;
  // Call CPP function
  libflame::gecon<T>(&norm, &n, abuff, &lda, &anorm, &rcond, workbuff,
                    iworkbuff, &info_cpp);
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(float)) {
    gecon_ref = (Fptr_NL_LAPACKE_gecon)dlsym(lapackModule, "sgecon_");
  } else if (typeid(T) == typeid(double)) {
    gecon_ref = (Fptr_NL_LAPACKE_gecon)dlsym(lapackModule, "dgecon_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (gecon_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  info_ref = -1;
  gecon_ref(&norm, &n, arefbuff, &lda, &anorm, &rcondref, workrefbuff,
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
      print_array<T>(arrayname, workbuff, 4 * n);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, 4 * n);
      
      // Prints IWORK array contents
      strncpy(arrayname, "IWORK output", arraysize);
      print_array<integer>(arrayname, iworkbuff, n);
      strncpy(arrayname, "IWORK ref output", arraysize);
      print_array<integer>(arrayname, iworkrefbuff, n);
    #endif
    double diff = computeError<T>(1, 1, &rcond, &rcondref);
    diff += computeError<T>(4, n, workbuff, workrefbuff);
    diff += computeError<integer>(1, n, iworkbuff, iworkrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
           " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] iworkbuff; delete[] iworkrefbuff;
  delete[] ipivbuff; delete[] ipivrefbuff;
}

/*! @brief  gecon_test_cmplx is function template for gecon() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  gecon_test_cmplx is function template for gecon() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  gecon_test_cmplx() function template calls C and CPP based lbrary APIs
	  with valid test values and returns the differences in output.
	  
    Complex reference:
	  http://www.netlib.org/lapack/explore-html/d4/d7e/group__complex_g_ecomputational_gaa2ad4e4b1c9cb56a23dd49a798aa9bc8.html#gaa2ad4e4b1c9cb56a23dd49a798aa9bc8
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d3/d01/group__complex16_g_ecomputational_gabe73145daeba3ec10e961054b75a07ce.html#gabe73145daeba3ec10e961054b75a07ce
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template<typename T, typename Ta>
void gecon_test_cmplx(int ip)
{
  typedef integer (*Fptr_NL_LAPACKE_gecon)(char* norm, integer* n, T* a,
                      integer* lda, Ta* anorm, Ta* rcond, T* work,
                      Ta* rwork, integer* info);
  Fptr_NL_LAPACKE_gecon gecon_ref;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
  /* NORM is CHARACTER*1
          Specifies whether the 1-norm condition number or the
          infinity-norm condition number is required:
          = '1' or 'O':  1-norm;
          = 'I':         Infinity-norm.*/
  char norm = lin_solver_paramslist[ip].norm;
  if ((norm != '1') && (norm != 'O') && (norm != 'I')) {
    PRINTF("norm should be 1 or O or I. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = lin_solver_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).*/
  integer lda = lin_solver_paramslist[ip].lda;
  if (lda < max(1, n)) {
    PRINTF("lda < max(1, n) but it should be: LDA >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // ANORM is REAL or DOUBLE PRECISION. The 1-norm of the original matrix A.
  Ta anorm =  (Ta)lin_solver_paramslist[ip].anorm;
  
  // RCOND is REAL or DOUBLE PRECISION.
  Ta rcond = 0.0, rcondref = 0.0;
  
  // WORK is COMPLEX or COMPLEX*16 array, dimension (2*N))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, 2*n, 0);
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (2*N))
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, 2*n, 0);
  
  // Call GETRF() to get A and IPIV buffers and pass these buffers to GECON()
  PRINTF("Calling GETRF() to get A, IPIV buffers and pass these buffers " \
         "to GECON()\n");
  
  // IPIV is INTEGER array, dimension (N)
  integer *ipivbuff, *ipivrefbuff;
  allocate_init_buffer(ipivbuff, ipivrefbuff, n, 0);
  
  /* M is INTEGER
          The number of rows of the matrix A.  M >= 0.*/
  integer m = lin_solver_paramslist[ip].m;
  if (m < 0) {
    PRINTF("m should be >= 0. Please correct the input data.\n");
  }
  
  // Call getrf_internal() to get buffers.
  integer info_cpp = -1;
  integer info_ref = -1;
  getrf_internal<T>(m, n, abuff, lda, ipivbuff,
              arefbuff, ipivrefbuff, &info_cpp, &info_ref);
  PRINTF ("getrf_internal() info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  if ((info_cpp < 0) && (info_ref < 0)) {
    PRINTF("Info returned by CPP or C API is not successful to get" \
           " A and IPIV buffers from GETRF(). Exiting....\n");
    closelibs();
    exit (-1);
  }
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("norm = %c\n", norm);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("anorm = %lf\n", anorm);
    PRINTF("Size of WORK array (2*n) = %d\n", 2*n);
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
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, 2 * n);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, 2 * n);
    
    // Prints RWORK array contents
    strncpy(arrayname, "RWORK input", arraysize);
    print_array<Ta>(arrayname, rworkbuff, 2 * n);
    strncpy(arrayname, "RWORK ref input", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, 2 * n);
  #endif
  
  info_cpp = -1;
  // Call CPP function
  libflame::gecon<T, Ta>(&norm, &n, abuff, &lda, &anorm, &rcond, workbuff,
                rworkbuff, &info_cpp);
  
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    gecon_ref = (Fptr_NL_LAPACKE_gecon)dlsym(lapackModule, "cgecon_");
  } else if (typeid(T) == typeid(dcomplex)) {
    gecon_ref = (Fptr_NL_LAPACKE_gecon)dlsym(lapackModule, "zgecon_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (gecon_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  info_ref = -1;
  gecon_ref(&norm, &n, arefbuff, &lda, &anorm, &rcondref, workrefbuff,
        rworkrefbuff, &info_ref);
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
      print_array<T>(arrayname, workbuff, 2 * n);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, 2 * n);
      
      // Prints RWORK array contents
      strncpy(arrayname, "RWORK output", arraysize);
      print_array<Ta>(arrayname, rworkbuff, 2 * n);
      strncpy(arrayname, "RWORK ref output", arraysize);
      print_array<Ta>(arrayname, rworkrefbuff, 2 * n);
    #endif
    double diff = computeError<Ta>(1, 1, &rcond, &rcondref);
    diff += computeError<T>(2, n, workbuff, workrefbuff);
    diff += computeError<Ta>(2, n, rworkbuff, rworkrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
           " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
  delete[] ipivbuff; delete[] ipivrefbuff;
}

/* Use TEST macro and call C++ test function template with
   float as typenames.*/
TEST(LAPACKCPP_gecon, SGECON) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    gecon_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double as typenames.*/
TEST(LAPACKCPP_gecon, DGECON) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    gecon_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_gecon, CGECON) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    gecon_test_cmplx<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_gecon, ZGECON) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    gecon_test_cmplx<dcomplex, double> (index);
  }
}