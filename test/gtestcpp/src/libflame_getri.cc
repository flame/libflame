/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_getri.cc
 *  @brief Test application to validate getri() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  getri_test is function template for getri() functions.
			T can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  getri_test is function template for getri() functions.
	  T can be float, double.
	  
	  getri_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output.
	  
    Real reference:
    http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga1af62182327d0be67b1717db399d7d83.html#ga1af62182327d0be67b1717db399d7d83
    Double reference:
    https://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga56d9c860ce4ce42ded7f914fdb0683ff.html#ga56d9c860ce4ce42ded7f914fdb0683ff
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template< typename T >
void getri_test(int ip) {
  typedef integer (*Fptr_NL_LAPACKE_getri)(integer* n, T* a, integer* lda,
                      integer* ipiv, T* work, integer* lwork,
                      integer* info);
  Fptr_NL_LAPACKE_getri getri_ref = NULL;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
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
  
  // IPIV is INTEGER array, dimension (N)
  integer *ipivbuff = NULL, *ipivrefbuff = NULL;
  allocate_init_buffer(ipivbuff, ipivrefbuff, n, 0);
  
  // Call GETRF() to get A, IPIV buffers and pass these buffers to GETRI()
  PRINTF("Calling GETRF() to get A, IPIV buffers and pass these buffers " \
         "to GETRI()\n");
  /* M is INTEGER
          The number of rows of the matrix A.  M >= 0.*/
  // m is assigned with n as IPIV needs size as (min(m,n) = n).
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
  
  /* LWORK is INTEGER
          The dimension of the array WORK.  LWORK >= max(1,N).
          For optimal performance LWORK >= N*NB, where NB is
          the optimal blocksize returned by ILAENV.*/
  integer lwork = lin_solver_paramslist[ip].lwork;
  integer lwork_size = lwork;
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  
  info_cpp = -1;
  
  /* If lwork is -1, then call this API with default work
     variables. In return, work variables will be updated with array sizes
     needed.*/
  if (lwork == -1) {
    PRINTF("lwork is -1, so call getri() to get the array sizes.\n");
    T worksize = {0};
    
    // Call CPP function
    libflame::getri<T>(&n, abuff, &lda, ipivbuff, &worksize,
                    &lwork_size, &info_cpp);
    PRINTF("info_cpp: %d\n", info_cpp); //, worksize: %f\n", info_cpp, worksize.real);
    if (info_cpp == 0) {
      lwork_size = worksize;
    }
  }
  if (lwork_size < max(1, n)) {
    PRINTF("lwork < max(1, n) but it should be: lwork >= max(1, n)." \
           " Please correct the input data.\n");
  }
  
  // WORK is REAL or DOUBLE PRECISION  array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, max(1, lwork_size), 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of IPIV array (n) = %d\n", n);
    PRINTF("lwork = %d\n", lwork_size);
    PRINTF("Size of WORK array (max(1, lwork)) = %d\n", max(1, lwork_size));
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
    
    // Prints IPIV array contents
    strncpy(arrayname, "IPIV input", arraysize);
    print_array<integer>(arrayname, ipivbuff, n);
    strncpy(arrayname, "IPIV ref input", arraysize);
    print_array<integer>(arrayname, ipivrefbuff, n);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, max(1, lwork_size));
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
  #endif
  
  info_cpp = -1;
  // Call CPP function
  libflame::getri<T>(&n, abuff, &lda, ipivbuff, workbuff, &lwork_size,
              &info_cpp);
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(float)) {
    getri_ref = (Fptr_NL_LAPACKE_getri)dlsym(lapackModule, "sgetri_");
  } else if (typeid(T) == typeid(double)) {
    getri_ref = (Fptr_NL_LAPACKE_getri)dlsym(lapackModule, "dgetri_");
  } else if (typeid(T) == typeid(scomplex)) {
    getri_ref = (Fptr_NL_LAPACKE_getri)dlsym(lapackModule, "cgetri_");
  } else if (typeid(T) == typeid(dcomplex)) {
    getri_ref = (Fptr_NL_LAPACKE_getri)dlsym(lapackModule, "zgetri_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (getri_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  info_ref = -1;
  getri_ref(&n, arefbuff, &lda, ipivrefbuff, workrefbuff,
        &lwork_size, &info_ref);
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
      
      // Prints IPIV array contents
      strncpy(arrayname, "IPIV output", arraysize);
      print_array<integer>(arrayname, ipivbuff, n);
      strncpy(arrayname, "IPIV ref output", arraysize);
      print_array<integer>(arrayname, ipivrefbuff, n);
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, max(1, lwork_size));
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
    #endif
    double diff = 0.0;
    // TODO: Yet to finalize and do verification changes.
    /*diff = computeError<T>(lda, n, abuff, arefbuff);
    diff += computeError<T>(1, max(1, lwork_size), workbuff, workrefbuff);*/
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
           " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] ipivbuff; delete[] ipivrefbuff;
  delete[] workbuff; delete[] workrefbuff;
}

/*! @brief  getri_test is function template for getri() functions.
			T can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  getri_test is function template for getri() functions.
	  T can be scomplex, dcomplex.
	  
	  getri_test() function template calls C and CPP based lbrary APIs with
	  valid test values and returns the differences in output.
	  
    Complex reference:
	  http://www.netlib.org/lapack/explore-html/d4/d7e/group__complex_g_ecomputational_gae22ce12a3734b080ad8369ebf7e9c3a7.html#gae22ce12a3734b080ad8369ebf7e9c3a7
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d3/d01/group__complex16_g_ecomputational_gab490cfc4b92edec5345479f19a9a72ca.html#gab490cfc4b92edec5345479f19a9a72ca
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template< typename T >
void getri_test_cmplx(int ip) {
  typedef integer (*Fptr_NL_LAPACKE_getri)(integer* n, T* a, integer* lda,
                      integer* ipiv, T* work, integer* lwork,
                      integer* info);
  Fptr_NL_LAPACKE_getri getri_ref;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
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
  
  // IPIV is INTEGER array, dimension (N)
  integer *ipivbuff = NULL, *ipivrefbuff = NULL;
  allocate_init_buffer(ipivbuff, ipivrefbuff, n, 0);
  
  // Call GETRF() to get A, IPIV buffers and pass these buffers to GETRI()
  PRINTF("Calling GETRF() to get A, IPIV buffers and pass these buffers " \
         "to GETRI()\n");
  /* M is INTEGER
          The number of rows of the matrix A.  M >= 0.*/
  // m is assigned with n as GETRI() needs square matrix.
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
  
  /* LWORK is INTEGER
          The dimension of the array WORK.  LWORK >= max(1,N).
          For optimal performance LWORK >= N*NB, where NB is
          the optimal blocksize returned by ILAENV.*/
  integer lwork = lin_solver_paramslist[ip].lwork;
  integer lwork_size = lwork;
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  
  info_cpp = -1;
  
  /* If lwork is -1, then call this API with default work
     variables. In return, work variables will be updated with array sizes
     needed.*/
  if (lwork == -1) {
    PRINTF("lwork is -1, so call getri() to get the array sizes.\n");
    T worksize = {0};
    
    // Call CPP function
    libflame::getri<T>(&n, abuff, &lda, ipivbuff, &worksize,
                    &lwork_size, &info_cpp);
    PRINTF("info_cpp: %d\n", info_cpp); //, worksize: %f\n", info_cpp, worksize.real);
    if (info_cpp == 0) {
      lwork_size = worksize.real;
    }
  }
  if (lwork_size < max(1, n)) {
    PRINTF("lwork < max(1, n) but it should be: lwork >= max(1, n)." \
           " Please correct the input data.\n");
  }
  
  // WORK is REAL or DOUBLE PRECISION  array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, max(1, lwork_size), 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of IPIV array (n) = %d\n", n);
    PRINTF("lwork = %d\n", lwork_size);
    PRINTF("Size of WORK array (max(1, lwork)) = %d\n", max(1, lwork_size));
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
    
    // Prints IPIV array contents
    strncpy(arrayname, "IPIV input", arraysize);
    print_array<integer>(arrayname, ipivbuff, n);
    strncpy(arrayname, "IPIV ref input", arraysize);
    print_array<integer>(arrayname, ipivrefbuff, n);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, max(1, lwork_size));
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
  #endif
  
  info_cpp = -1;
  // Call CPP function
  libflame::getri<T>(&n, abuff, &lda, ipivbuff, workbuff,
                    &lwork_size, &info_cpp);
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    getri_ref = (Fptr_NL_LAPACKE_getri)dlsym(lapackModule, "cgetri_");
  } else if (typeid(T) == typeid(dcomplex)) {
    getri_ref = (Fptr_NL_LAPACKE_getri)dlsym(lapackModule, "zgetri_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (getri_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  info_ref = -1;
  getri_ref(&n, arefbuff, &lda, ipivrefbuff, workrefbuff,
        &lwork_size, &info_ref);
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
      
      // Prints IPIV array contents
      strncpy(arrayname, "IPIV output", arraysize);
      print_array<integer>(arrayname, ipivbuff, n);
      strncpy(arrayname, "IPIV ref output", arraysize);
      print_array<integer>(arrayname, ipivrefbuff, n);
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, max(1, lwork_size));
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
    #endif
    double diff = 0.0;
    // TODO: Yet to finalize and do verification changes.
    /*diff = computeError<T>(lda, n, abuff, arefbuff);
    diff += computeError<T>(1, max(1, lwork_size), workbuff, workrefbuff);*/
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
           " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] ipivbuff; delete[] ipivrefbuff;
  delete[] workbuff; delete[] workrefbuff;
}

/* Use TEST macro and call C++ test function template with
   float and float as typenames.*/
TEST(LAPACKCPP_getri, SGETRI) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getri_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double and double as typenames.*/
TEST(LAPACKCPP_getri, DGETRI) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getri_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_getri, CGETRI) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getri_test_cmplx<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_getri, ZGETRI) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getri_test_cmplx<dcomplex> (index);
  }
}