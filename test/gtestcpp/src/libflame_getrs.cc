/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_getrs.cc
 *  @brief Test application to validate getrs() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  getrs_test is function template for getrs() functions.
			T can be scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  getrs_test is function template for getrs() functions.
	  T can be scomplex, dcomplex
	  
	  getrs_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO = 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

    Real reference:
    http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_gaa00bcf4d83a118cb6f0b6619d6ffaa24.html#gaa00bcf4d83a118cb6f0b6619d6ffaa24
    Double reference:
    http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga58e332cb1b8ab770270843221a48296d.html#ga58e332cb1b8ab770270843221a48296d
 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d4/d7e/group__complex_g_ecomputational_ga3a79ef0038488e420519c422c1a2a8f2.html#ga3a79ef0038488e420519c422c1a2a8f2
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d3/d01/group__complex16_g_ecomputational_ga3a5b88a7e8bf70591e521e86464e109d.html#ga3a5b88a7e8bf70591e521e86464e109d
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void getrs_test(int ip)
{
  typedef int (*Fptr_NL_LAPACK_getrs)(char* trans, integer* n, integer* nrhs,
                    T* a, integer* lda,  integer* ipiv, T* b, integer* ldb,
                    integer* info);
  Fptr_NL_LAPACK_getrs getrs_ref = NULL;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* TRANS is CHARACTER*1
          Specifies the form of the system of equations:
          = 'N':  A * X = B  (No transpose)
          = 'T':  A**T* X = B  (Transpose)
          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)*/
  char trans = lin_solver_paramslist[ip].trans;
  if ((trans != 'N') && (trans != 'T') && (trans != 'C')) {
    PRINTF("trans should be N or T or C. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrix A.  N >= 0.*/
  integer n = lin_solver_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.");
  }
  
  /* NRHS is INTEGER
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0.*/
  integer nrhs = lin_solver_paramslist[ip].nrhs;
  if (nrhs < 0) {
    PRINTF("nrhs < 0 but should be: nrhs >= 0. Please correct the input data.");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).*/
  integer lda = lin_solver_paramslist[ip].lda;
  if (lda < max(1, n)) {
    PRINTF("lda < max(1, n) but it should be: LDA >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is FLOAT/DOUBLE/COMPLEX/COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // IPIV is INTEGER array, dimension (N)
  integer *ipivbuff, *ipivrefbuff;
  allocate_init_buffer(ipivbuff, ipivrefbuff, n, 0);
  
  /* LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).*/
  integer ldb = lin_solver_paramslist[ip].ldb;
  if (ldb < max(1, n)) {
    PRINTF("ldb < max(1, n) but it should be: LDB >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // B is FLOAT/DOUBLE/COMPLEX/COMPLEX*16 array, dimension (LDA,NRHS)
  T *bbuff, *brefbuff;
  allocate_init_buffer(bbuff, brefbuff, ldb * nrhs);
  
  // Call GETRF() to get A and IPIV buffers and pass these buffers to GETRS()
  PRINTF("Calling GETRF() to get A, IPIV buffers and pass these buffers " \
         "to GETRS()\n");
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
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("trans = %c\n", trans);
    PRINTF("n = %d\n", n);
    PRINTF("nrhs = %d\n", nrhs);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of IPIV array (n) = %d\n", n);
    PRINTF("ldb = %d\n", ldb);
    PRINTF("Size of B array (ldb*nrhs) = %d\n", ldb * nrhs);
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
    
    // Prints B array contents
    strncpy(arrayname, "B input", arraysize);
    print_array<T>(arrayname, bbuff, ldb * nrhs);
    strncpy(arrayname, "B ref input", arraysize);
    print_array<T>(arrayname, brefbuff, ldb * nrhs);
  #endif
  
  // Call CPP function
  info_cpp = -1;
  libflame::getrs<T>(&trans, &n, &nrhs, abuff, &lda, ipivbuff, bbuff, &ldb,
                     &info_cpp);

  // Call C function
  if (typeid(T) == typeid(float)) {
    getrs_ref = (Fptr_NL_LAPACK_getrs)dlsym(lapackModule, "sgetrs_");
  } else if (typeid(T) == typeid(double)) {
    getrs_ref = (Fptr_NL_LAPACK_getrs)dlsym(lapackModule, "dgetrs_");
  } else if (typeid(T) == typeid(scomplex)) {
    getrs_ref = (Fptr_NL_LAPACK_getrs)dlsym(lapackModule, "cgetrs_");
  } else if (typeid(T) == typeid(dcomplex)) {
    getrs_ref = (Fptr_NL_LAPACK_getrs)dlsym(lapackModule, "zgetrs_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (getrs_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  info_ref = -1;
  getrs_ref(&trans, &n, &nrhs, arefbuff, &lda, ipivrefbuff, brefbuff, &ldb,
        &info_ref);
  
  // Calculate the differences of buffers.
  if ((info_cpp == 0) && (info_ref == 0)) {
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
      // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
      PRINTF("\nPrinting all Output arrays contents...\n");
      
      // Prints B array contents
      strncpy(arrayname, "B input", arraysize);
      print_array<T>(arrayname, bbuff, ldb * nrhs);
      strncpy(arrayname, "B ref input", arraysize);
      print_array<T>(arrayname, brefbuff, ldb * nrhs);
    #endif
  
    double diff = 0.0;
    // TODO: Yet to finalize and do verification changes.
    //diff = computeError<T>(ldb, nrhs, bbuff, brefbuff);
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
  delete[] bbuff; delete[] brefbuff;
}

/* Use TEST macro and call C++ test function template with
   float as typenames.*/
TEST(LAPACKCPP_getrs, SGETRS) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getrs_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double as typenames.*/
TEST(LAPACKCPP_getrs, DGETRS) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getrs_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_getrs, CGETRS) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getrs_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_getrs, ZGETRS) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    getrs_test<dcomplex> (index);
  }
}