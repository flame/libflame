/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hbgst.cc
 *  @brief Test application to validate hbgst() using CPP template interface
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hbgst_test is function template for hbgst() functions.
            T can be scomplex, dcomplex
            Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hbgst_test is function template for hbevx() functions.
	  T can be scomplex, dcomplex
	  
	  hbgst_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO is >= 0.
    And passses the test case if difference is <= threshold.
    Fails the test case if difference > threshold or INFO < 0.
	  
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d3/db9/group__complex_o_t_h_e_rcomputational_ga808bf06bc4d353a18ab94f5eaf7c67f0.html#ga808bf06bc4d353a18ab94f5eaf7c67f0
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d0/da6/group__complex16_o_t_h_e_rcomputational_ga4c139408320128b94a42695614ae2646.html#ga4c139408320128b94a42695614ae2646
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template< typename T, typename Ta >
void hbgst_test(int ip)
{
  typedef integer (*Fptr_NL_LAPACKE_hbgst)(char* vect, char* uplo, integer* n,
                      integer* ka, integer* kb, T* ab, integer* ldab,  T* bb,
                      integer* ldbb, T* x, integer* ldx, T* work, Ta* rwork,
                      integer* info);
  Fptr_NL_LAPACKE_hbgst HBGST;

  // Initialise random number generators with timestamp.
  srand (time(NULL));
  
  /* N is INTEGER
          The order of the matrices A and B.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.\n");
  }
  
  /* VECT is CHARACTER*1
          = 'N':  do not form the transformation matrix X;
          = 'V':  form X.*/
  char vect = eig_paramslist[ip].vect_hbgst;
  if ((vect != 'N') && (vect != 'V')) {
    PRINTF("vect should be N or V. Please correct the input data.\n");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }

  /* KA is INTEGER
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.*/
  integer ka = 0;
  if (uplo == 'U') {
	  ka = eig_paramslist[ip].sda;
  } else if (uplo == 'L') {
	  ka = eig_paramslist[ip].subda;
  }
  if (ka < 0) {
    PRINTF("ka < 0 but it should be: KA >= 0. Please correct the input" \
           " data.\n");
  }
  
  /* KB is INTEGER
          The number of superdiagonals of the matrix B if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KA >= KB >= 0.*/
  integer kb = 0;
  if (uplo == 'U') {
	  kb = eig_paramslist[ip].sdb;
  } else if (uplo == 'L') {
	  kb = eig_paramslist[ip].subdb;
  }
  if (kb < 0) {
    PRINTF("kb < 0 but it should be: KB >= 0. Please correct the input" \
           " data.\n");
  }
  
  /* LDAB is INTEGER
          The leading dimension of the array AB.  LDAB >= KA+1.*/
  integer ldab = eig_paramslist[ip].ldab_hbgv;
  if (ldab < (ka + 1)) {
    PRINTF("ldab < (ka + 1) but should be: LDAB >= (KA + 1). Please correct" \
           " the input data.\n");
  }
  
  // AB is COMPLEX or COMPLEX*16 array, dimension (LDAB,N)
  T *abbuff, *abrefbuff;
  allocate_init_buffer(abbuff, abrefbuff, ldab * n);
  
  /* LDBB is INTEGER
          The leading dimension of the array BB.  LDBB >= KB+1.*/
  integer ldbb = eig_paramslist[ip].ldbb;
  if (ldbb < (kb + 1)) {
    PRINTF("ldbb < (kb + 1) but it should be: LDBB >= (KB + 1). Please " \
           "correct the input data.\n");
  }
  
  /* BB is COMPLEX or COMPLEX*16 array, dimension (LDBB,N) */
  T *bbbuff, *bbrefbuff;
  allocate_init_buffer(bbbuff, bbrefbuff, ldbb * n);
  
  /* LDX is INTEGER
          The leading dimension of the array X.
          LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise.*/
  integer ldx = eig_paramslist[ip].ldx;
  if (ldx < max(1, n)) {
    PRINTF("ldx < max(1, n) but it should be: LDX >= max(1,N). Please " \
           "correct the input data.\n");
  }
  if ((vect == 'V') && (ldx < 1)) {
    PRINTF("When Vect = V, ldx < 1 but it should be: LDX >=1. Please " \
           "correct the input data.\n");
  }
  
  /* X is COMPLEX or COMPLEX*16 array, dimension (LDX,N)
          If VECT = 'V', the n-by-n matrix X.
          If VECT = 'N', the array X is not referenced.*/
  T *xbuff, *xrefbuff; // output buffer
  allocate_init_buffer(xbuff, xrefbuff, ldx * n, 0);
  
  // WORK is COMPLEX array, dimension (N)
  T *workbuff, *workrefbuff;
  allocate_init_buffer(workbuff, workrefbuff, n, 0);
  
  // RWORK is REAL array, dimension (N)
  Ta *rworkbuff, *rworkrefbuff;
  allocate_init_buffer(rworkbuff, rworkrefbuff, n, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("vect = %c\n", vect);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("ka = %d\n", ka);
    PRINTF("kb = %d\n", kb);
    PRINTF("ldab = %d\n", ldab);
    PRINTF("Size of AB array (ldab*n) = %d\n", ldab*n);
    PRINTF("ldbb = %d\n", ldbb);
    PRINTF("Size of BB array (ldbb*n) = %d\n", ldbb*n);
    PRINTF("ldx = %d\n", ldx);
    PRINTF("Size of X array (ldx*n) = %d\n", ldx*n);
    PRINTF("Size of WORK array (n) = %d\n", n);
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
    
    // Prints AB array contents
    strncpy(arrayname, "AB input", arraysize);
    print_array<T>(arrayname, abbuff, ldab * n);
    strncpy(arrayname, "AB ref input", arraysize);
    print_array<T>(arrayname, abrefbuff, ldab * n);
    
    // Prints BB array contents
    strncpy(arrayname, "BB input", arraysize);
    print_array<T>(arrayname, bbbuff, ldbb * n);
    strncpy(arrayname, "BB ref input", arraysize);
    print_array<T>(arrayname, bbrefbuff, ldbb * n);
    
    // Prints X array contents
    strncpy(arrayname, "X input", arraysize);
    print_array<T>(arrayname, xbuff, ldx * n);
    strncpy(arrayname, "X ref input", arraysize);
    print_array<T>(arrayname, xrefbuff, ldx * n);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, n);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, n);
    
    // Prints RWORK array contents
    strncpy(arrayname, "RWORK input", arraysize);
    print_array<Ta>(arrayname, rworkbuff, n);
    strncpy(arrayname, "RWORK ref input", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, n);
  #endif
  
  // Call CPP function
  integer info_cpp = -1;
  libflame::hbgst<T, Ta>(&vect, &uplo, &n, &ka, &kb, abbuff, &ldab,
              bbbuff, &ldbb, xbuff, &ldx, workbuff, rworkbuff, &info_cpp);

  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    HBGST = (Fptr_NL_LAPACKE_hbgst)dlsym(lapackModule, "chbgst_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HBGST = (Fptr_NL_LAPACKE_hbgst)dlsym(lapackModule, "zhbgst_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HBGST == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }
  integer info_ref = -1;
  HBGST(&vect, &uplo, &n, &ka, &kb, abrefbuff, &ldab, bbrefbuff, &ldbb,
    xrefbuff, &ldx, workrefbuff, rworkrefbuff, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  if ((info_cpp >= 0) && (info_ref >= 0)) {
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all Output arrays contents...\n");
    
    // Prints AB array contents
    strncpy(arrayname, "AB output", arraysize);
    print_array<T>(arrayname, abbuff, ldab * n);
    strncpy(arrayname, "AB ref output", arraysize);
    print_array<T>(arrayname, abrefbuff, ldab * n);
    
    // Prints X array contents
    strncpy(arrayname, "X output", arraysize);
    print_array<T>(arrayname, xbuff, ldx * n);
    strncpy(arrayname, "X ref output", arraysize);
    print_array<T>(arrayname, xrefbuff, ldx * n);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK output", arraysize);
    print_array<T>(arrayname, workbuff, n);
    strncpy(arrayname, "WORK ref output", arraysize);
    print_array<T>(arrayname, workrefbuff, n);
    
    // Prints RWORK array contents
    strncpy(arrayname, "RWORK output", arraysize);
    print_array<Ta>(arrayname, rworkbuff, n);
    strncpy(arrayname, "RWORK ref output", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, n);
    #endif
    double diff = computeError<T>(ldab, n, abrefbuff, abbuff);
    if (vect == 'V') {
      diff += computeError<T>(ldx, n, xbuff, xrefbuff);
    }
    diff += computeError<T>(1, n, workbuff, workrefbuff);
    diff += computeError<Ta>(1, n, rworkbuff, rworkrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, abs(diff), SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abbuff; delete[] abrefbuff;
  delete[] bbbuff; delete[] bbrefbuff;
  delete[] xbuff; delete[] xrefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_hbgst, CHBGST) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hbgst_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_hbgst, ZHBGST) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hbgst_test<dcomplex, double> (index);
  }
}
