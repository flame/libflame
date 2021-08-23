/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_heequb.cc
 *  @brief Test application to validate heequb() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  heequb_test is function template for heequb() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  heequb_test is function template for heequb() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  heequb_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO is >= 0.
    And passses the test case if difference is <= threshold.
    Fails the test case if difference > threshold or INFO < 0.
	  
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d4/d74/group__complex_h_ecomputational_ga986174490b3d9eb0d10502d96883e153.html#ga986174490b3d9eb0d10502d96883e153
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d3/d80/group__complex16_h_ecomputational_ga0626d54efa3610a40077cf6685df73f1.html#ga0626d54efa3610a40077cf6685df73f1
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
		      Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template< typename T, typename Ta >
void heequb_test(int ip)
{
  typedef int (*fptr_NL_LAPACK_heequb)(char* uplo, integer* n, T* a,
                  integer* lda, Ta* s, Ta* scond, Ta* amax, T* work,
                  integer* info);
  fptr_NL_LAPACK_heequb heequb_ref;
  
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
          The order of the matrix A. N >= 0.*/
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

  // S is REAL or DOUBLE PRECISION array, dimension (N).
  Ta *sbuff, *srefbuff;
  allocate_init_buffer(sbuff, srefbuff, n);

  // scond and amax are REAL or DOUBLE PRECISION.
  Ta scond = 0.0, scondref = 0.0;
  Ta amax = 0.0, amaxref = 0.0;
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (2*N))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, 2*n);
  
  // Print input values other than arrays.
  #if PRINT_INPUT_VALUES
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda*n);
    PRINTF("Size of S array (n) = %d\n", n);
    PRINTF("scond = %f\n", scond);
    PRINTF("amax = %f\n", amax);
    PRINTF("Size of WORK array (2*n) = %d\n", 2*n);
  #endif

  #if PRINT_ARRAYS
  // Array to store array name to print.
  char arrayname[20] = "";
  integer arraysize = sizeof(arrayname);
  #endif
  
  #if (PRINT_ARRAYS && PRINT_INPUT_ARRAYS)
    // Print all input arrays if PRINT_INPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all Input arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A input", arraysize);
    print_array<T>(arrayname, abuff, n);
    strncpy(arrayname, "A ref input", arraysize);
    print_array<T>(arrayname, arefbuff, n);
    
    // Prints S array contents
    strncpy(arrayname, "S input", arraysize);
    print_array<Ta>(arrayname, sbuff, n);
    strncpy(arrayname, "S ref input", arraysize);
    print_array<Ta>(arrayname, srefbuff, n);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, 2 * n);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, 2 * n);
  #endif
  
  integer info_cpp = -1;
  
  // Call CPP function
  libflame::heequb<T, Ta>(&uplo, &n, abuff, &lda, sbuff, &scond, &amax,
              workbuff, &info_cpp);

  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    heequb_ref = (fptr_NL_LAPACK_heequb)dlsym(lapackModule, "cheequb_");
  } else if (typeid(T) == typeid(dcomplex)) {
    heequb_ref = (fptr_NL_LAPACK_heequb)dlsym(lapackModule, "zheequb_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (heequb_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }
  
  integer info_ref = -1;
  // Call C function
  heequb_ref(&uplo, &n, arefbuff, &lda, srefbuff, &scondref, &amaxref,
        workrefbuff, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  if ((info_cpp >= 0) && (info_ref >= 0)) {
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
      // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
      PRINTF("\nPrinting all Output arrays contents...\n");
      
      // Prints S array contents
      strncpy(arrayname, "S output", arraysize);
      print_array<Ta>(arrayname, sbuff, n);
      strncpy(arrayname, "S ref output", arraysize);
      print_array<Ta>(arrayname, srefbuff, n);
      
      PRINTF ("scond: %lf, scondref: %lf\n", scond, scondref);
      PRINTF ("amax: %lf, amaxref: %lf\n", amax, amaxref);
  
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, 2 * n);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, 2 * n);
    #endif
    
    double diff = computeError<Ta>(1, 1, &amaxref, &amax);
    diff = computeError<Ta>(1, 1, &scondref, &scond);
	  diff += computeError<Ta>(1, n, srefbuff, sbuff);
    diff += computeError<T>(2, n, workbuff, workrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] sbuff; delete[] srefbuff;
  delete[] workbuff; delete[] workrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_heequb, CHEEQUB) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    heequb_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_heequb, ZHEEQUB) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    heequb_test<dcomplex, double> (index);
  }
}