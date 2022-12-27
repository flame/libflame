/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_hetf2_rook.cc
 *  @brief Test application to validate hetf2_rook() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hetf2_rook_test is function template for hetf2_rook() functions.
			T can be scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  hetf2_rook_test is function template for hetf2_rook() functions.
	  T can be scomplex, dcomplex
	  
	  hetf2_rook_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d4/d74/group__complex_h_ecomputational_ga068309e57c51f1fa0171ca3d93b5848f.html#ga068309e57c51f1fa0171ca3d93b5848f
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d3/d80/group__complex16_h_ecomputational_ga872b3466c66763af705ef67440976ce3.html#ga872b3466c66763af705ef67440976ce3
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void hetf2_rook_test(int ip)
{
  typedef int (*Fptr_NL_LAPACK_hetf2_rook)(char *uplo, integer *n, T *a,
                  integer *lda, integer *ipiv, integer *info);
  Fptr_NL_LAPACK_hetf2_rook HETF2_ROOK;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* UPLO is CHARACTER*1
          Specifies whether the upper or lower triangular part of the
          Hermitian matrix A is stored:
          = 'U':  Upper triangular
          = 'L':  Lower triangular*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }

  /* N is INTEGER
          The order of the matrix A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.");
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
  
  // IPIV is INTEGER array, dimension (N)
  integer *ipivbuff, *ipivrefbuff;
  allocate_init_buffer(ipivbuff, ipivrefbuff, n, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of IPIV array (n) = %d\n", n);
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
  #endif
  
  // Call CPP function
  integer info_cpp = -1;
  libflame::hetf2_rook<T>(&uplo, &n, abuff, &lda, ipivbuff, &info_cpp);

  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    HETF2_ROOK = (Fptr_NL_LAPACK_hetf2_rook)dlsym(lapackModule, "chetf2_rook_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HETF2_ROOK = (Fptr_NL_LAPACK_hetf2_rook)dlsym(lapackModule, "zhetf2_rook_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HETF2_ROOK == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  integer info_ref = -1;
  HETF2_ROOK(&uplo, &n, arefbuff, &lda, ipivrefbuff, &info_ref);
  
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
    #endif
  
    double diff = computeError<T>(lda, n, abuff, arefbuff);
    diff += computeError<integer>(1, n, ipivbuff, ipivrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] ipivbuff; delete[] ipivrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex, float as typenames.*/
TEST(LAPACKCPP_hetf2_rook, CHETF2_ROOK) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetf2_rook_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex, double as typenames.*/
TEST(LAPACKCPP_hetf2_rook, ZHETF2_ROOK) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetf2_rook_test<dcomplex> (index);
  }
}