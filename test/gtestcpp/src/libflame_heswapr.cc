/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_heswapr.cc
 *  @brief Test application to validate heswapr() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  heswapr_test is function template for heswapr() functions.
			T can be scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  heswapr_test is function template for heswapr() functions.
	  T can be scomplex, dcomplex
	  
	  heswapr_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/dd/d78/group__complex_h_eauxiliary_ga955ef3394562af1ae4ae2d113e4423bd.html#ga955ef3394562af1ae4ae2d113e4423bd
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d8/dab/group__complex16_h_eauxiliary_ga37927cbfd870be1c27dd7438c9b1e61f.html#ga37927cbfd870be1c27dd7438c9b1e61f
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template< typename T >
void heswapr_test(int ip)
{
  typedef int (*Fptr_NL_LAPACK_heswapr)(char* uplo, integer* n,
                  T* a, integer* lda, integer* i1, integer* i2);
  Fptr_NL_LAPACK_heswapr HESWAPR;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* UPLO is CHARACTER*1
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix.
          = 'U':  Upper triangular, form is A = U*D*U**T;
          = 'L':  Lower triangular, form is A = L*D*L**T.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }

  /* N is INTEGER
          The order of the matrices A and B.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.");
  }
  
  // LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
  integer lda = eig_paramslist[ip].lda;
  if (lda < max(1, n)) {
    PRINTF("lda < max(1, n) but it should be: LDA >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA, N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  // I1 is INTEGER. Index of the first row to swap
  integer i1 = eig_paramslist[ip].i1;
  
  // I2 is INTEGER. Index of the second row to swap
  integer i2 = eig_paramslist[ip].i2;
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("i1 = %d\n", i1);
    PRINTF("i2 = %d\n", i2);
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
  #endif
  
  // Call CPP function
  libflame::heswapr<T>(&uplo, &n, abuff, &lda, &i1, &i2);

  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    HESWAPR = (Fptr_NL_LAPACK_heswapr)dlsym(lapackModule, \
                        "cheswapr_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HESWAPR = (Fptr_NL_LAPACK_heswapr)dlsym(lapackModule, \
                        "zheswapr_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HESWAPR == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  
  HESWAPR(&uplo, &n, arefbuff, &lda, &i1, &i2);
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all Output arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A output", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref output", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
  #endif
  
  double diff = computeError<T>(lda, n, arefbuff, abuff);
  EXPECT_NEAR(0.0, abs(diff), SYM_EIGEN_THRESHOLD);
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_heswapr, CHESWAPR) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    heswapr_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_heswapr, ZHESWAPR) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    heswapr_test<dcomplex> (index);
  }
}