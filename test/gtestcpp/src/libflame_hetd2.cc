/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_hetd2.cc
 *  @brief Test application to validate hetd2() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hetd2_test is function template for hetd2() functions.
			T can be scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  hetd2_test is function template for hetd2() functions.
	  T can be scomplex, dcomplex
    Ta can be float, double
	  
	  hetd2_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d4/d74/group__complex_h_ecomputational_ga4a93e0522d4d3aa68c54c9b6ebdfbce9.html#ga4a93e0522d4d3aa68c54c9b6ebdfbce9
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d3/d80/group__complex16_h_ecomputational_ga89e75b4e5009c2b69a142563eeb9ea33.html#ga89e75b4e5009c2b69a142563eeb9ea33
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template< typename T,  typename Ta >
void hetd2_test(int ip)
{
  typedef int (*Fptr_NL_LAPACK_hetd2)(char* uplo, integer* n, T* a,
                  integer* lda, Ta*  d, Ta*  e, T* tau, integer* info);
  Fptr_NL_LAPACK_hetd2 HETD2;
  
  
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
          The order of the matrices A and B.  N >= 0.*/
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
  
  // D is REAL or DOUBLE PRECISION array, dimension (N)
  Ta *dbuff, *drefbuff;
  allocate_init_buffer(dbuff, drefbuff, n, 0);
  
  // E is REAL or DOUBLE PRECISION array, dimension (N-1)
  Ta *ebuff, *erefbuff;
  allocate_init_buffer(ebuff, erefbuff, n-1, 0);

  // TAU is COMPLEX or COMPLEX*16 array, dimension (N-1)
  T *taubuff, *taurefbuff;
  allocate_init_buffer(taubuff, taurefbuff, n-1, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of D array (n) = %d\n", n);
    PRINTF("Size of E array (n-1) = %d\n", n-1);
    PRINTF("Size of TAU array (n-1) = %d\n", n-1);
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
    
    // Prints D array contents
    strncpy(arrayname, "D input", arraysize);
    print_array<Ta>(arrayname, dbuff, n);
    strncpy(arrayname, "D ref input", arraysize);
    print_array<Ta>(arrayname, drefbuff, n);
    
    // Prints E array contents
    strncpy(arrayname, "E input", arraysize);
    print_array<Ta>(arrayname, ebuff, n-1);
    strncpy(arrayname, "E ref input", arraysize);
    print_array<Ta>(arrayname, erefbuff, n-1);
    
    // Prints TAU array contents
    strncpy(arrayname, "TAU input", arraysize);
    print_array<T>(arrayname, taubuff, n-1);
    strncpy(arrayname, "TAU ref input", arraysize);
    print_array<T>(arrayname, taurefbuff, n-1);
  #endif
  
  // Call CPP function
  integer info_cpp = -1;
  libflame::hetd2<T, Ta>(&uplo, &n, abuff, &lda, dbuff, ebuff, taubuff,
                        &info_cpp);

  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    HETD2 = (Fptr_NL_LAPACK_hetd2)dlsym(lapackModule, "chetd2_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HETD2 = (Fptr_NL_LAPACK_hetd2)dlsym(lapackModule, "zhetd2_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HETD2 == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  integer info_ref = -1;
  HETD2(&uplo, &n, arefbuff, &lda, drefbuff, erefbuff, taurefbuff, &info_ref);
  
  // Calculate the differences of buffers.
  if ((info_cpp == 0) && (info_ref == 0)) {
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
      // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
      PRINTF("\nPrinting all Output arrays contents...\n");
      // Prints A array contents
      strncpy(arrayname, "A output", arraysize);
      print_array<T>(arrayname, abuff, lda * n);
      strncpy(arrayname, "A ref output", arraysize);
      print_array<T>(arrayname, arefbuff, lda * n);
      
      // Prints D array contents
      strncpy(arrayname, "D output", arraysize);
      print_array<Ta>(arrayname, dbuff, n);
      strncpy(arrayname, "D ref output", arraysize);
      print_array<Ta>(arrayname, drefbuff, n);
      
      // Prints E array contents
      strncpy(arrayname, "E output", arraysize);
      print_array<Ta>(arrayname, ebuff, n-1);
      strncpy(arrayname, "E ref output", arraysize);
      print_array<Ta>(arrayname, erefbuff, n-1);
      
      // Prints TAU array contents
      strncpy(arrayname, "TAU output", arraysize);
      print_array<T>(arrayname, taubuff, n-1);
      strncpy(arrayname, "TAU ref output", arraysize);
      print_array<T>(arrayname, taurefbuff, n-1);
    #endif
  
    double diff = computeError<T>(lda, n, abuff, arefbuff);
    diff += computeError<Ta>(1, n, dbuff, drefbuff);
    diff += computeError<Ta>(1, n-1, ebuff, erefbuff);
    diff += computeError<T>(1, n-1, taubuff, taurefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] dbuff; delete[] drefbuff;
  delete[] ebuff; delete[] erefbuff;
  delete[] taubuff; delete[] taurefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex, float as typenames.*/
TEST(LAPACKCPP_hetd2, CHETD2) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetd2_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex, double as typenames.*/
TEST(LAPACKCPP_hetd2, ZHETD2) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetd2_test<dcomplex, double> (index);
  }
}