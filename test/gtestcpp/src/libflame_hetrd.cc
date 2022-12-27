/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_hetrd.cc
 *  @brief Test application to validate hetrd() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hetrd_test is function template for hetrd() functions.
			T can be scomplex, dcomplex
      Ta can be float, double
 * @details
 * \b Purpose:
    \verbatim
	  hetrd_test is function template for hetrd() functions.
	  T can be scomplex, dcomplex
    Ta can be float, double
	  
	  hetrd_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO = 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d4/d74/group__complex_h_ecomputational_gafacb74520e4816c134bc0b2ff61d25f1.html#gafacb74520e4816c134bc0b2ff61d25f1
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d3/d80/group__complex16_h_ecomputational_ga65f7a5eadb6a10738216bd47aafb49ad.html#ga65f7a5eadb6a10738216bd47aafb49ad
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template< typename T,  typename Ta >
void hetrd_test(int ip)
{
  typedef int (*Fptr_NL_LAPACK_hetrd)(char* uplo, integer* n, T* a,
                  integer* lda, Ta*  d, Ta*  e, T* tau, T* work,
                  integer* lwork, integer* info);
  Fptr_NL_LAPACK_hetrd HETRD;
  
  
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
  
  // D is REAL or DOUBLE PRECISION array, dimension (N)
  Ta *dbuff, *drefbuff;
  allocate_init_buffer(dbuff, drefbuff, n, 0);
  
  // E is REAL or DOUBLE PRECISION array, dimension (N-1)
  Ta *ebuff, *erefbuff;
  allocate_init_buffer(ebuff, erefbuff, n-1, 0);

  // TAU is COMPLEX or COMPLEX*16 array, dimension (N-1)
  T *taubuff, *taurefbuff;
  allocate_init_buffer(taubuff, taurefbuff, n-1, 0);
  
  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = eig_paramslist[ip].lwork_hetrd;
  integer lwork_size = lwork;
  
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  
  integer info_cpp = -1;
  
  /* If lwork is -1, then call this API with default work
     variables. In return, work variables will be updated with array sizes
     needed.*/
  if (lwork == -1) {
    PRINTF("lwork is -1, so call hetrd() to get the array sizes.\n");
    T worksize = {0};
    
    // Call CPP function
    libflame::hetrd<T, Ta>(&uplo, &n, abuff, &lda, dbuff, ebuff, taubuff,
                        &worksize, &lwork_size, &info_cpp);
    PRINTF("info_cpp: %d, worksize: %f\n", info_cpp, worksize.real);
    if (info_cpp == 0) {
      if (lwork == -1) {
        lwork_size = worksize.real;
      }
    }
  }
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, fla_max(1, lwork_size), 0);
  
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
    PRINTF("lwork = %d\n", lwork_size);
    PRINTF("Size of WORK array (MAX(1, LWORK)) = %d\n", fla_max(1, lwork_size));
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
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, fla_max(1, lwork_size));
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, fla_max(1, lwork_size));
  #endif
  
  // Call CPP function
  info_cpp = -1;
  libflame::hetrd<T, Ta>(&uplo, &n, abuff, &lda, dbuff, ebuff, taubuff,
                        workbuff, &lwork_size, &info_cpp);

  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    HETRD = (Fptr_NL_LAPACK_hetrd)dlsym(lapackModule, "chetrd_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HETRD = (Fptr_NL_LAPACK_hetrd)dlsym(lapackModule, "zhetrd_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HETRD == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  integer info_ref = -1;
  HETRD(&uplo, &n, arefbuff, &lda, drefbuff, erefbuff, taurefbuff, \
        workrefbuff, &lwork_size, &info_ref);
  
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
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, fla_max(1, lwork_size));
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, fla_max(1, lwork_size));
    #endif
  
    double diff = computeError<T>(lda, n, abuff, arefbuff);
    diff += computeError<Ta>(1, n, dbuff, drefbuff);
    diff += computeError<Ta>(1, n-1, ebuff, erefbuff);
    diff += computeError<T>(1, n-1, taubuff, taurefbuff);
    diff += computeError<T>(1, fla_max(1, lwork_size), workrefbuff, workbuff);
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
  delete[] workbuff; delete[] workrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex, float as typenames.*/
TEST(LAPACKCPP_hetrd, CHETRD) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetrd_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex, double as typenames.*/
TEST(LAPACKCPP_hetrd, ZHETRD) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetrd_test<dcomplex, double> (index);
  }
}