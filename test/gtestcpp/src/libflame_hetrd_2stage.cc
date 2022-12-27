/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_hetrd_2stage.cc
 *  @brief Test application to validate hetrd_2stage() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hetrd_2stage_test is function template for hetrd_2stage() functions.
			T can be scomplex, dcomplex
      Ta can be float, double
 * @details
 * \b Purpose:
    \verbatim
	  hetrd_2stage_test is function template for hetrd_2stage() functions.
	  T can be scomplex, dcomplex
    Ta can be float, double
	  
	  hetrd_2stage_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO = 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d4/d74/group__complex_h_ecomputational_gaf3e33440fb683b215f6c2569869d6965.html#gaf3e33440fb683b215f6c2569869d6965
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d3/d80/group__complex16_h_ecomputational_ga15264199d62f32abbd25a5b880b62209.html#ga15264199d62f32abbd25a5b880b62209
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template< typename T,  typename Ta >
void hetrd_2stage_test(int ip)
{
  typedef int (*Fptr_NL_LAPACK_hetrd_2stage)(char *vect, char* uplo,
                  integer* n, T* a, integer* lda, Ta*  d, Ta*  e, T* tau,
                  T *hous2, integer *lhous2, T* work, integer* lwork,
                  integer* info);
  Fptr_NL_LAPACK_hetrd_2stage HETRD_2STAGE;
  
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* VECT is CHARACTER*1
          = 'N':  No need for the Housholder representation, 
                  in particular for the second stage (Band to
                  tridiagonal) and thus LHOUS2 is of size fla_max(1, 4*N);
          = 'V':  the Householder representation is needed to 
                  either generate Q1 Q2 or to apply Q1 Q2, 
                  then LHOUS2 is to be queried and computed.
                  (NOT AVAILABLE IN THIS RELEASE).*/
  char vect = eig_paramslist[ip].vect_hetrd_2stage;
  if ((vect != 'N') && (vect != 'V')) {
    PRINTF("vect should be N or V. Please correct the input data.\n");
  }
  
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
  
  /* LHOUS2 is INTEGER
          The dimension of the array HOUS2.
          If LWORK = -1, or LHOUS2=-1,
          then a query is assumed; the routine
          only calculates the optimal size of the HOUS2 array, returns
          this value as the first entry of the HOUS2 array, and no error
          message related to LHOUS2 is issued by XERBLA.
          If VECT='N', LHOUS2 = fla_max(1, 4*n);
          if VECT='V', option not yet available.*/
  integer lhous2 = eig_paramslist[ip].lhous2;
  integer lhous2_size = lhous2;
  if ((lhous2 < -1) || (lhous2 == 0)) {
    PRINTF("lhous2 is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  
  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = eig_paramslist[ip].lwork_hetrd_2stage;
  integer lwork_size = lwork;
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  
  integer info_cpp = -1;
  
  /* If lwork/lhous2 is -1, then call this API with default work
     variables. In return, work variables will be updated with array sizes
     needed.*/
  if ((lwork == -1) || (lhous2 == -1)) {
    PRINTF("lwork/lhous2 is -1, so call hetrd_2stage() to get the array" \
          " sizes.\n");
    T lhous2size = {0};
    T worksize = {0};
    
    // Call CPP function
    libflame::hetrd_2stage<T, Ta>(&vect, &uplo, &n, abuff, &lda, dbuff,
                        ebuff, taubuff, &lhous2size, &lhous2_size, &worksize,
                        &lwork_size, &info_cpp);
    PRINTF("info_cpp: %d, lhous2size: %f, worksize: %f\n", info_cpp, \
          lhous2size.real, worksize.real);
    if (info_cpp == 0) {
      if (lwork == -1) {
        lwork_size = worksize.real;
      }
      if (lhous2 == -1) {
        lhous2_size = lhous2size.real;
      }
    }
  }
  
  // HOUS2 is COMPLEX or COMPLEX*16  array, dimension (LHOUS2))
  T *hous2buff = NULL, *hous2refbuff = NULL;
  allocate_init_buffer(hous2buff, hous2refbuff, lhous2_size, 0);
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (LWORK)
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, lwork_size, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("vect = %c\n", vect);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of D array (n) = %d\n", n);
    PRINTF("Size of E array (n-1) = %d\n", n-1);
    PRINTF("Size of TAU array (n-1) = %d\n", n-1);
    PRINTF("lhous2 = %d\n", lhous2_size);
    PRINTF("Size of HOUS2 array (lhous2) = %d\n", lhous2_size);
    PRINTF("lwork = %d\n", lwork_size);
    PRINTF("Size of WORK array (LWORK) = %d\n", lwork_size);
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
    
    // Prints HOUS2 array contents
    strncpy(arrayname, "HOUS2 input", arraysize);
    print_array<T>(arrayname, hous2buff, lhous2_size);
    strncpy(arrayname, "HOUS2 ref input", arraysize);
    print_array<T>(arrayname, hous2refbuff, lhous2_size);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, lwork_size);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, lwork_size);
  #endif
  
  // Call CPP function
  info_cpp = -1;
  libflame::hetrd_2stage<T, Ta>(&vect, &uplo, &n, abuff, &lda, dbuff, ebuff,
                        taubuff, hous2buff, &lhous2_size, workbuff,
                        &lwork_size, &info_cpp);

  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    HETRD_2STAGE = (Fptr_NL_LAPACK_hetrd_2stage)dlsym(lapackModule, \
                          "chetrd_2stage_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HETRD_2STAGE = (Fptr_NL_LAPACK_hetrd_2stage)dlsym(lapackModule, \
                          "zhetrd_2stage_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HETRD_2STAGE == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  integer info_ref = -1;
  HETRD_2STAGE(&vect, &uplo, &n, arefbuff, &lda, drefbuff, erefbuff,
        taurefbuff, hous2refbuff, &lhous2_size, workrefbuff, &lwork_size,
        &info_ref);
  
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
      
      // Prints HOUS2 array contents
      strncpy(arrayname, "HOUS2 output", arraysize);
      print_array<T>(arrayname, hous2buff, lhous2_size);
      strncpy(arrayname, "HOUS2 ref output", arraysize);
      print_array<T>(arrayname, hous2refbuff, lhous2_size);
    
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, lwork_size);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, lwork_size);
    #endif
  
    double diff = computeError<T>(lda, n, abuff, arefbuff);
    diff += computeError<Ta>(1, n, dbuff, drefbuff);
    diff += computeError<Ta>(1, n-1, ebuff, erefbuff);
    diff += computeError<T>(1, n-1, taubuff, taurefbuff);
    diff += computeError<T>(1, lhous2_size, hous2buff, hous2refbuff);
    diff += computeError<T>(1, lwork_size, workrefbuff, workbuff);
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
  delete[] hous2buff; delete[] hous2refbuff;
  delete[] workbuff; delete[] workrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex, float as typenames.*/
TEST(LAPACKCPP_hetrd_2stage, CHETRD_2STAGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetrd_2stage_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex, double as typenames.*/
TEST(LAPACKCPP_hetrd_2stage, ZHETRD_2STAGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetrd_2stage_test<dcomplex, double> (index);
  }
}