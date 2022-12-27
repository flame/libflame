/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_hegs2.cc
 *  @brief Test application to validate hegs2() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hegs2_test is function template for hegs2() functions.
			T can be scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  hegs2_test is function template for hegs2() functions.
	  T can be scomplex, dcomplex
	  
	  hegs2_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO = 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO <= 0.

 	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d4/d74/group__complex_h_ecomputational_ga853b4bda3e17fcec31189fa94aac6363.html#ga853b4bda3e17fcec31189fa94aac6363
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d3/d80/group__complex16_h_ecomputational_ga61ea380a543b4d76de88729a0ad42e60.html#ga61ea380a543b4d76de88729a0ad42e60
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void hegs2_test(int ip)
{
  typedef integer (*fptr_NL_LAPACK_hegs2)(integer* itype, char* uplo,
                      integer* n, T* a, integer* lda, T* b,
                      integer* ldb, integer* info);
  fptr_NL_LAPACK_hegs2 hegs2_ref = NULL;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* ITYPE is INTEGER
          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
          = 2 or 3: compute U*A*U**H or L**H *A*L.*/
  integer itype = lin_solver_paramslist[ip].itype;
  if ((itype < 1) || (itype > 3)) {
    PRINTF("itype should be from 1 to 3. Please correct the input data.\n");
  }
  
   /* UPLO is CHARACTER*1
          Specifies whether the upper or lower triangular part of the
          Hermitian matrix A is stored, and how B has been factorized.
          = 'U':  Upper triangular
          = 'L':  Lower triangular*/
  char uplo = lin_solver_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = lin_solver_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= fla_max(1,N).*/
  integer lda = lin_solver_paramslist[ip].lda;
  if (lda < fla_max(1, n)) {
    PRINTF("lda < fla_max(1, n) but it should be: LDA >= fla_max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  /* LDB is INTEGER
          The leading dimension of the array B.  LDB >= fla_max(1,N).*/
  integer ldb = lin_solver_paramslist[ip].ldb;
  if (ldb < fla_max(1, n)) {
    PRINTF("ldb < fla_max(1, n) but it should be: LDB >= fla_max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // B is COMPLEX or COMPLEX*16 array, dimension (LDB,N)
  T *bbuff, *brefbuff;
  allocate_init_buffer(bbuff, brefbuff, ldb * n);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("itype = %d\n", itype);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("ldb = %d\n", ldb);
    PRINTF("Size of B array (ldb*n) = %d\n", ldb * n);
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
    
    // Prints B array contents
    strncpy(arrayname, "B input", arraysize);
    print_array<T>(arrayname, bbuff, ldb * n);
    strncpy(arrayname, "B ref input", arraysize);
    print_array<T>(arrayname, brefbuff, ldb * n);
  #endif
  
  // Call CPP function
  integer info_cpp = -1;
  libflame::hegs2<T>(&itype, &uplo, &n, abuff, &lda, bbuff, &ldb, &info_cpp);
  
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    hegs2_ref = (fptr_NL_LAPACK_hegs2)dlsym(lapackModule, "chegs2_");
  } else if (typeid(T) == typeid(dcomplex)) {
    hegs2_ref = (fptr_NL_LAPACK_hegs2)dlsym(lapackModule, "zhegs2_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (hegs2_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  integer info_ref = -1;
  hegs2_ref(&itype, &uplo, &n, arefbuff, &lda, brefbuff, &ldb, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
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
      
      // Prints B array contents
      strncpy(arrayname, "B output", arraysize);
      print_array<T>(arrayname, bbuff, ldb * n);
      strncpy(arrayname, "B ref output", arraysize);
      print_array<T>(arrayname, brefbuff, ldb * n);
    #endif
    
    double diff = computeError<T>(lda, n, abuff, arefbuff);
    diff += computeError<T>(ldb, n, bbuff, brefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, LIN_SLVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  delete[] bbuff; delete[] brefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_hegs2, CHEGS2) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hegs2_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_hegs2, ZHEGS2) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hegs2_test<dcomplex> (index);
  }
}