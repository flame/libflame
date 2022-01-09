/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hesv_rk.cc
 *  @brief Test application to validate hesv_rk() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hesv_rk_test is function template for hesv_rk() functions.
			T can be scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  hesv_rk_test is function template for hesv_rk() functions.
	  T can be scomplex, dcomplex
	  
	  hesv_rk_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d8/d8e/group__complex_h_esolve_ga4c2a4eceb9e7f2e5068c3ec3c14a9e88.html#ga4c2a4eceb9e7f2e5068c3ec3c14a9e88
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d3/d9d/group__complex16_h_esolve_ga15080de6926fb2099b184fb5c8367453.html#ga15080de6926fb2099b184fb5c8367453
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template< typename T >
void hesv_rk_test(int ip)
{
  typedef int (*fptr_NL_LAPACK_hesv_rk)(char* uplo, integer* n,
                  integer* nrhs, T* a, integer* lda, T* e, integer* ipiv,
                  T* b, integer* ldb, T* work, integer* lwork, integer* info);
  fptr_NL_LAPACK_hesv_rk hesv_rk_ref = NULL;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = lin_driver_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }

  /* N is INTEGER
          The order of the matrix A. N >= 0.*/
  integer n = lin_driver_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.");
  }

  /* NRHS is INTEGER. The number of right hand sides,
      i.e., the number of columns of the matrix B.  NRHS >= 0.*/
  integer nrhs = lin_driver_paramslist[ip].nrhs;
  if (nrhs < 0) {
    PRINTF("nrhs < 0 but should be: nrhs >= 0. Please correct the input data.");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).*/
  integer lda = lin_driver_paramslist[ip].lda;
  if (lda < max(1, n)) {
    PRINTF("lda < max(1, n) but it should be: LDA >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);

  // E is COMPLEX array, dimension (N)
  T *ebuff, *erefbuff;
  allocate_init_buffer(ebuff, erefbuff, n, 0);
  
  // IPIV is INTEGER array, dimension (N)
  integer *ipivbuff, *ipivrefbuff;
  allocate_init_buffer(ipivbuff, ipivrefbuff, n, 0);
  
  /* LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).*/
  integer ldb = lin_driver_paramslist[ip].ldb;
  if (ldb < max(1, n)) {
    PRINTF("ldb < max(1, n) but it should be: LDB >= max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // B is COMPLEX or COMPLEX16 array, dimension (LDB,NRHS)
  T *bbuff, *brefbuff;
  allocate_init_buffer(bbuff, brefbuff, ldb * nrhs);

  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = lin_driver_paramslist[ip].lwork;
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
    PRINTF("lwork is -1, so call hesv_rk() to get the array sizes.\n");
    T worksize = {0};
    
    // Call CPP function
    libflame::hesv_rk<T>(&uplo, &n, &nrhs, abuff, &lda, ebuff, ipivbuff,
                bbuff, &ldb, &worksize, &lwork_size, &info_cpp);
    PRINTF("info_cpp: %d, worksize: %f\n", info_cpp, worksize.real);
    if (info_cpp == 0) {
      if (lwork == -1) {
        lwork_size = worksize.real;
      }
    }
  }
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, max(1, lwork_size), 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("nrhs = %d\n", nrhs);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of E array (n) = %d\n", n);
    PRINTF("Size of IPIV array (n) = %d\n", n);
    PRINTF("ldb = %d\n", ldb);
    PRINTF("Size of B array (ldb*nrhs) = %d\n", ldb * nrhs);
    PRINTF("Size of WORK array (MAX(1, LWORK)) = %d\n", max(1, lwork_size));
    PRINTF("lwork = %d\n", lwork_size);
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
    
    // Prints E array contents
    strncpy(arrayname, "E input", arraysize);
    print_array<T>(arrayname, ebuff, n);
    strncpy(arrayname, "E ref input", arraysize);
    print_array<T>(arrayname, erefbuff, n);
    
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
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, max(1, lwork_size));
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
  #endif
  
  // Call CPP function
  info_cpp = -1;
  libflame::hesv_rk<T>(&uplo, &n, &nrhs, abuff, &lda, ebuff, ipivbuff,
                bbuff, &ldb, workbuff, &lwork_size, &info_cpp);
  
  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    hesv_rk_ref = (fptr_NL_LAPACK_hesv_rk)dlsym(lapackModule, "chesv_rk_");
  } else if (typeid(T) == typeid(dcomplex)) {
    hesv_rk_ref = (fptr_NL_LAPACK_hesv_rk)dlsym(lapackModule, "zhesv_rk_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (hesv_rk_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  
  integer info_ref = -1;
  hesv_rk_ref(&uplo, &n, &nrhs, arefbuff, &lda, erefbuff, ipivrefbuff, brefbuff,
      &ldb, workrefbuff, &lwork_size, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
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
      
      // Prints E array contents
      strncpy(arrayname, "E output", arraysize);
      print_array<T>(arrayname, ebuff, n);
      strncpy(arrayname, "E ref output", arraysize);
      print_array<T>(arrayname, erefbuff, n);
    
      // Prints IPIV array contents
      strncpy(arrayname, "IPIV output", arraysize);
      print_array<integer>(arrayname, ipivbuff, n);
      strncpy(arrayname, "IPIV ref output", arraysize);
      print_array<integer>(arrayname, ipivrefbuff, n);
      
      // Prints B array contents
      strncpy(arrayname, "B output", arraysize);
      print_array<T>(arrayname, bbuff, ldb * nrhs);
      strncpy(arrayname, "B ref output", arraysize);
      print_array<T>(arrayname, brefbuff, ldb * nrhs);
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, max(1, lwork_size));
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, max(1, lwork_size));
    #endif
    
    double diff = computeError<T>(lda, n, arefbuff, abuff);
    diff += computeError<integer>(1, n, ipivrefbuff, ipivbuff);
    diff += computeError<T>(1, n, erefbuff, ebuff);
    diff += computeError<T>(ldb, nrhs, brefbuff, bbuff);
    diff += computeError<T>(1, max(1, lwork_size), workrefbuff, workbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, abs(diff), LIN_DRVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers.
  delete[] abuff; delete[] arefbuff;
  delete[] bbuff; delete[] brefbuff;
  delete[] ebuff; delete[] erefbuff;
  delete[] ipivbuff; delete[] ipivrefbuff;
  delete[] workbuff; delete[] workrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_hesv_rk, CHESV_RK) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hesv_rk_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_hesv_rk, ZHESV_RK) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hesv_rk_test<dcomplex> (index);
  }
}