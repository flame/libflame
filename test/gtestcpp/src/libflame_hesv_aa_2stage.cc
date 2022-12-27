/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hesv_aa_2stage.cc
 *  @brief Test application to validate hesv_aa_2stage() using CPP template
           interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hesv_aa_2stage_test is function template for hesv_aa_2stage() functions.
			T can be scomplex, dcomplex
 * @details
 * \b Purpose:
    \verbatim
	  hesv_aa_2stage_test is function template for hesv_aa_2stage() functions.
	  T can be scomplex, dcomplex
	  
	  hesv_aa_2stage_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/da/dff/group__complex_s_ycomputational_gabb66fb23be3a7311b71271a2717b35eb.html#gabb66fb23be3a7311b71271a2717b35eb
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d3/d9d/group__complex16_h_esolve_gada4828eb3ecee73a77548a48357e0879.html#gada4828eb3ecee73a77548a48357e0879
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array present in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void hesv_aa_2stage_test(int ip)
{
  typedef int (*fptr_NL_LAPACK_hesv_aa_2stage)(char* uplo, integer* n,
                  integer* nrhs, T* a, integer* lda, T* tb, integer* ltb,
                  integer* ipiv, integer* ipiv2, T* b, integer* ldb,
                  T* work, integer* lwork, integer* info);
  fptr_NL_LAPACK_hesv_aa_2stage hesv_aa_2stage_ref = NULL;
  
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
          The leading dimension of the array A.  LDA >= fla_max(1,N).*/
  integer lda = lin_driver_paramslist[ip].lda;
  if (lda < fla_max(1, n)) {
    PRINTF("lda < fla_max(1, n) but it should be: LDA >= fla_max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);

  // IPIV is INTEGER array, dimension (N)
  integer *ipivbuff, *ipivrefbuff;
  allocate_init_buffer(ipivbuff, ipivrefbuff, n, 0);
  
  // IPIV2 is INTEGER array, dimension (N)
  integer *ipiv2buff, *ipiv2refbuff;
  allocate_init_buffer(ipiv2buff, ipiv2refbuff, n, 0);
  
  /* LDB is INTEGER
          The leading dimension of the array B.  LDB >= fla_max(1,N).*/
  integer ldb = lin_driver_paramslist[ip].ldb;
  if (ldb < fla_max(1, n)) {
    PRINTF("ldb < fla_max(1, n) but it should be: LDB >= fla_max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // B is COMPLEX or COMPLEX16 array, dimension (LDB,NRHS)
  T *bbuff, *brefbuff;
  allocate_init_buffer(bbuff, brefbuff, ldb * nrhs);
  
  /* LTB is INTEGER
          The size of the array TB. LTB >= 4*N, internally
          used to select NB such that LTB >= (3*NB+1)*N.*/
  integer ltb = lin_driver_paramslist[ip].ltb;
  integer ltb_size = ltb;
  if ((ltb < -1) || (ltb == 0)) {
    PRINTF("ltb is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  
  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = lin_driver_paramslist[ip].lwork_hesv_aa_2stage;
  integer lwork_size = lwork;
  
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  
  integer info_cpp = -1;
  
  /* If lwork/lrwork/liwork is -1, then call this API with default work
     variables. In return, work variables will be updated with array sizes
     needed.*/
  if ((lwork == -1) || (ltb == -1)) {
    PRINTF("lwork or ltb is -1, so call hesv_aa_2stage() to get the array" \
           " sizes.\n");
    T worksize = {0};
    T tbsize = {0};
    
    // Call CPP function
    libflame::hesv_aa_2stage<T>(&uplo, &n, &nrhs, abuff, &lda, &tbsize, &ltb,
                ipivbuff, ipiv2buff, bbuff, &ldb, &worksize, &lwork,
                &info_cpp);
    PRINTF("info_cpp: %d, worksize: %f, tbsize: %f\n", info_cpp, \
            worksize.real, tbsize.real);
    if (info_cpp == 0) {
      if (lwork == -1) {
        lwork_size = worksize.real;
      }
      if (ltb == -1) {
        ltb_size = tbsize.real;
      }
    }
  }
  
  if (ltb_size < 4*n) {
    PRINTF("ltb < 4*n but it should be: ltb >= 4*n. Please " \
           "correct the input data.\n");
  }
  
  // TB is COMPLEX or COMPLEX16 array, dimension (LTB)
  T *tbbuff, *tbrefbuff;
  allocate_init_buffer(tbbuff, tbrefbuff, ltb_size, 0);
  
  if (lwork_size < n) {
    PRINTF("lwork < n but it should be: lwork >= n. Please " \
           "correct the input data.\n");
  }
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (LWORK)
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, lwork_size, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("nrhs = %d\n", nrhs);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("ltb = %d\n", ltb_size);
    PRINTF("Size of TB array (ltb) = %d\n", ltb);
    PRINTF("Size of IPIV array (n) = %d\n", n);
    PRINTF("Size of IPIV2 array (n) = %d\n", n);
    PRINTF("ldb = %d\n", ldb);
    PRINTF("Size of B array (ldb*nrhs) = %d\n", ldb * nrhs);
    PRINTF("Size of WORK array (LWORK) = %d\n", lwork_size);
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
    
    // Prints TB array contents
    strncpy(arrayname, "TB input", arraysize);
    print_array<T>(arrayname, tbbuff, ltb_size);
    strncpy(arrayname, "TB ref input", arraysize);
    print_array<T>(arrayname, tbrefbuff, ltb_size);
    
    // Prints IPIV array contents
    strncpy(arrayname, "IPIV input", arraysize);
    print_array<integer>(arrayname, ipivbuff, n);
    strncpy(arrayname, "IPIV ref input", arraysize);
    print_array<integer>(arrayname, ipivrefbuff, n);
    
    // Prints IPIV2 array contents
    strncpy(arrayname, "IPIV2 input", arraysize);
    print_array<integer>(arrayname, ipiv2buff, n);
    strncpy(arrayname, "IPIV2 ref input", arraysize);
    print_array<integer>(arrayname, ipiv2refbuff, n);
    
    // Prints B array contents
    strncpy(arrayname, "B input", arraysize);
    print_array<T>(arrayname, bbuff, ldb * nrhs);
    strncpy(arrayname, "B ref input", arraysize);
    print_array<T>(arrayname, brefbuff, ldb * nrhs);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, lwork_size);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, lwork_size);
  #endif

  // Call CPP function
  info_cpp = -1;
  libflame::hesv_aa_2stage<T>(&uplo, &n, &nrhs, abuff, &lda, tbbuff, &ltb_size,
              ipivbuff, ipiv2buff, bbuff, &ldb, workbuff, &lwork_size,
              &info_cpp);
  
  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    hesv_aa_2stage_ref = (fptr_NL_LAPACK_hesv_aa_2stage)dlsym(lapackModule, \
                          "chesv_aa_2stage_");
  } else if (typeid(T) == typeid(dcomplex)) {
    hesv_aa_2stage_ref = (fptr_NL_LAPACK_hesv_aa_2stage)dlsym(lapackModule, \
                          "zhesv_aa_2stage_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (hesv_aa_2stage_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  
  integer info_ref = -1;
  hesv_aa_2stage_ref(&uplo, &n, &nrhs, arefbuff, &lda, tbrefbuff, &ltb_size,
      ipivrefbuff, ipiv2refbuff, brefbuff, &ldb, workrefbuff, &lwork_size,
      &info_ref);
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
      
      // Prints TB array contents
      strncpy(arrayname, "TB output", arraysize);
      print_array<T>(arrayname, tbbuff, ltb_size);
      strncpy(arrayname, "TB ref output", arraysize);
      print_array<T>(arrayname, tbrefbuff, ltb_size);
      
      // Prints IPIV array contents
      strncpy(arrayname, "IPIV output", arraysize);
      print_array<integer>(arrayname, ipivbuff, n);
      strncpy(arrayname, "IPIV ref output", arraysize);
      print_array<integer>(arrayname, ipivrefbuff, n);
      
      // Prints IPIV2 array contents
      strncpy(arrayname, "IPIV2 output", arraysize);
      print_array<integer>(arrayname, ipiv2buff, n);
      strncpy(arrayname, "IPIV2 ref output", arraysize);
      print_array<integer>(arrayname, ipiv2refbuff, n);
      
      // Prints B array contents
      strncpy(arrayname, "B output", arraysize);
      print_array<T>(arrayname, bbuff, ldb * nrhs);
      strncpy(arrayname, "B ref output", arraysize);
      print_array<T>(arrayname, brefbuff, ldb * nrhs);
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, lwork_size);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, lwork_size);
    #endif
    
    double diff = computeError<T>(lda, n, arefbuff, abuff);
    diff += computeError<T>(1, ltb_size, tbrefbuff, tbbuff);
    diff += computeError<integer>(1, n, ipivrefbuff, ipivbuff);
    diff += computeError<integer>(1, n, ipiv2refbuff, ipiv2buff);
    diff += computeError<T>(ldb, nrhs, brefbuff, bbuff);
    diff += computeError<T>(1, lwork_size, workrefbuff, workbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, abs(diff), LIN_DRVR_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers.
  delete[] abuff; delete[] arefbuff;
  delete[] tbbuff; delete[] tbrefbuff;
  delete[] bbuff; delete[] brefbuff;
  delete[] ipivbuff; delete[] ipivrefbuff;
  delete[] ipiv2buff; delete[] ipiv2refbuff;
  delete[] workbuff; delete[] workrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_hesv_aa_2stage, CHESV_AA_2STAGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hesv_aa_2stage_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_hesv_aa_2stage, ZHESV_AA_2STAGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hesv_aa_2stage_test<dcomplex> (index);
  }
}