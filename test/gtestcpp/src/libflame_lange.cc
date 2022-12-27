/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_lange.cc
 *  @brief Test application to validate lange() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  lange_test is function template for lange() functions.
			T can be float, double
 * @details
 * \b Purpose:
    \verbatim
	  lange_test is function template for lange() functions.
	  T can be float, double
	  
	  lange_test() function template calls C and CPP based library APIs with
	  valid test values, calculate the differences in output.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold.

    Real Reference:
    http://www.netlib.org/lapack/explore-html/da/d23/group__real_g_eauxiliary_ga459d27829607393670ef7de8a6914933.html#ga459d27829607393670ef7de8a6914933
    Double Reference:
    http://www.netlib.org/lapack/explore-html/de/d39/group__double_g_eauxiliary_gaefa80dbd8cd1732740478618b8b622a1.html#gaefa80dbd8cd1732740478618b8b622a1
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void lange_test(int ip)
{
  typedef T (*fptr_NL_LAPACK_lange)(char* norm, integer* m,
                      integer* n, T* a, integer* lda, T* work);
  fptr_NL_LAPACK_lange lange_ref = NULL;
  T value = 0, valueref = 0;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* NORM is CHARACTER*1
          Specifies whether the 1-norm condition number or the
          infinity-norm condition number is required:
          CLANSY = ( fla_max(abs(A(i,j))), NORM = 'M' or 'm'
             (
             ( norm1(A),         NORM = '1', 'O' or 'o'
             (
             ( normI(A),         NORM = 'I' or 'i'
             (
             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'*/
  char norm = eig_paramslist[ip].norm;
  if ((norm != '1') && (norm != 'O') && (norm != 'I')) {
    PRINTF("norm should be 1 or O or I. Please correct the input data.\n");
  }
  
  /* M is INTEGER
          The number of rows of the matrix A.  M >= 0.  When M = 0,
          ZLANGE is set to zero.*/
  integer m = eig_paramslist[ip].m;
  if (m < 0) {
    PRINTF("m < 0 but it should be: m >= 0. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= fla_max(M,1).*/
  integer lda = eig_paramslist[ip].lda_lange;
  if (lda < fla_max(m, 1)) {
    PRINTF("lda < fla_max(m, 1) but it should be: LDA >= fla_max(M,1). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  /* LWORK >= M when NORM = 'I'; otherwise, WORK is not
          referenced.*/
  integer lwork = eig_paramslist[ip].lwork_lange;
  
  // WORK is REAL or DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  if (norm == 'I') {
    if (lwork >= m) {
      allocate_init_buffer(workbuff, workrefbuff, fla_max(1, lwork), 0);
    } else {
      PRINTF("lwork < m but it should be: lwork >= m. Please " \
             "correct the input data.\n");
    }
  }
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("norm = %c\n", norm);
    PRINTF("m = %d\n", m);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("lwork = %d\n", lwork);
    PRINTF("Size of WORK array (MAX(1,LWORK)) = %d\n", fla_max(1, lwork));
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
    
    if ((workbuff != NULL) && (workrefbuff != NULL)) {
      // Prints WORK array contents
      strncpy(arrayname, "WORK input", arraysize);
      print_array<T>(arrayname, workbuff, fla_max(1, lwork));
      strncpy(arrayname, "WORK ref input", arraysize);
      print_array<T>(arrayname, workrefbuff, fla_max(1, lwork));
    }
  #endif
  
  // Call CPP function
  value = libflame::lange<T>(&norm, &m, &n, abuff, &lda, workbuff);
  
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(float)) {
    lange_ref = (fptr_NL_LAPACK_lange)dlsym(lapackModule, "slange_");
  } else if (typeid(T) == typeid(double)) {
    lange_ref = (fptr_NL_LAPACK_lange)dlsym(lapackModule, "dlange_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (lange_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  valueref = lange_ref(&norm, &m, &n, arefbuff, &lda, workrefbuff);
  PRINTF("value = %lf, valueref = %lf\n", value, valueref);
  
  // Calculate the differences of buffers.
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all Output arrays contents...\n");
    if ((workbuff != NULL) && (workrefbuff != NULL)) {
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, fla_max(1, lwork));
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, fla_max(1, lwork));
    }
  #endif
  
  double diff = computeError<T>(1, 1, &value, &valueref);
  if ((workbuff != NULL) && (workrefbuff != NULL)) {
    diff = computeError<T>(1, fla_max(1, lwork), workbuff, workrefbuff);
  }
  PRINTF("diff: %lf\n", diff);
  EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  if ((workbuff != NULL) && (workrefbuff != NULL)) {
    delete[] workbuff; delete[] workrefbuff;
  }
}

/*! @brief  lange_test_cmplx is function template for lange() functions.
			T can be scomplex, dcomplex
      Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  lange_test_cmplx is function template for lange() functions.
	  T can be scomplex, dcomplex
    Ta can be float, double.
	  
	  lange_test_cmplx() function template calls C and CPP based library APIs
	  with valid test values, calculate the differences in output.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold.
    
    Complex reference:
	  http://www.netlib.org/lapack/explore-html/d3/d71/group__complex_g_eauxiliary_gaa4e1d57c726257bbbfe0c89ef5461c3b.html#gaa4e1d57c726257bbbfe0c89ef5461c3b
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d0/d9e/group__complex16_g_eauxiliary_ga7908bb12a6f02dbfa4d5a92a27c0e9b7.html#ga7908bb12a6f02dbfa4d5a92a27c0e9b7
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T, typename Ta>
void lange_test_cmplx(int ip)
{
  typedef Ta (*fptr_NL_LAPACK_lange)(char* norm, integer* m, integer* n,
                      T* a, integer* lda, Ta* work);
  fptr_NL_LAPACK_lange lange_ref = NULL;
  Ta value = 0, valueref = 0;
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* NORM is CHARACTER*1
          Specifies whether the 1-norm condition number or the
          infinity-norm condition number is required:
          CLANSY = ( fla_max(abs(A(i,j))), NORM = 'M' or 'm'
             (
             ( norm1(A),         NORM = '1', 'O' or 'o'
             (
             ( normI(A),         NORM = 'I' or 'i'
             (
             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'*/
  char norm = eig_paramslist[ip].norm;
  if ((norm != '1') && (norm != 'O') && (norm != 'I')) {
    PRINTF("norm should be 1 or O or I. Please correct the input data.\n");
  }
  
  /* M is INTEGER
          The number of rows of the matrix A.  M >= 0.  When M = 0,
          ZLANGE is set to zero.*/
  integer m = eig_paramslist[ip].m;
  if (m < 0) {
    PRINTF("m < 0 but it should be: m >= 0. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= fla_max(M,1).*/
  integer lda = eig_paramslist[ip].lda_lange;
  if (lda < fla_max(m, 1)) {
    PRINTF("lda < fla_max(m, 1) but it should be: LDA >= fla_max(M,1). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  /* LWORK >= N when NORM = 'I' or '1' or 'O';otherwise,
                          WORK is not referenced.*/
  integer lwork = eig_paramslist[ip].lwork_lange;
  
  // WORK is REAL or DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  Ta *workbuff = NULL, *workrefbuff = NULL;
  if (norm == 'I') {
    if (lwork >= m) {
      allocate_init_buffer(workbuff, workrefbuff, fla_max(1, lwork), 0);
    } else {
      PRINTF("lwork < m but it should be: lwork >= m. Please " \
             "correct the input data.\n");
    }
  }
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("norm = %c\n", norm);
    PRINTF("m = %d\n", m);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("lwork = %d\n", lwork);
    PRINTF("Size of WORK array (MAX(1,LWORK)) = %d\n", fla_max(1, lwork));
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
    if ((workbuff != NULL) && (workrefbuff != NULL)) {
      // Prints WORK array contents
      strncpy(arrayname, "WORK input", arraysize);
      print_array<Ta>(arrayname, workbuff, fla_max(1, lwork));
      strncpy(arrayname, "WORK ref input", arraysize);
      print_array<Ta>(arrayname, workrefbuff, fla_max(1, lwork));
    }
  #endif
  
  // Call CPP function
  value = libflame::lange<T, Ta>(&norm, &m, &n, abuff, &lda, workbuff);
  
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    lange_ref = (fptr_NL_LAPACK_lange)dlsym(lapackModule, "clange_");
  } else if (typeid(T) == typeid(dcomplex)) {
    lange_ref = (fptr_NL_LAPACK_lange)dlsym(lapackModule, "zlange_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (lange_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  valueref = lange_ref(&norm, &m, &n, arefbuff, &lda, workrefbuff);
  PRINTF("value = %lf, valueref = %lf\n", value, valueref);
  
  // Calculate the differences of buffers.
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all Output arrays contents...\n");
    
    if ((workbuff != NULL) && (workrefbuff != NULL)) {
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<Ta>(arrayname, workbuff, fla_max(1, lwork));
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<Ta>(arrayname, workrefbuff, fla_max(1, lwork));
    }
  #endif
  
  double diff = computeError<Ta>(1, 1, &value, &valueref);
  if ((workbuff != NULL) && (workrefbuff != NULL)) {
    diff = computeError<Ta>(1, fla_max(1, lwork), workbuff, workrefbuff);
  }
  PRINTF("diff: %lf\n", diff);
  EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  if ((workbuff != NULL) && (workrefbuff != NULL)) {
    delete[] workbuff; delete[] workrefbuff;
  }
}

/* Use TEST macro and call C++ test function template with
   float as typenames.*/
TEST(LAPACKCPP_lange, SLANGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    lange_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double as typenames.*/
TEST(LAPACKCPP_lange, DLANGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    lange_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_lange, CLANGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    lange_test_cmplx<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_lange, ZLANGE) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    lange_test_cmplx<dcomplex, double> (index);
  }
}