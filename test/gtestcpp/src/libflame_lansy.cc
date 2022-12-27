/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_lansy.cc
 *  @brief Test application to validate lansy() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  lansy_test is function template for lansy() functions.
			T can be float, double
 * @details
 * \b Purpose:
    \verbatim
	  lansy_test is function template for lansy() functions.
	  T can be float, double
	  
	  lansy_test() function template calls C and CPP based library APIs with
	  valid test values, calculate the differences in output.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold.

    Real Reference:
    http://www.netlib.org/lapack/explore-html/de/d0d/group__real_s_yauxiliary_ga611e1beaaad792e0753a47723c8380ed.html#ga611e1beaaad792e0753a47723c8380ed
    Double Reference:
    http://www.netlib.org/lapack/explore-html/d9/df5/group__double_s_yauxiliary_ga8e0d957efd6f93764d9bc98a7aa1927a.html#ga8e0d957efd6f93764d9bc98a7aa1927a
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void lansy_test(int ip)
{
  typedef T (*fptr_NL_LAPACK_lansy)(char* norm, char* uplo,
                      integer* n, T* a, integer* lda, T* work);
  fptr_NL_LAPACK_lansy lansy_ref = NULL;
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
  
   /* UPLO is CHARACTER*1
          Specifies whether the upper or lower triangular part of the
          Hermitian matrix A is stored, and how B has been factorized.
          = 'U':  Upper triangular
          = 'L':  Lower triangular*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= fla_max(1,N).*/
  integer lda = eig_paramslist[ip].lda;
  if (lda < fla_max(1, n)) {
    PRINTF("lda < fla_max(1, n) but it should be: LDA >= fla_max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  /* LWORK >= N when NORM = 'I' or '1' or 'O';otherwise,
                          WORK is not referenced.*/
  integer lwork = eig_paramslist[ip].lwork_lansy;
  
  // WORK is REAL or DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  if ((norm == '1') or (norm == 'I') || (norm == 'O')) {
    if (lwork >= n) {
      allocate_init_buffer(workbuff, workrefbuff, fla_max(1, lwork), 0);
    } else {
      PRINTF("lwork < n but it should be: lwork >= n. Please " \
             "correct the input data.\n");
    }
  }
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("norm = %c\n", norm);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
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
  value = libflame::lansy<T>(&norm, &uplo, &n, abuff, &lda, workbuff);
  
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(float)) {
    lansy_ref = (fptr_NL_LAPACK_lansy)dlsym(lapackModule, "slansy_");
  } else if (typeid(T) == typeid(double)) {
    lansy_ref = (fptr_NL_LAPACK_lansy)dlsym(lapackModule, "dlansy_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (lansy_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  valueref = lansy_ref(&norm, &uplo, &n, arefbuff, &lda, workrefbuff);
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
    diff += computeError<T>(1, fla_max(1, lwork), workbuff, workrefbuff);
  }
  PRINTF("diff: %lf\n", diff);
  EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  
  // Free up the buffers
  delete[] abuff; delete[] arefbuff;
  if ((workbuff != NULL) && (workrefbuff != NULL)) {
    delete[] workbuff; delete[] workrefbuff;
  }
}
/*! @brief  lansy_test_cmplx is function template for lansy() functions.
			T can be scomplex, dcomplex
      Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  lansy_test_cmplx is function template for lansy() functions.
	  T can be scomplex, dcomplex
    Ta can be float, double.
	  
	  lansy_test_cmplx() function template calls C and CPP based library APIs
	  with valid test values, calculate the differences in output.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold.
    
    Complex reference:
	  http://www.netlib.org/lapack/explore-html/d2/d15/group__complex_s_yauxiliary_gad2c86a28190eb12c91cda1c4faef5df7.html#gad2c86a28190eb12c91cda1c4faef5df7
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/de/d14/group__complex16_s_yauxiliary_gae1d67e9c7403f3d6e2c5db6073b014d3.html#gae1d67e9c7403f3d6e2c5db6073b014d3
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T, typename Ta>
void lansy_test_cmplx(int ip)
{
  typedef Ta (*fptr_NL_LAPACK_lansy)(char* norm, char* uplo,
                      integer* n, T* a, integer* lda, Ta* work);
  fptr_NL_LAPACK_lansy lansy_ref = NULL;
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
  
   /* UPLO is CHARACTER*1
          Specifies whether the upper or lower triangular part of the
          Hermitian matrix A is stored, and how B has been factorized.
          = 'U':  Upper triangular
          = 'L':  Lower triangular*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but it should be: n >= 0. Please correct the input data.\n");
  }
  
  /* LDA is INTEGER
          The leading dimension of the array A.  LDA >= fla_max(1,N).*/
  integer lda = eig_paramslist[ip].lda;
  if (lda < fla_max(1, n)) {
    PRINTF("lda < fla_max(1, n) but it should be: LDA >= fla_max(1,N). Please " \
           "correct the input data.\n");
  }
  
  // A is COMPLEX or COMPLEX*16 array, dimension (LDA,N)
  T *abuff, *arefbuff;
  allocate_init_buffer(abuff, arefbuff, lda * n);
  
  /* LWORK >= N when NORM = 'I' or '1' or 'O';otherwise,
                          WORK is not referenced.*/
  integer lwork = eig_paramslist[ip].lwork_lansy;
  
  // WORK is REAL or DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  Ta *workbuff = NULL, *workrefbuff = NULL;
  if ((norm == '1') or (norm == 'I') || (norm == 'O')) {
    if (lwork >= n) {
      allocate_init_buffer(workbuff, workrefbuff, fla_max(1, lwork), 0);
    } else {
      PRINTF("lwork < n but it should be: lwork >= n. Please " \
             "correct the input data.\n");
    }
  }
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("norm = %c\n", norm);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
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
  value = libflame::lansy<T, Ta>(&norm, &uplo, &n, abuff, &lda, workbuff);
  
  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    lansy_ref = (fptr_NL_LAPACK_lansy)dlsym(lapackModule, "clansy_");
  } else if (typeid(T) == typeid(dcomplex)) {
    lansy_ref = (fptr_NL_LAPACK_lansy)dlsym(lapackModule, "zlansy_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  if (lansy_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }

  valueref = lansy_ref(&norm, &uplo, &n, arefbuff, &lda, workrefbuff);
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
    diff += computeError<Ta>(1, fla_max(1, lwork), workbuff, workrefbuff);
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
TEST(LAPACKCPP_lansy, SLANSY) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    lansy_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double as typenames.*/
TEST(LAPACKCPP_lansy, DLANSY) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    lansy_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_lansy, CLANSY) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    lansy_test_cmplx<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_lansy, ZLANSY) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    lansy_test_cmplx<dcomplex, double> (index);
  }
}