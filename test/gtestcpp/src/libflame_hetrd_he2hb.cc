/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_hetrd_he2hb.cc
 *  @brief Test application to validate hetrd_he2hb() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hetrd_he2hb_test is function template for hetrd_he2hb() functions.
			T can be scomplex, dcomplex
      Ta can be float, double
 * @details
 * \b Purpose:
    \verbatim
	  hetrd_he2hb_test is function template for hetrd_he2hb() functions.
	  T can be scomplex, dcomplex
    Ta can be float, double
	  
	  hetrd_he2hb_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO = 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d4/d74/group__complex_h_ecomputational_gad8c7862093b3ac5727a6e2a3b1df1b73.html#gad8c7862093b3ac5727a6e2a3b1df1b73
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d0/da6/group__complex16_o_t_h_e_rcomputational_gacb962c8889c699b70126db722d2c05d0.html#gacb962c8889c699b70126db722d2c05d0
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void hetrd_he2hb_test(int ip)
{
  typedef int (*Fptr_NL_LAPACK_hetrd_he2hb)(char *uplo, integer *n,
                    integer *kd, T *a, integer *lda, T *ab,
                    integer *ldab, T *tau, T *work,
                    integer *lwork, integer *info);
  Fptr_NL_LAPACK_hetrd_he2hb HETRD_HE2HB;
  
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
  
  /* KD is INTEGER
          The number of superdiagonals of the reduced matrix if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
          The reduced matrix is stored in the array AB.*/
  integer kd = eig_paramslist[ip].kd;
  if (kd < 0) {
    PRINTF("kd < 0 but should be: kd >= 0. Please correct the input data.");
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
  
  // LDAB is INTEGER. The leading dimension of the array AB.  LDAB >= KD+1.
  integer ldab = eig_paramslist[ip].ldab;
  if (ldab < (kd + 1)) {
    PRINTF("ldab < (kd+1) but it should be: LDAB >= (kd+1). Please " \
           "correct the input data.\n");
  }
  
  // AB is COMPLEX or COMPLEX*16 array, dimension (LDAB, N)
  T *abbuff, *abrefbuff;
  allocate_init_buffer(abbuff, abrefbuff, ldab * n, 0);
  
  // TAU is COMPLEX or COMPLEX*16 array, dimension (N-KD)
  T *taubuff, *taurefbuff;
  allocate_init_buffer(taubuff, taurefbuff, n-kd, 0);
  
  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = eig_paramslist[ip].lwork_hbtrd_he2hb;
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
    PRINTF("lwork is -1, so call hetrd_2stage() to get the array sizes.\n");
    T worksize = {0};
    
    // Call CPP function
    libflame::hetrd_he2hb<T>(&uplo, &n, &kd, abuff, &lda, abbuff, &ldab,
                        taubuff, &worksize, &lwork_size, &info_cpp);
    PRINTF("info_cpp: %d, worksize: %f\n", info_cpp, worksize.real);
    if (info_cpp == 0) {
      if (lwork == -1) {
        lwork_size = worksize.real;
      }
    }
  }

  // WORK is COMPLEX or COMPLEX*16  array, dimension (LWORK)
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, lwork_size, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("kd = %d\n", kd);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("ldab = %d\n", ldab);
    PRINTF("Size of AB array (ldab*n) = %d\n", ldab * n);
    PRINTF("Size of TAU array (n-kd) = %d\n", n-kd);
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
    
    // Prints AB array contents
    strncpy(arrayname, "AB input", arraysize);
    print_array<T>(arrayname, abbuff, ldab * n);
    strncpy(arrayname, "AB ref input", arraysize);
    print_array<T>(arrayname, abrefbuff, ldab * n);
    
    // Prints TAU array contents
    strncpy(arrayname, "TAU input", arraysize);
    print_array<T>(arrayname, taubuff, n-kd);
    strncpy(arrayname, "TAU ref input", arraysize);
    print_array<T>(arrayname, taurefbuff, n-kd);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, lwork_size);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, lwork_size);
  #endif
  
  // Call CPP function
  info_cpp = -1;
  libflame::hetrd_he2hb<T>(&uplo, &n, &kd, abuff, &lda, abbuff, &ldab,
                        taubuff, workbuff, &lwork_size, &info_cpp);

  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    HETRD_HE2HB = (Fptr_NL_LAPACK_hetrd_he2hb)dlsym(lapackModule, \
                          "chetrd_he2hb_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HETRD_HE2HB = (Fptr_NL_LAPACK_hetrd_he2hb)dlsym(lapackModule, \
                          "zhetrd_he2hb_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HETRD_HE2HB == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  
  integer info_ref = -1;
  HETRD_HE2HB(&uplo, &n, &kd, arefbuff, &lda, abrefbuff, &ldab,
              taurefbuff, workrefbuff, &lwork_size, &info_ref);
  
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
      
      // Prints AB array contents
      strncpy(arrayname, "AB output", arraysize);
      print_array<T>(arrayname, abbuff, ldab * n);
      strncpy(arrayname, "AB ref output", arraysize);
      print_array<T>(arrayname, abrefbuff, ldab * n);
      
      // Prints TAU array contents
      strncpy(arrayname, "TAU output", arraysize);
      print_array<T>(arrayname, taubuff, n-kd);
      strncpy(arrayname, "TAU ref output", arraysize);
      print_array<T>(arrayname, taurefbuff, n-kd);
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, lwork_size);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, lwork_size);
    #endif
  
    double diff = computeError<T>(lda, n, abuff, arefbuff);
    diff += computeError<T>(ldab, n, abbuff, abrefbuff);
    diff += computeError<T>(1, n-kd, taubuff, taurefbuff);
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
  delete[] abbuff; delete[] abrefbuff;
  delete[] taubuff; delete[] taurefbuff;
  delete[] workbuff; delete[] workrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex, float as typenames.*/
TEST(LAPACKCPP_hetrd_he2hb, CHETRD_HE2HB) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetrd_he2hb_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex, double as typenames.*/
TEST(LAPACKCPP_hetrd_he2hb, ZHETRD_HE2HB) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetrd_he2hb_test<dcomplex> (index);
  }
}