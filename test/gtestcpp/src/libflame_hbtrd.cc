/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hbtrd.cc
 *  @brief Test application to validate hbtrd() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hbtrd_test is function template for hbtrd() functions.
			T can be scomplex, dcomplex
			Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hbtrd_test is function template for hbtrd() functions.
	  T can be scomplex, dcomplex
	  Ta can be float, double.
	  
	  hbtrd_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO is >= 0.
    And passses the test case if difference is <= threshold.
    Fails the test case if difference > threshold or INFO < 0.
	  
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d3/db9/group__complex_o_t_h_e_rcomputational_ga7de86c95768cba8a2168ee787f18f9f4.html#ga7de86c95768cba8a2168ee787f18f9f4
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d0/da6/group__complex16_o_t_h_e_rcomputational_gae10651c17f5235233e41c53bfc4f9f93.html#gae10651c17f5235233e41c53bfc4f9f93
    \endverbatim
	
 * @param[in] IP
      IP is INTEGER
		  Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template< typename T, typename Ta >
void hbtrd_test(int ip)
{
  typedef integer (*fptr_NL_LAPACK_hbtrd)(char* vect, char* uplo, integer* n,
                        integer* kd, T* ab, integer* ldab, Ta* d, Ta* e, T* q,
                        integer* ldq, T* work, integer* info);
  fptr_NL_LAPACK_hbtrd hbtrd_ref = NULL;
  
  // Initialise random number generators with timestamp
  srand (time(NULL));
  
  /* VECT is CHARACTER*1
          = 'N':  do not form Q;
          = 'V':  form Q;
          = 'U':  update a matrix X, by forming X*Q.*/
  char vect = eig_paramslist[ip].vect;
  if ((vect != 'N') && (vect != 'V') && (vect != 'U')) {
    PRINTF("vect should be N or V or U. Please correct the input data.\n");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("jobz should be N or V. Please correct the input data.");
  }
  
  /* N is INTEGER
          The order of the matrices A.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.");
  }
  
  /* KD is INTEGER
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.*/
  integer kd = 0;
  if (uplo == 'U') {
	  kd = eig_paramslist[ip].sda;
  } else if (uplo == 'L') {
  	kd = eig_paramslist[ip].subda;
  }
  if (kd < 0) {
    PRINTF("kd is 0 but should be: KD >= 0. Please correct the input data.");
  }
  
  /* LDAB is INTEGER
          The leading dimension of the array AB.  LDAB >= KD+1.*/
  integer ldab = eig_paramslist[ip].ldab;
  if (ldab < (kd+1)) {
    PRINTF("ldab < (kd+1) but it should be: LDAB >= KD + 1. Please correct" \
          " the input data.\n");
  }
  
  // AB is COMPLEX or COMPLEX*16 array, dimension (LDAB,N)
  T *abbuff, *abrefbuff;
  allocate_init_buffer(abbuff, abrefbuff, ldab * n);
  
  // D is REAL or DOUBLE PRECISION array, dimension (N)
  Ta *dbuff, *drefbuff;
  allocate_init_buffer(dbuff, drefbuff, n, 0);
  
  // E is REAL or DOUBLE PRECISION array, dimension (N-1)
  Ta *ebuff, *erefbuff;
  allocate_init_buffer(ebuff, erefbuff, n-1, 0);
  
  /* LDQ is INTEGER
          The leading dimension of the array Q.
          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.*/
  integer ldq = eig_paramslist[ip].ldq;
  if (ldq < 1) {
    PRINTF("ldq < 1 but should be: LDQ >= 1. Please correct the input data.");
  }
  if ((vect == 'V') || (vect == 'U')) {
    if (ldq < n) {
      PRINTF("When vect is V or U, ldq < n but it should be: LDQ >= N." \
              "Please correct the input data.\n");
    }
  }
  
  // Q is COMPLEX or COMPLEX*16 array, dimension (LDQ,N)
  T *qbuff, *qrefbuff; // output buffer
  allocate_init_buffer(qbuff, qrefbuff, ldq * n);
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (N)
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, n, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("vect = %c\n", vect);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("kd = %d\n", kd);
    PRINTF("ldab = %d\n", ldab);
    PRINTF("Size of AB array (ldab*n) = %d\n", ldab*n);
    PRINTF("Size of D array (n) = %d\n", n);
    PRINTF("Size of E array (n) = %d\n", n-1);
    PRINTF("ldq = %d\n", ldq);
    PRINTF("Size of Q array (ldq*n) = %d\n", ldq*n);
    PRINTF("Size of WORK array (n) = %d\n", n);
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
    
    // Prints AB array contents
    strncpy(arrayname, "AB input", arraysize);
    print_array<T>(arrayname, abbuff, ldab * n);
    strncpy(arrayname, "AB ref input", arraysize);
    print_array<T>(arrayname, abrefbuff, ldab * n);
    
    // Prints D array contents
    strncpy(arrayname, "D input", arraysize);
    print_array<Ta>(arrayname, dbuff, n);
    strncpy(arrayname, "D ref input", arraysize);
    print_array<Ta>(arrayname, drefbuff, n);
    
    // Prints E array contents
    strncpy(arrayname, "E input", arraysize);
    print_array<Ta>(arrayname, ebuff, n);
    strncpy(arrayname, "E ref input", arraysize);
    print_array<Ta>(arrayname, erefbuff, n);
    
    // Prints Q array contents
    strncpy(arrayname, "Q input", arraysize);
    print_array<T>(arrayname, qbuff, ldq * n);
    strncpy(arrayname, "Q ref input", arraysize);
    print_array<T>(arrayname, qrefbuff, ldq * n);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, n);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, n);
    #endif
  
  // Call CPP function
  integer info_cpp = libflame::hbtrd<T, Ta>(&vect, &uplo, &n, &kd, abbuff,
                        &ldab, dbuff, ebuff, qbuff, &ldq, workbuff, &info_cpp);

  // Call C function
  
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    hbtrd_ref = (fptr_NL_LAPACK_hbtrd)dlsym(lapackModule, "chbtrd_");
  } else if (typeid(T) == typeid(dcomplex)) {
    hbtrd_ref = (fptr_NL_LAPACK_hbtrd)dlsym(lapackModule, "zhbtrd_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (hbtrd_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }
  
  integer info_ref = -1;
  hbtrd_ref(&vect, &uplo, &n, &kd, abrefbuff, &ldab, drefbuff, erefbuff,
    qrefbuff, &ldq, workrefbuff, &info_ref);
  PRINTF ("info_cpp: %d, info_ref: %d\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  if ((info_cpp >= 0) && (info_ref >= 0)) {
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
        defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
      // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
      PRINTF("\nPrinting all Output arrays contents...\n");
      
      // Prints AB array contents
      strncpy(arrayname, "AB output", arraysize);
      print_array<T>(arrayname, abbuff, ldab * n);
      strncpy(arrayname, "AB ref output", arraysize);
      print_array<T>(arrayname, abrefbuff, ldab * n);
      
      // Prints D array contents
      strncpy(arrayname, "D output", arraysize);
      print_array<Ta>(arrayname, dbuff, n);
      strncpy(arrayname, "D ref output", arraysize);
      print_array<Ta>(arrayname, drefbuff, n);
      
      // Prints E array contents
      strncpy(arrayname, "E output", arraysize);
      print_array<Ta>(arrayname, ebuff, n);
      strncpy(arrayname, "E ref output", arraysize);
      print_array<Ta>(arrayname, erefbuff, n);
      
      // Prints Q array contents
      strncpy(arrayname, "Q output", arraysize);
      print_array<T>(arrayname, qbuff, ldq * n);
      strncpy(arrayname, "Q ref output", arraysize);
      print_array<T>(arrayname, qrefbuff, ldq * n);
      
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, n);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, n);
    #endif
    
    double diff = computeError<T>(ldab, n, abrefbuff, abbuff);
    diff += computeError<Ta>(1, n, dbuff, drefbuff);
    diff += computeError<Ta>(1, n, ebuff, erefbuff);
    if (vect != 'N') {
      diff += computeError<T>(ldq, n, qbuff, qrefbuff);
    }
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, abs(diff), SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abbuff; delete[] abrefbuff;
  delete[] dbuff; delete[] drefbuff;
  delete[] ebuff; delete[] erefbuff;
  delete[] qbuff; delete[] qrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_hbtrd, CHBTRD) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hbtrd_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_hbtrd, ZHBTRD) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hbtrd_test<dcomplex, double> (index);
  }
}