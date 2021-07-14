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
	  valid test values and returns the differences in output.
	  
	  Complex reference:
	  http://www.netlib.org/lapack/explore-html/d3/db9/group__complex_o_t_h_e_rcomputational_ga7de86c95768cba8a2168ee787f18f9f4.html#ga7de86c95768cba8a2168ee787f18f9f4
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/d0/da6/group__complex16_o_t_h_e_rcomputational_gae10651c17f5235233e41c53bfc4f9f93.html#gae10651c17f5235233e41c53bfc4f9f93
    \endverbatim
	
 * @param[in] IP
      IP is INTEGER
		  Used to pass Index of Eigen Parameters array present in config file.

 * @return DOUBLE
      Returns Differences value after comparing output of C and CPP based
		  library APIs.
 * */
template< typename T, typename Ta >
double hbtrd_test(int ip)
{
  typedef integer (*Fptr_NL_LAPACKE_hbtrd)(char* vect, char* uplo, integer* n,
                        integer* kd, T* ab, integer* ldab, Ta* d, Ta* e, T* q,
                        integer* ldq, T* work, integer* info);
  Fptr_NL_LAPACKE_hbtrd HBTRD = NULL;
  
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
  allocate_init_buffer(dbuff, drefbuff, n);
  
  // E is REAL or DOUBLE PRECISION array, dimension (N-1)
  Ta *ebuff, *erefbuff;
  allocate_init_buffer(ebuff, erefbuff, n-1);
  
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
  allocate_init_buffer(workbuff, workrefbuff, n);
  
  // Call CPP function
  integer info_cpp = libflame::hbtrd<T, Ta>(&vect, &uplo, &n, &kd, abbuff,
                        &ldab, dbuff, ebuff, qbuff, &ldq, workbuff, &info_cpp);

  // Call C function
  
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    HBTRD = (Fptr_NL_LAPACKE_hbtrd)dlsym(lapackModule, "chbtrd_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HBTRD = (Fptr_NL_LAPACKE_hbtrd)dlsym(lapackModule, "zhbtrd_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HBTRD == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }
  
  integer info_ref = -1;
  HBTRD(&vect, &uplo, &n, &kd, abrefbuff, &ldab, drefbuff, erefbuff,
    qrefbuff, &ldq, workrefbuff, &info_ref);
  PRINTF ("info_cpp: %u, info_ref: %u\n", info_cpp, info_ref);
  
  // Calculate the differences of buffers.
  double diff = 0.0;
  if ((info_cpp == 0) && (info_ref == 0)) {
    diff =  computeError<T>(ldab, n, abrefbuff, abbuff);
    diff +=  computeError<Ta>(1, n, dbuff, drefbuff);
    diff +=  computeError<Ta>(1, n, ebuff, erefbuff);
    if (vect != 'N') {
      diff +=  computeError<T>(ldq, n, qbuff, qrefbuff);
    }
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
  }
  
  // Free up the buffers
  delete[] abbuff; delete[] abrefbuff;
  delete[] dbuff; delete[] drefbuff;
  delete[] ebuff; delete[] erefbuff;
  delete[] qbuff; delete[] qrefbuff;
  
  // Return the difference.
  return abs(diff);
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_hbtrd, CHBTRD) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    diff = hbtrd_test<scomplex, float> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
    PRINTF("diff: %lf\n", diff);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_hbtrd, ZHBTRD) {
  double diff = 0.0;
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    diff = hbtrd_test<dcomplex, double> (index);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
    PRINTF("diff: %lf\n", diff);
  }
}