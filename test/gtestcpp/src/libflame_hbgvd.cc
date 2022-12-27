/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_hbgvd.cc
 *  @brief Test application to validate hbgvd() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hbgvd_test is function template for hbgvd() functions.
            T can be scomplex, dcomplex
            Ta can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  hbgvd_test is function template for hbgvd() functions.
	  T can be scomplex, dcomplex
    Ta can be float, double.
	  
	  hbgvd_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO >= 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.
	  
    Complex reference:
	  http://www.netlib.org/lapack/explore-html/d9/d98/group__complex_o_t_h_e_reigen_ga77b1c171ee971c0ff72107e4aa8b5376.html#ga77b1c171ee971c0ff72107e4aa8b5376
	  Complex double reference:
	  http://www.netlib.org/lapack/explore-html/db/d61/group__complex16_o_t_h_e_reigen_ga597ea234c22684386ad82c7515285514.html#ga597ea234c22684386ad82c7515285514
    \endverbatim
	
 * @param[in] IP
          IP is INTEGER
          Used to pass Index of Eigen Parameters array present in config file.

 * @return VOID
           Nothing.
 * */
template< typename T, typename Ta >
void hbgvd_test(int ip)
{
  typedef integer (*fptr_NL_LAPACK_hbgvd)(char* jobz, char* uplo, integer* n,
                      integer* ka, integer* kb, T* ab, integer* ldab, T* bb,
                      integer* ldbb, Ta* w, T* z, integer* ldz, T* work,
                      integer* lwork, Ta* rwork, integer* lrwork,
                      integer* iwork, integer* liwork, integer* info);
  fptr_NL_LAPACK_hbgvd hbgvd_ref = NULL;
  
  // Initialise random number generators with timestamp.
  srand (SRAND_SEED_VALUE);
  
  /* JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.*/
  char jobz = eig_paramslist[ip].jobz;
  if ((jobz != 'N') && (jobz != 'V')) {
    PRINTF("jobz should be N or V. Please correct the input data.\n");
  }
  
  /* UPLO is CHARACTER*1
          = 'U':  Upper triangles of A and B are stored;
          = 'L':  Lower triangles of A and B are stored.*/
  char uplo = eig_paramslist[ip].uplo;
  if ((uplo != 'U') && (uplo != 'L')) {
    PRINTF("uplo should be U or L. Please correct the input data.\n");
  }
  
  /* N is INTEGER
          The order of the matrices A and B.  N >= 0.*/
  integer n = eig_paramslist[ip].n;
  if (n < 0) {
    PRINTF("n < 0 but should be: n >= 0. Please correct the input data.\n");
  }

  /* KA is INTEGER
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'. KA >= 0.*/
  integer ka = 0;
  if (uplo == 'U') {
	  ka = eig_paramslist[ip].sda;
  } else if (uplo == 'L') {
	  ka = eig_paramslist[ip].subda;
  }
  if (ka < 0) {
    PRINTF("ka < 0 but should be: KA >= 0. Please correct the input data.\n");
  }
  
  /* KB is INTEGER
          The number of superdiagonals of the matrix B if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'. KB >= 0.*/
  integer kb = 0;
  if (uplo == 'U') {
	  kb = eig_paramslist[ip].sdb;
  } else if (uplo == 'L') {
	  kb = eig_paramslist[ip].subdb;
  }
  if (kb < 0) {
    PRINTF("kb < 0 but should be: KB >= 0. Please correct the input data.\n");
  }
  
  /* LDAB is INTEGER
          The leading dimension of the array AB.  LDAB >= KA+1.*/
  integer ldab = eig_paramslist[ip].ldab_hbgv;
  if (ldab < (ka + 1)) {
    PRINTF("ldab < (ka + 1) but should be: LDAB >= (KA + 1). Please correct" \
           " the input data.\n");
  }
  
  // AB is COMPLEX or COMPLEX*16 array, dimension (LDAB,N)
  T *abbuff, *abrefbuff;
  allocate_init_buffer(abbuff, abrefbuff, ldab * n);
  
  /* LDBB is INTEGER
          The leading dimension of the array BB.  LDBB >= KB+1.*/
  integer ldbb = eig_paramslist[ip].ldbb;
  if (ldbb < (kb + 1)) {
    PRINTF("ldbb < (kb + 1) but should be: LDBB >= (KB + 1). Please correct" \
           " the input data.\n");
  }
  
  // BB is COMPLEX or COMPLEX*16 array, dimension (LDBB,N)
  T *bbbuff, *bbrefbuff;
  allocate_init_buffer(bbbuff, bbrefbuff, ldbb * n);
  
  // W is REAL or DOUBLE PRECISION array, dimension (N)
  Ta *wbuff, *wrefbuff;
  allocate_init_buffer(wbuff, wrefbuff, n, 0);
  
  /* LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= N.*/
  integer ldz = eig_paramslist[ip].ldz_hbgv;
  if (ldz < 1) {
    PRINTF("ldz < 1 but should be: LDZ >= 1. Please correct the input data.\n");
  }
  if ((jobz == 'V') && (ldz < n)) {
    PRINTF("If jobz = V and ldz < n but should be: LDZ >= N. Please correct" \
           " the input data.\n");
  }
  
  // Z is COMPLEX or COMPLEX*16 array, dimension (LDZ, N)
  T *zbuff, *zrefbuff; // output buffer
  allocate_init_buffer(zbuff, zrefbuff, ldz * n, 0);
  
  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = eig_paramslist[ip].lwork_hbevd;
  integer lwork_size = 0;
  
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with " \
           "this size. Please change the input data.\n");
  }
  
  /* LRWORK is INTEGER.
     The length of the array RWORK.*/
  integer lrwork = eig_paramslist[ip].lrwork_hbevd;
  integer lrwork_size = 0;
  if ((lrwork < -1) || (lrwork == 0)) {
    PRINTF("lrwork is 0 or less than -1 and array cannot be allocated with " \
           "this size. Please change the input data.\n");
  }
  
  /* LIWORK is INTEGER.
     The length of the array IWORK.*/
  integer liwork = eig_paramslist[ip].liwork_hbevd;
  integer liwork_size = 0;
  if ((liwork < -1) || (liwork == 0)) {
    PRINTF("liwork is 0 or less than -1 and array cannot be allocated with " \
           "this size. Please change the input data.\n");
  }
  
  integer info_cpp = -1;
  
  /* If lwork/lrwork/liwork is -1, then call this API with default work
     variables. In return, work variables will be updated with array sizes
     needed.*/
  if ((lwork == -1) || (lrwork == -1) || (liwork == -1)) {
    PRINTF("lwork/lrwork/liwork is -1, so call hbgvd() to get the array " \
           "sizes\n");
    T worksize = {0};
    Ta rworksize = 0.0;
    integer iworksize = 0;
    
    // Call CPP function
    libflame::hbgvd<T, Ta>(&jobz, &uplo, &n, &ka, &kb, abbuff, &ldab, bbbuff,
              &ldbb, wbuff, zbuff, &ldz, &worksize, &lwork, &rworksize,
              &lrwork, &iworksize, &liwork, &info_cpp);
    PRINTF("info_cpp: %d, worksize: %f, rworksize: %f, iworksize: %d\n", 
            info_cpp, worksize.real, rworksize, iworksize);
    if (info_cpp == 0) {
      if (lwork == -1) {
        lwork_size = worksize.real;
      }
      if (lrwork == -1) {
        lrwork_size = rworksize;
      }
      if (liwork == -1) {
        liwork_size = iworksize;
      }
    }
  }
  
  // Check if lwork = -1 for intializing random buffer size.
  if (lwork != -1) {
    if (n <= 1) {
      if (lwork < 1) {
        PRINTF("lwork must be atleast 1, Please correct the input data.\n");
      }
    } else {
      if (jobz == 'N') {
        if (lwork < n) {
          PRINTF("lwork must be atleast n, Please correct the input data.\n");
        }
      } else if (jobz == 'V') {
        if (lwork < 2*n*2) {
          PRINTF("lwork must be atleast 2*n*2. Please correct the input" \
                 " data.\n");
        }
      }
    }
    lwork_size = lwork;
  }

  // WORK is COMPLEX or COMPLEX*16  array, dimension (MAX(1,LWORK))
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, fla_max(1, lwork_size), 0);
  
  // Check if lrwork = -1 for intializing random buffer size.
  if (lrwork != -1) {
    if (n <= 1) {
      if (lrwork < 1) {
        PRINTF("lrwork must be atleast 1, Please correct the input data.\n");
      }
    } else {
      if (jobz == 'N') {
        if (lrwork < n) {
          PRINTF("lrwork must be atleast n, Please correct the input data.\n");
        }
      } else if (jobz == 'V') {
        if (lrwork < (1 + 5*n + 2*n*2)) {
          PRINTF("lrwork must be atleast 1 + 5*n + 2*n*2, Please correct the" \
                  " input data.\n");
        }
      }
    }
    lrwork_size = lrwork;
  }
  
  // RWORK is REAL or DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
  Ta *rworkbuff = NULL, *rworkrefbuff = NULL;
  allocate_init_buffer(rworkbuff, rworkrefbuff, fla_max(1, lrwork_size), 0);
  
  // Check if liwork = -1 for intializing random buffer size.
  if (liwork != -1) {
    if ((jobz == 'N') && (n <= 1)) {
      if (liwork < 1) {
        PRINTF("liwork must be atleast 1, Please correct the input data.\n");
      }
    } else if ((jobz == 'V') && (n > 1)) {
      if (liwork < (3+5*n)) {
        PRINTF("liwork must be atleast 3 + 5*n, Please correct the" \
                " input data.\n");
      }
    }
    liwork_size = liwork;
  }
  
  // IWORK is INTEGER array, dimension (MAX(1,LIWORK))
  integer *iworkbuff = NULL, *iworkrefbuff = NULL;
  allocate_init_buffer(iworkbuff, iworkrefbuff, fla_max(1, liwork_size), 0);
  
    // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("jobz = %c\n", jobz);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("ka = %d\n", ka);
    PRINTF("kb = %d\n", kb);
    PRINTF("ldab = %d\n", ldab);
    PRINTF("Size of AB array (ldab*n) = %d\n", ldab * n);
    PRINTF("ldbb = %d\n", ldbb);
    PRINTF("Size of BB array (ldbb*n) = %d\n", ldbb * n);
    PRINTF("Size of W array (n) = %d\n", n);
    PRINTF("ldz = %d\n", ldz);
    PRINTF("Size of Z array (ldz*n) = %d\n", ldz * n);
    PRINTF("Size of WORK array (MAX(1, LWORK)) = %d\n", fla_max(1, lwork_size));
    PRINTF("LWORK = %d\n", lwork_size);
    PRINTF("Size of RWORK array (MAX(1, LRWORK)) = %d\n", fla_max(1, lrwork_size));
    PRINTF("LRWORK = %d\n", lrwork_size);
    PRINTF("Size of IWORK array (MAX(1, LIWORK)) = %d\n", fla_max(1, liwork_size));
    PRINTF("LIWORK = %d\n", liwork_size);
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
    
    // Prints BB array contents
    strncpy(arrayname, "BB input", arraysize);
    print_array<T>(arrayname, bbbuff, ldbb * n);
    strncpy(arrayname, "BB ref input", arraysize);
    print_array<T>(arrayname, bbrefbuff, ldbb * n);
    
    // Prints W array contents
    strncpy(arrayname, "W input", arraysize);
    print_array<Ta>(arrayname, wbuff, n);
    strncpy(arrayname, "W ref input", arraysize);
    print_array<Ta>(arrayname, wrefbuff, n);
    
    // Prints Z array contents
    strncpy(arrayname, "Z input", arraysize);
    print_array<T>(arrayname, zbuff, ldz * n);
    strncpy(arrayname, "Z ref input", arraysize);
    print_array<T>(arrayname, zrefbuff, ldz * n);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, fla_max(1, lwork_size));
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, fla_max(1, lwork_size));
    
    // Prints RWORK array contents
    strncpy(arrayname, "RWORK input", arraysize);
    print_array<Ta>(arrayname, rworkbuff, fla_max(1, lrwork_size));
    strncpy(arrayname, "RWORK ref input", arraysize);
    print_array<Ta>(arrayname, rworkrefbuff, fla_max(1, lrwork_size));
    
    // Prints IWORK array contents
    strncpy(arrayname, "IWORK input", arraysize);
    print_array<integer>(arrayname, iworkbuff, fla_max(1, liwork_size));
    strncpy(arrayname, "IWORK ref input", arraysize);
    print_array<integer>(arrayname, iworkrefbuff, fla_max(1, liwork_size));
  #endif

  // Call CPP function
  libflame::hbgvd<T, Ta>(&jobz, &uplo, &n, &ka, &kb, abbuff, &ldab, bbbuff,
              &ldbb, wbuff, zbuff, &ldz, workbuff, &lwork_size, rworkbuff,
              &lrwork_size, iworkbuff, &liwork_size, &info_cpp);

  // Call C function
  /* Check the typename T passed to this function template and call respective
     function.*/
  if (typeid(T) == typeid(scomplex)) {
    hbgvd_ref = (fptr_NL_LAPACK_hbgvd)dlsym(lapackModule, "chbgvd_");
  } else if (typeid(T) == typeid(dcomplex)) {
    hbgvd_ref = (fptr_NL_LAPACK_hbgvd)dlsym(lapackModule, "zhbgvd_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (hbgvd_ref == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
    closelibs();
    exit (-1);
  }
  
  integer info_ref = -1;
  hbgvd_ref(&jobz, &uplo, &n, &ka, &kb, abrefbuff, &ldab, bbrefbuff, &ldbb,
    wrefbuff, zrefbuff, &ldz, workrefbuff, &lwork_size, rworkrefbuff,
    &lrwork_size, iworkrefbuff, &liwork_size, &info_ref);
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
      
      // Prints BB array contents
      strncpy(arrayname, "BB output", arraysize);
      print_array<T>(arrayname, bbbuff, ldbb * n);
      strncpy(arrayname, "BB ref output", arraysize);
      print_array<T>(arrayname, bbrefbuff, ldbb * n);
      
      // Prints W array contents
      strncpy(arrayname, "W output", arraysize);
      print_array<Ta>(arrayname, wbuff, n);
      strncpy(arrayname, "W ref output", arraysize);
      print_array<Ta>(arrayname, wrefbuff, n);
      
      // Prints Z array contents
      strncpy(arrayname, "Z output", arraysize);
      print_array<T>(arrayname, zbuff, ldz * n);
      strncpy(arrayname, "Z ref output", arraysize);
      print_array<T>(arrayname, zrefbuff, ldz * n);
    
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, fla_max(1, lwork_size));
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, fla_max(1, lwork_size));
      
      // Prints RWORK array contents
      strncpy(arrayname, "RWORK output", arraysize);
      print_array<Ta>(arrayname, rworkbuff, fla_max(1, lrwork_size));
      strncpy(arrayname, "RWORK ref output", arraysize);
      print_array<Ta>(arrayname, rworkrefbuff, fla_max(1, lrwork_size));
      
      // Prints IWORK array contents
      strncpy(arrayname, "IWORK output", arraysize);
      print_array<integer>(arrayname, iworkbuff, fla_max(1, liwork_size));
      strncpy(arrayname, "IWORK ref output", arraysize);
      print_array<integer>(arrayname, iworkrefbuff, fla_max(1, liwork_size));
    #endif
    
    double diff = computeError<T>(ldab, n, abrefbuff, abbuff);
    diff += computeError<T>(ldbb, n, bbrefbuff, bbbuff);
    if (jobz == 'V') {
      diff += computeError<T>(ldz, n, zbuff, zrefbuff);
    }
    diff += computeError<Ta>(1, n, wbuff, wrefbuff);
    diff += computeError<T>(1, fla_max(1, lwork_size), workbuff, workrefbuff);
    diff += computeError<Ta>(1, fla_max(1, lrwork_size), rworkbuff, rworkrefbuff);
    diff += computeError<integer>(1, fla_max(1, liwork_size), iworkbuff,
                                  iworkrefbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, abs(diff), SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  if ((info_cpp > 0) || (info_ref > 0)) {
    if ((info_cpp <= n) || (info_ref <= n)) {
      PRINTF("the algorithm failed to converge: i off-diagonal elements of" \
             " an intermediate tridiagonal form did not converge to zero\n");
    }
    else if ((info_cpp > n) || (info_ref > n)) {
      PRINTF("if INFO = N + i, for 1 <= i <= N, then CPBSTF" \
             " returned INFO = i: B is not positive definite." \
             " The factorization of B could not be completed and" \
             " no eigenvalues or eigenvectors were computed.\n");
    }
  }
  
  // Free up the buffers.
  delete[] abbuff; delete[] abrefbuff;
  delete[] bbbuff; delete[] bbrefbuff;
  delete[] wbuff; delete[] wrefbuff;
  delete[] zbuff; delete[] zrefbuff;
  delete[] workbuff; delete[] workrefbuff;
  delete[] rworkbuff; delete[] rworkrefbuff;
  delete[] iworkbuff; delete[] iworkrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex and float as typenames.*/
TEST(LAPACKCPP_hbgvd, CHBGVD) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hbgvd_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex and double as typenames.*/
TEST(LAPACKCPP_hbgvd, ZHBGVD) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
	  PRINTF("index: %d\n", index);
    hbgvd_test<dcomplex, double> (index);
  }
}