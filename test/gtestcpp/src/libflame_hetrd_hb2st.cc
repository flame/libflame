/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_hetrd_hb2st.cc
 *  @brief Test application to validate hetrd_hb2st() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

/*! @brief  hetrd_hb2st_test is function template for hetrd_hb2st() functions.
			T can be scomplex, dcomplex
      Ta can be float, double
 * @details
 * \b Purpose:
    \verbatim
	  hetrd_hb2st_test is function template for hetrd_hb2st() functions.
	  T can be scomplex, dcomplex
    Ta can be float, double
	  
	  hetrd_hb2st_test() function template calls C and CPP based lbrary APIs with
	  valid test values, calculate the differences in output if INFO = 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

 	  Complex reference:
    http://www.netlib.org/lapack/explore-html/d3/db9/group__complex_o_t_h_e_rcomputational_ga040adda5091e3515acc9d869c3179418.html#ga040adda5091e3515acc9d869c3179418
    Complex double reference:
    http://www.netlib.org/lapack/explore-html/d0/da6/group__complex16_o_t_h_e_rcomputational_gacb962c8889c699b70126db722d2c05d0.html#gacb962c8889c699b70126db722d2c05d0
    \endverbatim
	
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template< typename T,  typename Ta >
void hetrd_hb2st_test(int ip)
{
  typedef int (*Fptr_NL_LAPACK_hetrd_hb2st)(char *stage1, char *vect,
                    char *uplo, integer *n, integer *kd, T *ab,
                    integer *ldab, Ta *d, Ta * e, T *hous,
                    integer *lhous, T *work, integer *lwork,
                    integer *info);
  Fptr_NL_LAPACK_hetrd_hb2st HETRD_HB2ST;
  
  
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* STAGE1 is CHARACTER*1
          = 'N':  "No": to mention that the stage 1 of the reduction  
                  from dense to band using the chetrd_he2hb routine
                  was not called before this routine to reproduce AB. 
                  In other term this routine is called as standalone. 
          = 'Y':  "Yes": to mention that the stage 1 of the 
                  reduction from dense to band using the chetrd_he2hb 
                  routine has been called to produce AB (e.g., AB is
                  the output of chetrd_he2hb.*/
  char stage1 = eig_paramslist[ip].stage1;
  if ((stage1 != 'N') && (stage1 != 'Y')) {
    PRINTF("stage1 should be N or Y. Please correct the input data.\n");
  }
  
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
  
  /* KD is INTEGER
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.*/
  integer kd = eig_paramslist[ip].kd;
  if (kd < 0) {
    PRINTF("kd < 0 but should be: kd >= 0. Please correct the input data.");
  }
  
  // LDAB is INTEGER. The leading dimension of the array AB.  LDAB >= KD+1.
  integer ldab = eig_paramslist[ip].ldab;
  if (ldab < (kd + 1)) {
    PRINTF("ldab < (kd+1) but it should be: LDAB >= (kd+1). Please " \
           "correct the input data.\n");
  }
  
  // AB is COMPLEX or COMPLEX*16 array, dimension (LDAB, N)
  T *abbuff, *abrefbuff;
  allocate_init_buffer(abbuff, abrefbuff, ldab * n);
  
  // D is REAL or DOUBLE PRECISION array, dimension (N)
  Ta *dbuff, *drefbuff;
  allocate_init_buffer(dbuff, drefbuff, n, 0);
  
  // E is REAL or DOUBLE PRECISION array, dimension (N-1)
  Ta *ebuff, *erefbuff;
  allocate_init_buffer(ebuff, erefbuff, n-1, 0);

  /* LHOUS is INTEGER
          The dimension of the array HOUS. LHOUS = MAX(1, dimension)
          If LWORK = -1, or LHOUS=-1,
          then a query is assumed; the routine
          only calculates the optimal size of the HOUS array, returns
          this value as the first entry of the HOUS array, and no error
          message related to LHOUS is issued by XERBLA.
          LHOUS = MAX(1, dimension) where
          dimension = 4*N if VECT='N'
          not available now if VECT='H' */
  integer lhous = eig_paramslist[ip].lhous;
  integer lhous_size = lhous;
  if ((lhous < -1) || (lhous == 0)) {
    PRINTF("lhous is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  
  /* LWORK is INTEGER.
     The length of the array WORK.*/
  integer lwork = eig_paramslist[ip].lwork_hetrd_hb2st;
  integer lwork_size = lwork;
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  
  integer info_cpp = -1;
  
  /* If lwork/lhous is -1, then call this API with default work
     variables. In return, work variables will be updated with array sizes
     needed.*/
  if ((lwork == -1) || (lhous == -1)) {
    PRINTF("lwork/lhous is -1, so call hetrd_2stage() to get the array" \
          " sizes.\n");
    T lhoussize = {0};
    T worksize = {0};
    
    // Call CPP function
    libflame::hetrd_hb2st<T, Ta>(&stage1, &vect, &uplo, &n, &kd, abbuff, &ldab,
                        dbuff, ebuff, &lhoussize, &lhous_size, &worksize,
                        &lwork_size, &info_cpp);
    PRINTF("info_cpp: %d, lhoussize: %f, worksize: %f\n", info_cpp, \
          lhoussize.real, worksize.real);
    if (info_cpp == 0) {
      if (lwork == -1) {
        lwork_size = worksize.real;
      }
      if (lhous == -1) {
        lhous_size = lhoussize.real;
      }
    }
  }
  
  // HOUS is COMPLEX or COMPLEX*16  array, dimension (LHOUS))
  T *housbuff = NULL, *housrefbuff = NULL;
  allocate_init_buffer(housbuff, housrefbuff, lhous_size, 0);
  
  // WORK is COMPLEX or COMPLEX*16  array, dimension (LWORK)
  T *workbuff = NULL, *workrefbuff = NULL;
  allocate_init_buffer(workbuff, workrefbuff, lwork_size, 0);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("stage1 = %c\n", stage1);
    PRINTF("vect = %c\n", vect);
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("kd = %d\n", kd);
    PRINTF("ldab = %d\n", ldab);
    PRINTF("Size of AB array (ldab*n) = %d\n", ldab * n);
    PRINTF("Size of D array (n) = %d\n", n);
    PRINTF("Size of E array (n-1) = %d\n", n-1);
    PRINTF("lhous = %d\n", lhous_size);
    PRINTF("Size of HOUS array (lhous) = %d\n", lhous_size);
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
    print_array<Ta>(arrayname, ebuff, n-1);
    strncpy(arrayname, "E ref input", arraysize);
    print_array<Ta>(arrayname, erefbuff, n-1);
    
    // Prints HOUS array contents
    strncpy(arrayname, "HOUS input", arraysize);
    print_array<T>(arrayname, housbuff, lhous_size);
    strncpy(arrayname, "HOUS ref input", arraysize);
    print_array<T>(arrayname, housrefbuff, lhous_size);
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, lwork_size);
    strncpy(arrayname, "WORK ref input", arraysize);
    print_array<T>(arrayname, workrefbuff, lwork_size);
  #endif
  
  // Call CPP function
  info_cpp = -1;
  libflame::hetrd_hb2st<T, Ta>(&stage1, &vect, &uplo, &n, &kd, abbuff, &ldab,
                        dbuff, ebuff, housbuff, &lhous_size, workbuff,
                        &lwork_size, &info_cpp);

  // Call C function
  if (typeid(T) == typeid(scomplex)) {
    HETRD_HB2ST = (Fptr_NL_LAPACK_hetrd_hb2st)dlsym(lapackModule, \
                          "chetrd_hb2st_");
  } else if (typeid(T) == typeid(dcomplex)) {
    HETRD_HB2ST = (Fptr_NL_LAPACK_hetrd_hb2st)dlsym(lapackModule, \
                          "zhetrd_hb2st_");
  } else {
	  PRINTF("Invalid typename is passed to %s() function template.\n",
           __FUNCTION__);
  }
  
  if (HETRD_HB2ST == NULL) {
    PRINTF("Could not get the symbol. Exiting...\n");
	  closelibs();
    exit(-1);
  }
  
  integer info_ref = -1;
  HETRD_HB2ST(&stage1, &vect, &uplo, &n, &kd, abrefbuff, &ldab, drefbuff,
              erefbuff, housrefbuff, &lhous_size, workrefbuff, &lwork_size,
              &info_ref);
  
  // Calculate the differences of buffers.
  if ((info_cpp == 0) && (info_ref == 0)) {
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
      print_array<Ta>(arrayname, ebuff, n-1);
      strncpy(arrayname, "E ref output", arraysize);
      print_array<Ta>(arrayname, erefbuff, n-1);
      
      // Prints HOUS array contents
      strncpy(arrayname, "HOUS output", arraysize);
      print_array<T>(arrayname, housbuff, lhous_size);
      strncpy(arrayname, "HOUS ref output", arraysize);
      print_array<T>(arrayname, housrefbuff, lhous_size);
    
      // Prints WORK array contents
      strncpy(arrayname, "WORK output", arraysize);
      print_array<T>(arrayname, workbuff, lwork_size);
      strncpy(arrayname, "WORK ref output", arraysize);
      print_array<T>(arrayname, workrefbuff, lwork_size);
    #endif
  
    double diff = computeError<T>(ldab, n, abbuff, abrefbuff);
    diff += computeError<Ta>(1, n, dbuff, drefbuff);
    diff += computeError<Ta>(1, n-1, ebuff, erefbuff);
    diff += computeError<T>(1, lhous_size, housbuff, housrefbuff);
    diff += computeError<T>(1, lwork_size, workrefbuff, workbuff);
    PRINTF("diff: %lf\n", diff);
    EXPECT_NEAR(0.0, diff, SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP or C API is not successful to compare" \
            " differences.\n");
    EXPECT_FALSE((info_cpp < 0) || (info_ref < 0));
  }
  
  // Free up the buffers
  delete[] abbuff; delete[] abrefbuff;
  delete[] dbuff; delete[] drefbuff;
  delete[] ebuff; delete[] erefbuff;
  delete[] housbuff; delete[] housrefbuff;
  delete[] workbuff; delete[] workrefbuff;
}

/* Use TEST macro and call C++ test function template with
   scomplex, float as typenames.*/
TEST(LAPACKCPP_hetrd_hb2st, CHETRD_HB2ST) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetrd_hb2st_test<scomplex, float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex, double as typenames.*/
TEST(LAPACKCPP_hetrd_hb2st, ZHETRD_HB2ST) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    hetrd_hb2st_test<dcomplex, double> (index);
  }
}