/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file liblame_geqrf.cc
 *  @brief Test application to validate geqrf() using CPP template interface.
 *  */

#include <gtest/gtest.h>
#include "main.h"
#include "libflame_test.hh"

#if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1))
// Array to store array name to print.
char arrayname[20] = "";
integer arraysize = sizeof(arrayname);
#endif

template<typename T, typename Ta>
void validate_geqrf(integer m_A,
  integer n_A,
  T *A,
  T *A_test,
  T *T_test,
  double* residual)
{
  T *Q = NULL, *R = NULL, *Ibuff = NULL, *work = NULL;
  integer cs_A, min_A;
  integer lwork = -1, tinfo;

  PRINTF("%s() entry\n", __FUNCTION__);
  cs_A = m_A;
  min_A = fla_min(m_A, n_A);

  // Create Q and R matrices.
  create_matrix<T>(&Q, m_A, m_A);
  create_matrix<T>(&R, m_A, n_A);
  reset_matrix<T>(m_A, m_A, Q, m_A);
  reset_matrix<T>(m_A, n_A, R, cs_A);

  // Create Identity matrix to validate orthogonal property of matrix Q
  create_matrix<T>(&Ibuff, m_A, m_A);

  /* Extract R matrix and elementary reflectors from the input/output
     matrix parameter A_test.*/
  if (m_A <= n_A) {
    copy_matrix<T>("full", m_A, m_A, A_test, m_A, Q, m_A);
    copy_matrix<T>("Upper", m_A, n_A, A_test, m_A, R, m_A);
  }  else {
    copy_matrix<T>("full", m_A, n_A, get_m_ptr<T>(A_test, 1, 0, m_A),
                    m_A, get_m_ptr<T>(Q, 1, 0, m_A), m_A);
    copy_matrix<T>("Upper", n_A, n_A, A_test, m_A, R, m_A);
  }
  
  T twork;
  Ta norm, norm_A, eps, resid1, resid2;

  if ((typeid(T) == typeid(float)) || (typeid(T) == typeid(double))) {
    /* orgqr api generates the Q martrix using the elementary reflectors and scalar
       factor values*/
    libflame::orgqr<T>(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork,
                      &tinfo);
    PRINTF("orgqr() tinfo: %d, lwork:%d\n", tinfo, get_work_value<T>(&twork));
  } else if ((typeid(T) == typeid(scomplex)) || (typeid(T) == typeid(dcomplex))) {
    /* ungqr api generates the Q martrix using the elementary reflectors and scalar
       factor values*/
    libflame::ungqr<T>(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork,
                      &tinfo);
    PRINTF("ungqr() tinfo: %d, lwork:%d\n", tinfo, get_work_value<T>(&twork));
  }

  lwork = get_work_value<T>(&twork);
  create_vector<T>(&work, lwork);
  
  if ((typeid(T) == typeid(float)) || (typeid(T) == typeid(double))) {
    libflame::orgqr<T>(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &tinfo);
    PRINTF("orgqr() tinfo: %d\n", tinfo);
  } else if ((typeid(T) == typeid(scomplex)) || (typeid(T) == typeid(dcomplex))) {
    libflame::ungqr<T>(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &tinfo);
    PRINTF("ungqr() tinfo: %d\n", tinfo);
  }
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    // Prints Q array contents
    strncpy(arrayname, "Q output", arraysize);
    print_array_data(arrayname, Q, m_A * m_A);
  #endif
  alpha_beta_init<T>();
  char transa = 'T';
  char transb = 'N';
  if ((typeid(T) == typeid(scomplex)) || (typeid(T) == typeid(dcomplex))) {
    transa = 'C';
  }
  
  libflame_utils::gemm<T>(&transa, &transb, &m_A, &n_A, &m_A, (T *)n_one, Q,
                    &m_A, A, &m_A, (T *)one, R, &m_A);
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    // Prints R array contents
    strncpy(arrayname, "R output", arraysize);
    print_array<T>(arrayname, R, m_A * n_A);
  #endif
  
  norm_A = libflame::lange<T, Ta>("1", &m_A, &n_A, A, &m_A, NULL);
  norm = libflame::lange<T, Ta>("1", &m_A, &n_A, R, &m_A, NULL);
  
  PRINTF("norm:%lf norm_A: %lf\n", norm, norm_A);

  eps = libflame::lamch<Ta>("P");
  PRINTF("eps: %lf\n", eps);
  resid1 = norm/(eps * norm_A * n_A);
  PRINTF("resid1: %lf\n", resid1);

  /* Test 2
     compute norm(I - Q*Q') / (N * EPS)*/
  libflame::laset<T>("full", &m_A, &m_A, (T *)zero, (T *)one, Ibuff, &m_A);
  
  libflame_utils::gemm<T>(&transb, &transa, &m_A, &m_A, &m_A, (T *)n_one,
                    Q, &m_A, Q, &m_A, (T *)one, Ibuff, &m_A);
  norm = libflame::lange<T, Ta>("1", &m_A, &m_A, Ibuff, &m_A, NULL);
  PRINTF("norm: %lf\n", norm);
  resid2 = norm/(eps * n_A);
  PRINTF("resid2: %lf\n", resid2);

  *residual = fla_max(resid1, resid2);

  // Free up buffers
  free_matrix(R);
  free_matrix(Q);
  free_matrix(Ibuff);
  free_vector(work);
  PRINTF("%s() exit\n", __FUNCTION__);
}

/*! @brief  geqrf_test is function template for geqrf() functions.
      T can be float, double
 * @details
 * \b Purpose:
    \verbatim
    geqrf_test is function template for geqrf() functions.
    T can be float, double
    
    geqrf_test() function template calls C and CPP based library APIs with
    valid test values, calculate the differences in output if INFO == 0.
    And passses the test case if difference <= threshold.
    Fails the test case if difference > threshold or INFO < 0.

    Real Reference:
    http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga7cb54fa1727bf0166523036f4948bc56.html#ga7cb54fa1727bf0166523036f4948bc56
    Double Reference:
    https://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga3766ea903391b5cf9008132f7440ec7b.html#ga3766ea903391b5cf9008132f7440ec7b
    \endverbatim
  
 * @params[in] IP
          IP is INTEGER
          Index of array in config file.

 * @return VOID
           Nothing.
 * */
template<typename T>
void geqrf_test(int ip)
{
  // Initialise random number generators with timestamp
  srand (SRAND_SEED_VALUE);
  
  /* M is INTEGER
          The number of rows of the matrix A. M >= 0.*/
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
  if (lda < fla_max(1, m)) {
    PRINTF("lda < fla_max(1, m) but it should be: LDA >= fla_max(1,M). Please " \
           "correct the input data.\n");
  }
  
  /* A is REAL or DOUBLE PRECISION or COMPLEX or COMPLEX*16 array,
      dimension (LDA,N)*/
  T *abuff, *arefbuff;
  create_matrix<T>(&abuff, m, n);
  create_matrix<T>(&arefbuff, m, n);
  rand_matrix<T>(abuff, m, n, lda);
  copy_matrix<T>("full", m, n, abuff, lda, arefbuff, lda);
  
  /* TAU is REAL or DOUBLE PRECISION or COMPLEX or COMPLEX*16 array,
      dimension (fla_min(M,N))*/
  T *taubuff;
  integer tausize = fla_min(m, n);
  create_vector<T>(&taubuff, tausize);
  
  /* LWORK is INTEGER
          The dimension of the array WORK.  LWORK >= fla_max(1,N).
          For optimum performance LWORK >= N*NB, where NB is
          the optimal blocksize.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.*/
  integer lwork = eig_paramslist[ip].lwork_geqrf;
  integer lwork_size = lwork;
  if ((lwork < -1) || (lwork == 0)) {
    PRINTF("lwork is 0 or less than -1 and array cannot be allocated with" \
           " this size. Please change the input data.\n");
  }
  integer info_cpp = -1;
  
  /* If lwork is -1, then call this API with default work variables.
     In return, work variable will be updated with array sizes needed.*/
  if (lwork == -1) {
    PRINTF("lwork is -1, so call geqrf() to get the array sizes.\n");
    T worksize = {0};
    libflame::geqrf<T>(&m, &n, abuff, &lda, taubuff, &worksize, &lwork_size,
                &info_cpp);
    PRINTF("info_cpp: %d, lwork: %d\n", info_cpp, get_work_value<T>(&worksize));
    if (info_cpp == 0) {
      lwork_size = get_work_value<T>(&worksize);
    }
  }
  
  /* WORK is REAL or DOUBLE PRECISION or COMPLEX or COMPLEX*16 array,
      dimension (MAX(1,LWORK))*/
  T *workbuff = NULL;
  create_vector<T>(&workbuff, lwork_size);
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("m = %d\n", m);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of TAU array (fla_min(m,n)) = %d\n", fla_min(m, n));
    PRINTF("lwork = %d\n", lwork_size);
    PRINTF("Size of WORK array (MAX(1,LWORK)) = %d\n", fla_max(1, lwork_size));
  #endif
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_INPUT_ARRAYS) && (PRINT_INPUT_ARRAYS == 1))
    // Print all input arrays if PRINT_INPUT_ARRAYS macro is enabled
    PRINTF("Printing all Input arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A input", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref input", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
    
    // Prints TAU array contents
    strncpy(arrayname, "TAU input", arraysize);
    print_array<T>(arrayname, taubuff, fla_min(m, n));
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK input", arraysize);
    print_array<T>(arrayname, workbuff, fla_max(1, lwork));
  #endif

  // Call C++ function.
  info_cpp = -1;
  libflame::geqrf<T>(&m, &n, abuff, &lda, taubuff, workbuff, &lwork_size,
                &info_cpp);
  PRINTF ("geqrf() info_cpp: %d\n", info_cpp);
  
  // Calculate the differences of buffers.
  if (info_cpp == 0) {
    double residual = 0.0;
    #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all output arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A output", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref output", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
    
    // Prints TAU array contents
    strncpy(arrayname, "TAU output", arraysize);
    print_array<T>(arrayname, taubuff, fla_min(m, n));
    
    // Prints WORK array contents
    strncpy(arrayname, "WORK output", arraysize);
    print_array<T>(arrayname, workbuff, fla_max(1, lwork));
  #endif
    if ((typeid(T) == typeid(float)) || (typeid(T) == typeid(scomplex))) {
      validate_geqrf<T, float>(m, n, arefbuff, abuff, taubuff, &residual);
    } else if ((typeid(T) == typeid(double)) || (typeid(T) == typeid(dcomplex))) {
      validate_geqrf<T, double>(m, n, arefbuff, abuff, taubuff, &residual);
    }
    PRINTF("residual: %lf\n", residual);
    EXPECT_NEAR(0.0, residual, SYM_EIGEN_THRESHOLD);
  } else {
    PRINTF("Info returned by CPP API is not successful to validate" \
           " the changes.\n");
    EXPECT_FALSE(info_cpp < 0);
  }
  
  // Free up the buffers
  delete[] abuff;
  delete[] arefbuff;
  delete[] taubuff;
  delete[] workbuff;
}

/* Use TEST macro and call C++ test function template with
   float as typenames.*/
TEST(LAPACKCPP_geqrf, SGEQRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    geqrf_test<float> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   double as typenames.*/
TEST(LAPACKCPP_geqrf, DGEQRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    geqrf_test<double> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   scomplex as typenames.*/
TEST(LAPACKCPP_geqrf, CGEQRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    geqrf_test<scomplex> (index);
  }
}

/* Use TEST macro and call C++ test function template with
   dcomplex as typenames.*/
TEST(LAPACKCPP_geqrf, ZGEQRF) {
  for (short int index = 0; index < NUM_SUB_TESTS; index++) {
    PRINTF("index: %d\n", index);
    geqrf_test<dcomplex> (index);
  }
}