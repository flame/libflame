/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file main.h
 *  @brief Defines structures, macros, function declarations to use in APIs
           of CPP template interface.
 *  */

#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <dlfcn.h>

extern "C" {
#include <FLA_f2c.h>
}

#define NUM_SUB_TESTS (2)

/* Global variables declaration */

#define NETLIB_LAPACK_LIB_PATH "../../lapack/lib/liblapack.so"
#define NETLIB_BLAS_LIB_PATH "../../blis/lib/libblis-mt.so"

extern char NETLIB_LAPACK_LIB[60];
extern char NETLIB_BLAS_LIB[60];

extern void *blasModule, *lapackModule;

/* structure to hold eigen parameters */
typedef struct EIG_paramlist_t {
  char jobz;  /* JOBZ is CHARACTER*1
                  = 'N':  Compute eigenvalues only;
                  = 'V':  Compute eigenvalues and eigenvectors.*/
  char uplo;  /*UPLO is CHARACTER*1
                  = 'U':  Upper triangle of A is stored;
                  = 'L':  Lower triangle of A is stored.*/
  char job;         // Must be 'N', 'P', 'S' or 'B'
  integer n;  // N is INTEGER. The order of the matrix A.  N >= 0.
  integer subda;  // The number of subdiagonals of matrix A- hbev
  integer sda;  // The number of superdiagonals of the matrix A
  integer ldab;  /* LDAB is INTEGER
                    The leading dimension of the array AB.  LDAB >= KD + 1.*/
  integer ldz;  /* LDZ is INTEGER
                    The leading dimension of the array Z.  LDZ >= 1, and if
                    JOBZ = 'V', LDZ >= max(1,N).*/
  // Added for hbev_2stage
  char jobz_2stage;  /* JOBZ is CHARACTER*1 ('V' is not supported)
                  = 'N':  Compute eigenvalues only;
                  = 'V':  Compute eigenvalues and eigenvectors.
                          Not available in this release.*/
  integer lwork;  /* LWORK is INTEGER
                      The length of the array WORK. LWORK >= 1, when N <= 1;
                      otherwise  
                      If JOBZ = 'N' and N > 1, LWORK must be queried.
                                               LWORK = MAX(1, dimension) where
                                               dimension = (2KD+1)*N + KD*NTHREADS
                                               where KD is the size of the band.
                                               NTHREADS is the number of threads used when
                                               openMP compilation is enabled, otherwise =1.
                      If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available.

                      If LWORK = -1, then a workspace query is assumed; the routine
                      only calculates the optimal sizes of the WORK, RWORK and
                      IWORK arrays, returns these values as the first entries of
                      the WORK, RWORK and IWORK arrays, and no error message
                      related to LWORK or LRWORK or LIWORK is issued by XERBLA.*/
  // Added for hbevd
  integer lwork_hbevd;  /* LWORK is INTEGER
                            The length of the array WORK. LWORK >= 1, when N <= 1;
                            otherwise  
                            If JOBZ = 'N' and N > 1, LWORK must be queried.
                                                     LWORK = MAX(1, dimension) where
                                                     dimension = (2KD+1)*N + KD*NTHREADS
                                                     where KD is the size of the band.
                                                     NTHREADS is the number of threads used when
                                                     openMP compilation is enabled, otherwise =1.
                            If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available.

                            If LWORK = -1, then a workspace query is assumed; the routine
                            only calculates the optimal sizes of the WORK, RWORK and
                            IWORK arrays, returns these values as the first entries of
                            the WORK, RWORK and IWORK arrays, and no error message
                            related to LWORK or LRWORK or LIWORK is issued by XERBLA.*/
  integer lrwork_hbevd;  /* LRWORK is INTEGER
                              The dimension of array RWORK.
                              If N <= 1,               LRWORK must be at least 1.
                              If JOBZ = 'N' and N > 1, LRWORK must be at least N.
                              If JOBZ = 'V' and N > 1, LRWORK must be at least
                                            1 + 5*N + 2*N**2.

                              If LRWORK = -1, then a workspace query is assumed; the
                              routine only calculates the optimal sizes of the WORK, RWORK
                              and IWORK arrays, returns these values as the first entries
                              of the WORK, RWORK and IWORK arrays, and no error message
                              related to LWORK or LRWORK or LIWORK is issued by XERBLA.*/
  integer liwork_hbevd;  /* LIWORK is INTEGER
                              The dimension of array IWORK.
                              If JOBZ = 'N' or N <= 1, LIWORK must be at least 1.
                              If JOBZ = 'V' and N > 1, LIWORK must be at least 3 + 5*N .

                              If LIWORK = -1, then a workspace query is assumed; the
                              routine only calculates the optimal sizes of the WORK, RWORK
                              and IWORK arrays, returns these values as the first entries
                              of the WORK, RWORK and IWORK arrays, and no error message
                              related to LWORK or LRWORK or LIWORK is issued by XERBLA.*/
  // Added for hbevx
  char range;  /* RANGE is CHARACTER*1
                    = 'A': all eigenvalues will be found;
                    = 'V': all eigenvalues in the half-open interval (VL,VU]
                           will be found;
                    = 'I': the IL-th through IU-th eigenvalues will be found.*/
  integer ldq;  /* Q is COMPLEX array, dimension (LDQ, N)
                    If JOBZ = 'V', the N-by-N unitary matrix used in the
                                    reduction to tridiagonal form.
                    If JOBZ = 'N', the array Q is not referenced.*/
  double vl;  /* VL is REAL or DOUBLE PRECISION
                  If RANGE='V', the lower bound of the interval to
                  be searched for eigenvalues. VL < VU.
                  Not referenced if RANGE = 'A' or 'I'.*/
  double vu;  /* VU is REAL or DOUBLE PRECISION
                  If RANGE='V', the upper bound of the interval to
                  be searched for eigenvalues. VL < VU.
                  Not referenced if RANGE = 'A' or 'I'.*/
  integer il;  /* IL is INTEGER
                    If RANGE='I', the index of the
                    smallest eigenvalue to be returned.
                    1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
                    Not referenced if RANGE = 'A' or 'V'.*/
  integer iu; /* IU is INTEGER
                  If RANGE='I', the index of the
                  largest eigenvalue to be returned.
                  1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
                  Not referenced if RANGE = 'A' or 'V'.*/
  double abstol;  /* ABSTOL is REAL or DOUBLE PRECISION
                    The absolute error tolerance for the eigenvalues.*/
  // Added for hbgv
  integer ldz_hbgv; /* Z is COMPLEX array, dimension (LDZ, N)
                        If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
                        eigenvectors, with the i-th column of Z holding the
                        eigenvector associated with W(i). The eigenvectors are
                        normalized so that Z**H*B*Z = I.
                        If JOBZ = 'N', then Z is not referenced.*/
  integer sdb; // The number of superdiagonals of the matrix B
  integer subdb; // The number of subdiagonals of the matrix B
  integer ldab_hbgv; // LDAB is INTEGER. The leading dimension of the array AB.  LDAB >= KA+1.
  integer ldbb; // LDBB is INTEGER. The leading dimension of the array BB.  LDBB >= KB+1.
  
  // Added for hbgvd
  integer lwork_hbgvd; /* LWORK is INTEGER
                            The dimension of the array WORK.
                            If N <= 1,               LWORK >= 1.
                            If JOBZ = 'N' and N > 1, LWORK >= N.
                            If JOBZ = 'V' and N > 1, LWORK >= 2*N**2.

                            If LWORK = -1, then a workspace query is assumed; the routine
                            only calculates the optimal sizes of the WORK, RWORK and
                            IWORK arrays, returns these values as the first entries of
                            the WORK, RWORK and IWORK arrays, and no error message
                            related to LWORK or LRWORK or LIWORK is issued by XERBLA.*/
  // Added for hbgvx
  integer ldq_hbgvx; /* LDQ is INTEGER
                          The leading dimension of the array Q.  If JOBZ = 'N',
                          LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).*/
  char vect; /* VECT is CHARACTER*1
                  = 'N':  do not form Q;
                  = 'V':  form Q;
                  = 'U':  update a matrix X, by forming X*Q.*/
  char vect_hbgst; /* VECT is CHARACTER*1
                        = 'N':  do not form the transformation matrix X;
                        = 'V':  form X.*/
  integer ldx; /* LDX is INTEGER
                    The leading dimension of the array X.
                    LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise.*/
  
} EIG_paramlist;

/* structure to hold Linear solver parameters */
typedef struct Lin_solver_paramlist_t {
  // Added for hecon
  char uplo; /* UPLO is CHARACTER*1
                  Specifies whether the details of the factorization are stored
                  as an upper or lower triangular matrix.
                  = 'U':  Upper triangular, form is A = U*D*U**H;
                  = 'L':  Lower triangular, form is A = L*D*L**H.*/
  int n;         // The order of matrix A; N >= 0.
  int lda;       // The leading dimension of the array A. LDA >= max(1,N).
  double anorm; // ANORM is DOUBLE PRECISION. The 1-norm of the original matrix A.
} Lin_solver_paramlist;

extern EIG_paramlist eig_paramslist[NUM_SUB_TESTS];
extern Lin_solver_paramlist lin_solver_paramslist[NUM_SUB_TESTS];

/*! @brief  Read_EIG_params is function used to read and initialize 
      Symmetric Eigen values/vectors routines input data into global array.
 * @details
 * \b Purpose:
    \verbatim
    This function reads parameters needed for Eigen APIs 
    from the config settings file 'EIG_PARAMS.dat' and saves in the 
    'eig_paramslist' structure array.
    \endverbatim
	
 * @param[in] file_name
          file_name is charater array.
          Used to specify the name of the file to read the data.

 * @return void
          Nothing.
 * */
void Read_EIG_params(const char *file_name);

/*! @brief  Read_Lin_solver_params is function used to read and initialize 
      Linear solver routines input data into global array.
 * @details
 * \b Purpose:
    \verbatim
    This function reads parameters needed for Linear solver APIs 
   from the config settings file 'LIN_SLVR.dat' and saves in the 
   'lin_solver_paramslist' structure array.
    \endverbatim
	
 * @param[in] file_name
          file_name is charater array.
          Used to specify the name of the file to read the data.

 * @return void
          Nothing.
 * */
void Read_Lin_solver_params(const char *file_name);

/*! @brief  closelibs is function used to close the blasModule and lapackModule
            dynamic library.
 * @details
 * \b Purpose:
    \verbatim
    This function is used to close the blasModule and lapackModule
    dynamic library.
    \endverbatim
	
 * @param void
          Nothing.

 * @return void
          Nothing.
 * */
void closelibs(void);

// Macro to enable status/error messages of test APIs.
#define PRINT_MSGS 1

#if PRINT_MSGS 
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...) 
#endif

/* Macro to enable print of contents of arrays.
   0 - disables the print of array contents.
   1 - enables the print of array contents.
   Note: This macro need to be enabled along with PRINT_INPUT_ARRAYS or
   PRINT_OUTPUT_ARRAYS macros to print the array contents.*/
#define PRINT_ARRAYS 0

/* Macro to enable print of contents of input arrays.
   0 - disables the print of input array contents.
   1 - enables the print of input array contents.
   Note: This macro need to be enabled along with PRINT_ARRAYS to
   print the array contents.*/
#define PRINT_INPUT_ARRAYS 1

/* Macro to enable print of contents of output arrays.
   0 - disables the print of output array contents.
   1 - enables the print of output array contents.
   Note: This macro need to be enabled along with PRINT_ARRAYS to
   print the array contents.*/
#define PRINT_OUTPUT_ARRAYS 1

/* Macro to print input values other than array contents.
   0 - disables the print of input values.
   1 - enables the print of input values.*/
#define PRINT_INPUT_VALUES 0

/* Macro to enable print of all contents of each array.
   0 - considers arrays of ARRAY_PRINT_SIZE size for PRINT_ARRAYS.
   1 - considers maximum size of arrays allocated for PRINT_ARRAYS.*/
#define MAX_ARRAY_PRINT_SIZE 0

/* Macro to print specific size of contents of the buffer.
   Increase or decrease the size as needed.*/
#define ARRAY_PRINT_SIZE 10

// Threshold values of different APIs

// Threshold value from dsg.in for Symmetric Eigenvalues and eigenvectors APIs
#define SYM_EIGEN_THRESHOLD 20.0

// Threshold value from ctest.in or ztest.in for Linear solver APIs
#define LIN_SLVR_THRESHOLD 30.0

#endif // MAIN_H