/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

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


} EIG_paramlist;

extern EIG_paramlist eig_paramslist[NUM_SUB_TESTS];

/*! @brief  Read_EIG_params is function used to read and initialize 
      Symmetric Eigen values/vectors routines input data into global array.
 * @details
 * \b Purpose:
    \verbatim
    This function reads parameters needed for Eigen APIs 
    from the config settings file 'EIG_PARAMS.dat' and saves in the 
    'eig_paramslist' structure array
    
    \endverbatim
	
 * @param[in] file_name
    file_name is charater array.
    Used to specify the name of the file to read the data.

 * @return void
    Nothing.
 * */
void Read_EIG_params(const char *file_name);

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

// Threshold values of different APIs

// Threshold value from dsg.in for Symmetric Eigenvalues and eigenvectors APIs
#define SYM_EIGEN_THRESHOLD 20.0

#endif // MAIN_H
