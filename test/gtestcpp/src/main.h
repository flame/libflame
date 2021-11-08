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
  integer lda; /* LDA is INTEGER
                  The leading dimension of the array A.  LDA >= max(1,N).*/
  integer lwork_heev_2stage; /* LWORK is INTEGER
                                The length of the array WORK. LWORK >= 1, when N <= 1;
                                otherwise  
                                If JOBZ = 'N' and N > 1, LWORK must be queried.
                                                         LWORK = MAX(1, dimension) where
                                                         dimension = max(stage1,stage2) + (KD+1)*N + N
                                                                   = N*KD + N*max(KD+1,FACTOPTNB) 
                                                                     + max(2*KD*KD, KD*NTHREADS) 
                                                                     + (KD+1)*N + N
                                                         where KD is the blocking size of the reduction,
                                                         FACTOPTNB is the blocking used by the QR or LQ
                                                         algorithm, usually FACTOPTNB=128 is a good choice
                                                         NTHREADS is the number of threads used when
                                                         openMP compilation is enabled, otherwise =1.
                                If JOBZ = 'V' and N > 1, LWORK must be queried. Not yet available

                                If LWORK = -1, then a workspace query is assumed; the routine
                                only calculates the optimal size of the WORK array, returns
                                this value as the first entry of the WORK array, and no error
                                message related to LWORK is issued by XERBLA.*/
  integer lwork_heevd_2stage;  /* LWORK is INTEGER
                                  The dimension of the array WORK.
                                  If N <= 1,               LWORK must be at least 1.
                                  If JOBZ = 'N' and N > 1, LWORK must be queried.
                                                           LWORK = MAX(1, dimension) where
                                                           dimension = max(stage1,stage2) + (KD+1)*N + N+1
                                                                     = N*KD + N*max(KD+1,FACTOPTNB) 
                                                                       + max(2*KD*KD, KD*NTHREADS) 
                                                                       + (KD+1)*N + N+1
                                                           where KD is the blocking size of the reduction,
                                                           FACTOPTNB is the blocking used by the QR or LQ
                                                           algorithm, usually FACTOPTNB=128 is a good choice
                                                           NTHREADS is the number of threads used when
                                                           openMP compilation is enabled, otherwise =1.
                                  If JOBZ = 'V' and N > 1, LWORK must be at least 2*N + N**2

                                  If LWORK = -1, then a workspace query is assumed; the routine
                                  only calculates the optimal sizes of the WORK, RWORK and
                                  IWORK arrays, returns these values as the first entries of
                                  the WORK, RWORK and IWORK arrays, and no error message
                                  related to LWORK or LRWORK or LIWORK is issued by XERBLA.*/
  integer lwork_heevr; /* LWORK is INTEGER
                          The length of the array WORK.  LWORK >= max(1,2*N).
                          For optimal efficiency, LWORK >= (NB+1)*N,
                          where NB is the max of the blocksize for CHETRD and for
                          CUNMTR as returned by ILAENV.

                          If LWORK = -1, then a workspace query is assumed; the routine
                          only calculates the optimal sizes of the WORK, RWORK and
                          IWORK arrays, returns these values as the first entries of
                          the WORK, RWORK and IWORK arrays, and no error message
                          related to LWORK or LRWORK or LIWORK is issued by XERBLA.*/
  integer lrwork_heevr;  /* LRWORK is INTEGER
                            The length of the array RWORK.  LRWORK >= max(1,24*N).

                            If LRWORK = -1, then a workspace query is assumed; the
                            routine only calculates the optimal sizes of the WORK, RWORK
                            and IWORK arrays, returns these values as the first entries
                            of the WORK, RWORK and IWORK arrays, and no error message
                            related to LWORK or LRWORK or LIWORK is issued by XERBLA.*/
  
  integer liwork_heevr;  /* LIWORK is INTEGER
                            The dimension of the array IWORK.  LIWORK >= max(1,10*N).

                            If LIWORK = -1, then a workspace query is assumed; the
                            routine only calculates the optimal sizes of the WORK, RWORK
                            and IWORK arrays, returns these values as the first entries
                            of the WORK, RWORK and IWORK arrays, and no error message
                            related to LWORK or LRWORK or LIWORK is issued by XERBLA.*/
  // Added for hegst()
  integer itype; /* ITYPE is INTEGER
                    = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
                    = 2 or 3: compute U*A*U**H or L**H*A*L.*/
  integer ldb; /* LDB is INTEGER
                  The leading dimension of the array B.  LDB >= max(1,N).*/
  integer lwork_hegv;  /* LWORK is INTEGER
                          The length of the array WORK.  LWORK >= max(1,2*N-1).
                          For optimal efficiency, LWORK >= (NB+1)*N,
                          where NB is the blocksize for CHETRD returned by ILAENV.

                          If LWORK = -1, then a workspace query is assumed; the routine
                          only calculates the optimal size of the WORK array, returns
                          this value as the first entry of the WORK array, and no error
                          message related to LWORK is issued by XERBLA.*/
  // Added for heswapr()
  integer i1;  /* I1 is INTEGER
                  Index of the first row to swap*/
  integer i2;  /* I2 is INTEGER
                  Index of the second row to swap*/
  // Added for hetrd()
  integer lwork_hetrd; /* LWORK is INTEGER
                          The dimension of the array WORK.  LWORK >= 1.
                          For optimum performance LWORK >= N*NB, where NB is the
                          optimal blocksize.

                          If LWORK = -1, then a workspace query is assumed; the routine
                          only calculates the optimal size of the WORK array, returns
                          this value as the first entry of the WORK array, and no error
                          message related to LWORK is issued by XERBLA.*/
  // Added for hetrd_2stage()
  integer lhous2;  /* LHOUS2 is INTEGER
                      The dimension of the array HOUS2.
                      If LWORK = -1, or LHOUS2=-1,
                      then a query is assumed; the routine
                      only calculates the optimal size of the HOUS2 array, returns
                      this value as the first entry of the HOUS2 array, and no error
                      message related to LHOUS2 is issued by XERBLA.
                      If VECT='N', LHOUS2 = max(1, 4*n);
                      if VECT='V', option not yet available.*/
  integer lwork_hetrd_2stage;  /* LWORK is INTEGER
                                  The dimension of the array WORK. LWORK = MAX(1, dimension)
                                  If LWORK = -1, or LHOUS2 = -1,
                                  then a workspace query is assumed; the routine
                                  only calculates the optimal size of the WORK array, returns
                                  this value as the first entry of the WORK array, and no error
                                  message related to LWORK is issued by XERBLA.
                                  LWORK = MAX(1, dimension) where
                                  dimension   = max(stage1,stage2) + (KD+1)*N
                                              = N*KD + N*max(KD+1,FACTOPTNB) 
                                                + max(2*KD*KD, KD*NTHREADS) 
                                                + (KD+1)*N 
                                  where KD is the blocking size of the reduction,
                                  FACTOPTNB is the blocking used by the QR or LQ
                                  algorithm, usually FACTOPTNB=128 is a good choice
                                  NTHREADS is the number of threads used when
                                  openMP compilation is enabled, otherwise =1.*/
  char vect_hetrd_2stage;  /* VECT is CHARACTER*1
                              = 'N':  No need for the Housholder representation, 
                                      in particular for the second stage (Band to
                                      tridiagonal) and thus LHOUS2 is of size max(1, 4*N);
                              = 'V':  the Householder representation is needed to 
                                      either generate Q1 Q2 or to apply Q1 Q2, 
                                      then LHOUS2 is to be queried and computed.
                                      (NOT AVAILABLE IN THIS RELEASE).*/
  // Added for hetrd_hb2st
  char stage1; /* STAGE1 is CHARACTER*1
                  = 'N':  "No": to mention that the stage 1 of the reduction  
                          from dense to band using the chetrd_he2hb routine
                          was not called before this routine to reproduce AB. 
                          In other term this routine is called as standalone. 
                  = 'Y':  "Yes": to mention that the stage 1 of the 
                          reduction from dense to band using the chetrd_he2hb 
                          routine has been called to produce AB (e.g., AB is
                          the output of chetrd_he2hb.*/
  integer kd;  /* KD is INTEGER
                  The number of superdiagonals of the matrix A if UPLO = 'U',
                  or the number of subdiagonals if UPLO = 'L'.  KD >= 0.*/
  integer lwork_hetrd_hb2st;   /* LWORK is INTEGER
                                  The dimension of the array WORK. LWORK = MAX(1, dimension)
                                  If LWORK = -1, or LHOUS=-1,
                                  then a workspace query is assumed; the routine
                                  only calculates the optimal size of the WORK array, returns
                                  this value as the first entry of the WORK array, and no error
                                  message related to LWORK is issued by XERBLA.
                                  LWORK = MAX(1, dimension) where
                                  dimension   = (2KD+1)*N + KD*NTHREADS
                                  where KD is the blocking size of the reduction,
                                  FACTOPTNB is the blocking used by the QR or LQ
                                  algorithm, usually FACTOPTNB=128 is a good choice
                                  NTHREADS is the number of threads used when
                                  openMP compilation is enabled, otherwise =1.*/
  integer lhous; /* LHOUS is INTEGER
                    The dimension of the array HOUS. LHOUS = MAX(1, dimension)
                    If LWORK = -1, or LHOUS=-1,
                    then a query is assumed; the routine
                    only calculates the optimal size of the HOUS array, returns
                    this value as the first entry of the HOUS array, and no error
                    message related to LHOUS is issued by XERBLA.
                    LHOUS = MAX(1, dimension) where
                    dimension = 4*N if VECT='N'
                    not available now if VECT='H' */
  // Added for hbtrd_he2hb
  integer lwork_hbtrd_he2hb;   /* LWORK is INTEGER
                                  The dimension of the array WORK which should be calculated
                                  by a workspace query. LWORK = MAX(1, LWORK_QUERY)
                                  If LWORK = -1, then a workspace query is assumed; the routine
                                  only calculates the optimal size of the WORK array, returns
                                  this value as the first entry of the WORK array, and no error
                                  message related to LWORK is issued by XERBLA.
                                  LWORK_QUERY = N*KD + N*max(KD,FACTOPTNB) + 2*KD*KD
                                  where FACTOPTNB is the blocking used by the QR or LQ
                                  algorithm, usually FACTOPTNB=128 is a good choice otherwise
                                  putting LWORK=-1 will provide the size of WORK.*/
} EIG_paramlist;

/* structure to hold Linear solver parameters */
typedef struct Lin_solver_paramlist_t {
  // Added for hecon
  char uplo; /* UPLO is CHARACTER*1
                  Specifies whether the details of the factorization are stored
                  as an upper or lower triangular matrix.
                  = 'U':  Upper triangular, form is A = U*D*U**H;
                  = 'L':  Lower triangular, form is A = L*D*L**H.*/
  integer n;         // The order of matrix A; N >= 0.
  integer lda;       // The leading dimension of the array A. LDA >= max(1,N).
  double anorm; // ANORM is DOUBLE PRECISION. The 1-norm of the original matrix A.
  // Added for hegs2()
  integer itype; /* ITYPE is INTEGER
                    = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
                    = 2 or 3: compute U*A*U**H or L**H *A*L.*/
  integer ldb; /* LDB is INTEGER
                  The leading dimension of the array B.  LDB >= max(1,N).*/
  // Added for herfs()
  integer nrhs;  /* NRHS is INTEGER
                    The number of right hand sides, i.e., the number of columns
                    of the matrices B and X.  NRHS >= 0.*/
  integer ldaf;  /* LDAF is INTEGER
                    The leading dimension of the array AF.  LDAF >= max(1,N).*/
  integer ldx; /* LDX is INTEGER
                  The leading dimension of the array X.  LDX >= max(1,N).*/
  
} Lin_solver_paramlist;

/* structure to hold Linear driver parameters */
typedef struct Lin_driver_paramlist_t {
  // Added for hesv()
  char uplo; /* UPLO is CHARACTER*1
                = 'U':  Upper triangle of A is stored;
                = 'L':  Lower triangle of A is stored.*/
  integer n; /* N is INTEGER
                The number of linear equations, i.e., the order of the
                matrix A.  N >= 0.*/
  integer nrhs;  /* NRHS is INTEGER
                    The number of right hand sides, i.e., the number of columns
                    of the matrix B.  NRHS >= 0.*/
  integer lda; /* LDA is INTEGER
                  The leading dimension of the array A.  LDA >= max(1,N).*/
  integer ldb; /* LDB is INTEGER
                  The leading dimension of the array B.  LDB >= max(1,N).*/
  integer lwork; /* LWORK is INTEGER
                    The length of WORK.  LWORK >= 1, and for best performance
                    LWORK >= max(1,N*NB), where NB is the optimal blocksize for
                    CHETRF.
                    for LWORK < N, TRS will be done with Level BLAS 2
                    for LWORK >= N, TRS will be done with Level BLAS 3

                    If LWORK = -1, then a workspace query is assumed; the routine
                    only calculates the optimal size of the WORK array, returns
                    this value as the first entry of the WORK array, and no error
                    message related to LWORK is issued by XERBLA.*/
  // Added for hesv_aa()
  integer lwork_hesv_aa; /* LWORK is INTEGER
                            The length of WORK.  LWORK >= MAX(1,2*N,3*N-2), and for best 
                            performance LWORK >= MAX(1,N*NB), where NB is the optimal
                            blocksize for CHETRF.

                            If LWORK = -1, then a workspace query is assumed; the routine
                            only calculates the optimal size of the WORK array, returns
                            this value as the first entry of the WORK array, and no error
                            message related to LWORK is issued by XERBLA.*/
  // Added for hesv_aa_2stage()
  integer ltb; /* LTB is INTEGER
                  The size of the array TB. LTB >= 4*N, internally
                  used to select NB such that LTB >= (3*NB+1)*N.

                  If LTB = -1, then a workspace query is assumed; the
                  routine only calculates the optimal size of LTB, 
                  returns this value as the first entry of TB, and
                  no error message related to LTB is issued by XERBLA.*/
  integer lwork_hesv_aa_2stage;  /* LWORK is INTEGER
                                    The size of WORK. LWORK >= N, internally used to select NB
                                    such that LWORK >= N*NB.

                                    If LWORK = -1, then a workspace query is assumed; the
                                    routine only calculates the optimal size of the WORK array,
                                    returns this value as the first entry of the WORK array, and
                                    no error message related to LWORK is issued by XERBLA.*/
  // Added for hesvx()
  char fact;   /* FACT is CHARACTER*1
                  Specifies whether or not the factored form of A has been
                  supplied on entry.
                  = 'F':  On entry, AF and IPIV contain the factored form
                          of A.  A, AF and IPIV will not be modified.
                  = 'N':  The matrix A will be copied to AF and factored.*/
  integer ldaf; /* LDAF is INTEGER
                   The leading dimension of the array AF.  LDAF >= max(1,N).*/
  integer ldx; /* LDX is INTEGER
                   The leading dimension of the array X.  LDX >= max(1,N).*/
  // Added for hesvxx()
  char fact_hesvxx; /* FACT is CHARACTER*1
                       Specifies whether or not the factored form of the matrix A is
                       supplied on entry, and if not, whether the matrix A should be
                       equilibrated before it is factored.
                         = 'F':  On entry, AF and IPIV contain the factored form of A.
                                 If EQUED is not 'N', the matrix A has been
                                 equilibrated with scaling factors given by S.
                                 A, AF, and IPIV are not modified.
                         = 'N':  The matrix A will be copied to AF and factored.
                         = 'E':  The matrix A will be equilibrated if necessary, then
                                 copied to AF and factored.*/
  char equed; /* EQUED is CHARACTER*1
                 Specifies the form of equilibration that was done.
                   = 'N':  No equilibration (always true if FACT = 'N').
                   = 'Y':  Both row and column equilibration, i.e., A has been
                           replaced by diag(S) * A * diag(S).
                 EQUED is an input argument if FACT = 'F'; otherwise, it is an
                 output argument.*/
  integer n_err_bnds; /* N_ERR_BNDS is INTEGER
                         Number of error bounds to return for each right hand side
                         and each type (normwise or componentwise).  See ERR_BNDS_NORM and
                         ERR_BNDS_COMP below.*/
  integer nparams;  /* NPARAMS is INTEGER
                       Specifies the number of parameters set in PARAMS.  If <= 0, the
                       PARAMS array is never referenced and default values are used.*/
  
} Lin_driver_paramlist;

extern EIG_paramlist eig_paramslist[NUM_SUB_TESTS];
extern Lin_solver_paramlist lin_solver_paramslist[NUM_SUB_TESTS];
extern Lin_driver_paramlist lin_driver_paramslist[NUM_SUB_TESTS];

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

/* Macro to store seed value for srand().
   1 - default value to get deterministic random values on multiple runs
   of rand()*/
#define SRAND_SEED_VALUE 1

// Macro to enable status/error messages of test APIs.
#define PRINT_MSGS 0

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

// Threshold value from glm.in for Generalized linear driver APIs
#define LIN_DRVR_THRESHOLD 20.0

#endif // MAIN_H