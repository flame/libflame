/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
/*
 *  Copyright (c) 2020-2023 Advanced Micro Devices, Inc.  All rights reserved.
 */
#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  GESDD computes the singular value decomposition (SVD) of a
  M-by-N matrix A, optionally computing the left and right singular
  vectors.  If singular vectors are desired, it uses a
  divide-and-conquer algorithm.

  The SVD is written
      A = U * SIGMA * transpose(V)
  where SIGMA is an M-by-N matrix which is zero except for its
  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
  are the singular values of A; they are real and non-negative, and
  are returned in descending order.  The first min(m,n) columns of
  U and V are the left and right singular vectors of A.

  Note that the routine returns VT = V**T, not V.

  At this moment, this routine is redirected to GESVD.
*/

#define LAPACK_gesdd_real(prefix)                                       \
  int F77_ ## prefix ## gesdd( char* jobz,                              \
                               integer*  m,                                 \
                               integer*  n,                                 \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A,  integer* ldim_A, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_s,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_U,  integer* ldim_U, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_Vh, integer* ldim_Vh, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,  integer *lwork, \
                               integer *buff_i,                             \
                               integer *info )

#define LAPACK_gesdd_complex(prefix)                                    \
  int F77_ ## prefix ## gesdd( char* jobz,                              \
                               integer*  m,                                 \
                               integer*  n,                                 \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A,  integer* ldim_A, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_s,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_U,  integer* ldim_U, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_Vh, integer* ldim_Vh, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,  integer* lwork, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_r,   \
                               integer *buff_i,                             \
                               integer *info )

#define LAPACK_gesdd_real_body(prefix)                                  \
                                                                        \
  F77_ ## prefix ## gesvd( jobu, jobv,                                  \
                           m, n,                                        \
                           buff_A,  ldim_A,                             \
                           buff_s,                                      \
                           buff_U,  ldim_U,                             \
                           buff_Vh, ldim_Vh,                            \
                           buff_w,  lwork,                              \
                           info );                                      \


#define LAPACK_gesdd_complex_body(prefix)                               \
  char jobu[1], jobv[1];                                                \
                                                                        \
  if ( *jobz == 'O' ) {                                                 \
    if ( *m >= *n ) {                                                   \
      jobu[0] = 'O'; jobv[0] = 'A';                                     \
    } else {                                                            \
      jobu[0] = 'A'; jobv[0] = 'O';                                     \
    }                                                                   \
  } else {                                                              \
    jobu[0] = *jobz; jobv[0] = *jobz;                                   \
  }                                                                     \
                                                                        \
  F77_ ## prefix ## gesvd( jobu, jobv,                                  \
                           m, n,                                        \
                           buff_A,  ldim_A,                             \
                           buff_s,                                      \
                           buff_U,  ldim_U,                             \
                           buff_Vh, ldim_Vh,                            \
                           buff_w,  lwork,                              \
                           buff_r,                                      \
                           info );                                      \


LAPACK_gesdd_real(s)
{
  int fla_error = LAPACK_SUCCESS;
  AOCL_DTL_TRACE_LOG_INIT
  AOCL_DTL_SNPRINTF("sgesdd inputs: jobz %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldu %" FLA_IS ", ldvt %" FLA_IS "", *jobz, *m, *n, *ldim_A, *ldim_U, *ldim_Vh);
  extern int sgesdd_fla_check(char *jobu, char *jobvt, integer *m, integer *n, float *a, integer *lda, float *s, float *u, integer *ldu, float *vt, integer *ldvt, float *work, integer *lwork, integer *info);
  extern int lapack_sgesdd(char *jobz, integer *m, integer *n, real *a, integer *lda, real *s, real *u, integer *ldu, real *vt, integer *ldvt, real *work, integer *lwork, integer *iwork, integer *info);
 
#if FLA_AMD_OPT
  {
    if(*m > 750 || *n > 750){
      char jobu[1], jobv[1];

      if ( *jobz == 'O' ) {
        if ( *m >= *n ) {
          jobu[0] = 'O'; jobv[0] = 'A';
        }
        else {
          jobu[0] = 'A'; jobv[0] = 'O';
        }
      }
      else {
        jobu[0] = *jobz; jobv[0] = *jobz;
      }
      LAPACK_RETURN_CHECK_VAR1( sgesdd_fla_check( jobu, jobv,
                                          m, n,
                                          buff_A,  ldim_A,
                                          buff_s,
                                          buff_U,  ldim_U,
                                          buff_Vh, ldim_Vh,
                                          buff_w,  lwork,
                                          info ),fla_error )
      if(fla_error==LAPACK_SUCCESS)
      {
          LAPACK_gesdd_real_body(s)
          /** fla_error set to 0 on LAPACK_SUCCESS */
          fla_error = 0;
      }
    }
    else{
      LAPACK_RETURN_CHECK_VAR1( sgesdd_check( jobz,
                                      m, n,
                                      buff_A, ldim_A,
                                      buff_s,
                                      buff_U, ldim_U,
                                      buff_Vh, ldim_Vh,
                                      buff_w, lwork,
                                      buff_i, info), fla_error )

      if (fla_error == LAPACK_SUCCESS) {
        lapack_sgesdd(jobz,
                m, n,
                buff_A, ldim_A,
                buff_s,
                buff_U, ldim_U,
                buff_Vh, ldim_Vh,
                buff_w, lwork,
                buff_i, info);
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
      }
    }
  }
#else
  {
    char jobu[1], jobv[1];

    if ( *jobz == 'O' ) {
      if ( *m >= *n ) {
        jobu[0] = 'O'; jobv[0] = 'A';
      }
      else {
        jobu[0] = 'A'; jobv[0] = 'O';
      }
    }
    else {
      jobu[0] = *jobz; jobv[0] = *jobz;
    }
    LAPACK_RETURN_CHECK_VAR1( sgesdd_fla_check( jobu, jobv,
                                        m, n,
                                        buff_A,  ldim_A,
                                        buff_s,
                                        buff_U,  ldim_U,
                                        buff_Vh, ldim_Vh,
                                        buff_w,  lwork,
                                        info ),fla_error )
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_gesdd_real_body(s)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
  }
#endif

  AOCL_DTL_TRACE_LOG_EXIT
  return fla_error;
}


LAPACK_gesdd_real(d)
{
  int fla_error = LAPACK_SUCCESS;
  AOCL_DTL_TRACE_LOG_INIT
  AOCL_DTL_SNPRINTF("dgesdd inputs: jobz %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldu %" FLA_IS ", ldvt %" FLA_IS "", *jobz, *m, *n, *ldim_A, *ldim_U, *ldim_Vh);
  extern int lapack_dgesdd(char *jobz, integer *m, integer *n, doublereal * a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *iwork, integer *info);

#if FLA_AMD_OPT
  {
    LAPACK_RETURN_CHECK_VAR1( dgesdd_check( jobz,
                                    m, n,
                                    buff_A, ldim_A,
                                    buff_s,
                                    buff_U, ldim_U,
                                    buff_Vh, ldim_Vh,
                                    buff_w, lwork,
                                    buff_i, info), fla_error )

    if (fla_error == LAPACK_SUCCESS) {

      /* Initialize global context data */
      aocl_fla_init();

      lapack_dgesdd(jobz,
            m, n,
            buff_A, ldim_A,
            buff_s,
            buff_U, ldim_U,
            buff_Vh, ldim_Vh,
            buff_w, lwork,
            buff_i, info);
      /** fla_error set to 0 on LAPACK_SUCCESS */
      fla_error = 0;
    }
  }
#else
  {
    char jobu[1], jobv[1];

    if ( *jobz == 'O' ) {
      if ( *m >= *n ) {
        jobu[0] = 'O'; jobv[0] = 'A';
      }
      else {
        jobu[0] = 'A'; jobv[0] = 'O';
      }
    }
    else {
      jobu[0] = *jobz; jobv[0] = *jobz;
    }
    LAPACK_RETURN_CHECK_VAR1( dgesdd_fla_check( jobu, jobv,
                                    m, n,
                                    buff_A,  ldim_A,
                                    buff_s,
                                    buff_U,  ldim_U,
                                    buff_Vh, ldim_Vh,
                                    buff_w,  lwork,
                                    info ), fla_error )

    if (fla_error == LAPACK_SUCCESS) {
      LAPACK_gesdd_real_body(d)
      /** fla_error set to 0 on LAPACK_SUCCESS */
      fla_error = 0;
    }
  }
#endif

  AOCL_DTL_TRACE_LOG_EXIT
  return fla_error;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gesdd_complex(c)
{
  int fla_error = LAPACK_SUCCESS;
  AOCL_DTL_TRACE_LOG_INIT
  AOCL_DTL_SNPRINTF("cgesdd inputs: jobz %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldu %" FLA_IS ", ldvt %" FLA_IS "", *jobz, *m, *n, *ldim_A, *ldim_U, *ldim_Vh);
  {
    LAPACK_RETURN_CHECK_VAR1(cgesdd_check(jobz,
                                          m, n,
                                          buff_A, ldim_A,
                                          buff_s,
                                          buff_U, ldim_U,
                                          buff_Vh, ldim_Vh,
                                          buff_w, lwork,
                                          buff_r,
                                          buff_i,
                                          info),
                             fla_error);
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gesdd_complex_body(c)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_gesdd_complex(z)
{
  int fla_error = LAPACK_SUCCESS;
  AOCL_DTL_TRACE_LOG_INIT
  AOCL_DTL_SNPRINTF("zgesdd inputs: jobz %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldu %" FLA_IS ", ldvt %" FLA_IS "", *jobz, *m, *n, *ldim_A, *ldim_U, *ldim_Vh);
  {
    LAPACK_RETURN_CHECK_VAR1(zgesdd_check(jobz,
                                          m, n,
                                          buff_A, ldim_A,
                                          buff_s,
                                          buff_U, ldim_U,
                                          buff_Vh, ldim_Vh,
                                          buff_w, lwork,
                                          buff_r,
                                          buff_i,
                                          info),
                             fla_error);
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gesdd_complex_body(z)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#endif
