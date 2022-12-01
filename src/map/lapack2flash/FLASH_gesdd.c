/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLASH

#include "FLASH_lapack2flash_util_defs.h"
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
                               int*  m,                                 \
                               int*  n,                                 \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A,  int* ldim_A, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_s,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_U,  int* ldim_U, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_Vh, int* ldim_Vh, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,  int *lwork, \
                               int *buff_i,                             \
                               int *info )

#define LAPACK_gesdd_complex(prefix)                                    \
  int F77_ ## prefix ## gesdd( char* jobz,                              \
                               int*  m,                                 \
                               int*  n,                                 \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A,  int* ldim_A, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_s,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_U,  int* ldim_U, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_Vh, int* ldim_Vh, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,  int* lwork, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_r,   \
                               int *buff_i,                             \
                               int *info )

#define LAPACK_gesdd_real_body(prefix)                                  \
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
                           info );                                      \
  return 0;

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
  return 0;

LAPACK_gesdd_real(s)
{
    {
        LAPACK_RETURN_CHECK( sgesdd_check( jobz,
                                           m, n,
                                           buff_A,  ldim_A,
                                           buff_s,
                                           buff_U,  ldim_U,
                                           buff_Vh, ldim_Vh,
                                           buff_w,  lwork,
                                           buff_i,
                                           info ) )
    }
    {
        LAPACK_gesdd_real_body(s)
    }
}
LAPACK_gesdd_real(d)
{
    {
        LAPACK_RETURN_CHECK( dgesdd_check( jobz,
                                           m, n,
                                           buff_A,  ldim_A,
                                           buff_s,
                                           buff_U,  ldim_U,
                                           buff_Vh, ldim_Vh,
                                           buff_w,  lwork,
                                           buff_i,
                                           info ) )
    }
    {
        LAPACK_gesdd_real_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gesdd_complex(c)
{
    {
        LAPACK_RETURN_CHECK( cgesdd_check( jobz,
                                           m, n,
                                           buff_A,  ldim_A,
                                           buff_s,
                                           buff_U,  ldim_U,
                                           buff_Vh, ldim_Vh,
                                           buff_w,  lwork,
                                           buff_r,
                                           buff_i,
                                           info ) );
    }
    {
        LAPACK_gesdd_complex_body(c)
    }
}
LAPACK_gesdd_complex(z)
{
    {
        LAPACK_RETURN_CHECK( zgesdd_check( jobz,
                                           m, n,
                                           buff_A,  ldim_A,
                                           buff_s,
                                           buff_U,  ldim_U,
                                           buff_Vh, ldim_Vh,
                                           buff_w,  lwork,
                                           buff_r,
                                           buff_i,
                                           info ) );
    }
    {
        LAPACK_gesdd_complex_body(z)
    }
}
#endif

#endif
