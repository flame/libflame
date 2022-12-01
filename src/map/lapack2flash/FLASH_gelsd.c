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
  GELSD computes the minimum-norm solution to a real linear least
  squares problem:
  minimize 2-norm(| b - A*x |)
  using the singular value decomposition (SVD) of A. A is an M-by-N
  matrix which may be rank-deficient.

  Several right hand side vectors b and solution vectors x can be
  handled in a single call; they are stored as the columns of the
  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
  matrix X.

  The problem is solved in three steps:
  (1) Reduce the coefficient matrix A to bidiagonal form with
  Householder transformations, reducing the original problem
  into a "bidiagonal least squares problem" (BLS)
  (2) Solve the BLS using a divide and conquer approach.
  (3) Apply back all the Householder tranformations to solve
  the original least squares problem.

  The effective rank of A is determined by treating as zero those
  singular values which are less than RCOND times the largest singular
  value.

  The divide and conquer algorithm makes very mild assumptions about
  floating point arithmetic. It will work on machines with a guard
  digit in add/subtract, or on those binary machines without guard
  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
  Cray-2. It could conceivably fail on hexadecimal or decimal machines
  without guard digits, but we know of none.


  At this moment, this routine is redirected to GELSS.
*/

#define LAPACK_gelsd_real(prefix)                                       \
  int F77_ ## prefix ## gelsd( int* m,                                  \
                               int* n,                                  \
                               int* nrhs,                               \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, int* ldim_B, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_s,   \
                               PREFIX2LAPACK_REALDEF(prefix)* rcond,    \
                               int* rank,                               \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                               int* iwork,                              \
                               int* info )

#define LAPACK_gelsd_real_body(prefix)                                  \
  F77_ ## prefix ## gelss( m, n, nrhs,                                  \
                           buff_A, ldim_A,                              \
                           buff_B, ldim_B,                              \
                           buff_s, rcond, rank,                         \
                           buff_w, lwork, info);                        \
  return 0;


LAPACK_gelsd_real(s)
{
    {
        LAPACK_RETURN_CHECK( sgelsd_check( m, n, nrhs,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           buff_s, rcond, rank,
                                           buff_w, lwork,
                                           iwork,
                                           info ) )
    }
    {
        LAPACK_gelsd_real_body(s)
    }
}
LAPACK_gelsd_real(d)
{
    {
        LAPACK_RETURN_CHECK( dgelsd_check( m, n, nrhs,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           buff_s, rcond, rank,
                                           buff_w, lwork,
                                           iwork,
                                           info ) )
    }
    {
        LAPACK_gelsd_real_body(d)
    }
}

#define LAPACK_gelsd_complex(prefix)                                    \
  int F77_ ## prefix ## gelsd( int* m,                                  \
                               int* n,                                  \
                               int* nrhs,                               \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, int* ldim_B, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_s,   \
                               PREFIX2LAPACK_REALDEF(prefix)* rcond,    \
                               int* rank,                               \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_r,   \
                               int* iwork,                              \
                               int* info )

#define LAPACK_gelsd_complex_body(prefix)                               \
  F77_ ## prefix ## gelss( m, n, nrhs,                                  \
                           buff_A, ldim_A,                              \
                           buff_B, ldim_B,                              \
                           buff_s, rcond, rank,                         \
                           buff_w, lwork, buff_r, info);                \
  return 0;

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gelsd_complex(c)
{
    {
        LAPACK_RETURN_CHECK( cgelsd_check( m, n, nrhs,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           buff_s, rcond, rank,
                                           buff_w, lwork,
                                           buff_r, iwork,
                                           info ) )
    }
    {
        LAPACK_gelsd_complex_body(c)
    }
}
LAPACK_gelsd_complex(z)
{
    {
        LAPACK_RETURN_CHECK( zgelsd_check( m, n, nrhs,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           buff_s, rcond, rank,
                                           buff_w, lwork,
                                           buff_r, iwork,
                                           info ) )
    }
    {
        LAPACK_gelsd_complex_body(z)
    }
}
#endif

#endif
