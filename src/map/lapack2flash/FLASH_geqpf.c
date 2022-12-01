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
  - GEQPF computes a QR factorization with column pivoting of a
  real M-by-N matrix A: A*P = Q*R.
  This routine is deprecated and has been replaced by routine GEQP3.

  GEQPF is used for generalized eigenvalue decomposition;
  JPVT is not only output but also input. Tricky to handle this
  in flame pivoting. So, use the original implementation. It still
  use our QR.

  - GEQP3 computes a QR factorization with column pivoting of a
  matrix A: A*P = Q*R using Level 3 BLAS.

  FLA_QR_UT_piv uses level 3 BLAS.
*/

// GEQPF
#define LAPACK_geqpf(prefix)                                            \
  int F77_ ## prefix ## geqpf( int* m,                                   \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               int* buff_p,                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,   \
                               int* info )

// Notation for LAPACK column pvioting is not consistent to pivoting in LU.
// This does not perform pre-ordering when jpiv include non-zero pivots.
//
#define LAPACK_geqpf_body(prefix)                                       \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                \
  FLA_Obj      A, t, T, w, p, jpiv;                                     \
  dim_t        min_m_n  = min( *m, *n );                                \
  FLA_Error    init_result;                                             \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &t );            \
  FLA_Obj_attach_buffer( buff_t, 1, min_m_n, &t );                      \
  FLA_Set( FLA_ZERO, t );                                               \
                                                                        \
  FLA_Obj_create_without_buffer( PREFIX2FLAME_REALTYPE(prefix), *n, 1, &w ); \
  FLA_Obj_attach_buffer( buff_w, 1, *n, &w );                           \
                                                                        \
  /* LAPACK pivot storage */                                            \
  FLA_Obj_create_without_buffer( FLA_INT, *n, 1, &jpiv );               \
  FLA_Obj_attach_buffer( buff_p, 1, *n, &jpiv );                        \
                                                                        \
  /* FLAME pivot storage */                                             \
  FLA_Obj_create( FLA_INT, *n, 1, 0, 0, &p );                           \
  FLA_Set( FLA_ZERO, p );                                               \
                                                                        \
  /* QR_UT_piv */                                                       \
  FLA_QR_UT_create_T( A, &T );                                          \
  FLA_Set( FLA_ZERO, T );                                               \
                                                                        \
  FLA_QR_UT_piv( A, T, w, p );                                          \
  FLA_QR_UT_recover_tau( T, t );                                        \
  PREFIX2FLAME_INVERT_TAU(prefix,t);                                    \
                                                                        \
  /* Transform FLAME column pivots to LAPACK pivots */                  \
  FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, jpiv );              \
                                                                        \
  /* Cleaning */                                                        \
  FLA_Obj_free_without_buffer( &A );                                    \
  FLA_Obj_free_without_buffer( &t );                                    \
  FLA_Obj_free_without_buffer( &w );                                    \
  FLA_Obj_free_without_buffer( &jpiv );                                 \
  FLA_Obj_free( &p );                                                   \
  FLA_Obj_free( &T );                                                   \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;


LAPACK_geqpf(s)
{
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK( sgeqpf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_geqpf_body(s)
    }
}
LAPACK_geqpf(d)
{
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK( dgeqpf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_geqpf_body(d)
    }
}

#define LAPACK_geqpf_complex(prefix)                                    \
  int F77_ ## prefix ## geqpf( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               int* buff_p,                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,   \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_r,   \
                               int* info )

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqpf_complex(c)
{
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK( cgeqpf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w,
                                           buff_r,
                                           info ) )
    }
    {
        LAPACK_geqpf_body(c)
    }
}
LAPACK_geqpf_complex(z)
{
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK( zgeqpf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w,
                                           buff_r,
                                           info ) )
    }
    {
        LAPACK_geqpf_body(z)
    }
}
#endif

// GEQP3
#define LAPACK_geqp3(prefix)                                            \
  int F77_ ## prefix ## geqp3( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               int* buff_p,                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                               int* info )
LAPACK_geqp3(s)
{
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK( sgeqp3_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_geqpf_body(s)
    }
}
LAPACK_geqp3(d)
{
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK( dgeqp3_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_geqpf_body(d)
    }
}


#define LAPACK_geqp3_complex(prefix)                                            \
  int F77_ ## prefix ## geqp3( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               int* buff_p,                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_r,   \
                               int* info )


#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqp3_complex(c)
{
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK( cgeqp3_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w, lwork,
                                           buff_r,
                                           info ) )
    }
    {
        LAPACK_geqpf_body(c)
    }
}
LAPACK_geqp3_complex(z)
{
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK( zgeqp3_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w, lwork,
                                           buff_r,
                                           info ) )
    }
    {
        LAPACK_geqpf_body(z)
    }
}
#endif

#endif
