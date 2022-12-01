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
  GEQRF computes a QR factorization of a M-by-N matrix A: A = Q * R.
*/

// GEQRF and GEQR2
#define LAPACK_geqrf(prefix)                                            \
  int F77_ ## prefix ## geqrf(int* m,                                   \
                              int* n,                                   \
                              PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                              PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,          \
                              PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                              int* info )

#define LAPACK_geqrf_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Obj      A, t, T;                                         \
  int          min_m_n  = min( *m, *n );                        \
  FLA_Error    init_result;                                     \
                                                                \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &t );            \
  FLA_Obj_attach_buffer( buff_t, 1, min_m_n, &t );                      \
                                                                        \
  FLA_Set( FLA_ZERO, t );                                               \
                                                                        \
  FLA_QR_UT_create_T( A, &T );                                          \
  FLA_QR_UT( A, T );                                                    \
  FLA_QR_UT_recover_tau( T, t );                                        \
  PREFIX2FLAME_INVERT_TAU(prefix,t);                                    \
                                                                        \
  FLA_Obj_free_without_buffer( &A );                                    \
  FLA_Obj_free_without_buffer( &t );                                    \
  FLA_Obj_free( &T );                                                   \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;

LAPACK_geqrf(s)
{
    {
        LAPACK_RETURN_CHECK( sgeqrf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_geqrf_body(s)
    }
}
LAPACK_geqrf(d)
{
    {
        LAPACK_RETURN_CHECK( dgeqrf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_geqrf_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqrf(c)
{
    {
        LAPACK_RETURN_CHECK( cgeqrf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_geqrf_body(c)
    }
}
LAPACK_geqrf(z)
{
    {
        LAPACK_RETURN_CHECK( zgeqrf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_geqrf_body(z)
    }
}
#endif

#define LAPACK_geqr2(prefix)                                            \
  int F77_ ## prefix ## geqr2( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,         \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,         \
                               int* info )

LAPACK_geqr2(s)
{
    {
        LAPACK_RETURN_CHECK( sgeqr2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_geqrf_body(s)
    }
}
LAPACK_geqr2(d)
{
    {
        LAPACK_RETURN_CHECK( dgeqr2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_geqrf_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqr2(c)
{
    {
        LAPACK_RETURN_CHECK( cgeqr2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_geqrf_body(c)
    }
}
LAPACK_geqr2(z)
{
    {
        LAPACK_RETURN_CHECK( zgeqr2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_geqrf_body(z)
    }
}
#endif

// GEQRFP and GEQR2P
#define LAPACK_geqrfp(prefix)                                            \
  int F77_ ## prefix ## geqrfp(int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                               int* info )
LAPACK_geqrfp(s)
{
    {
        LAPACK_RETURN_CHECK( sgeqrfp_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w, lwork,
                                            info ) )
    }
    {
        LAPACK_geqrf_body(s)
    }
}
LAPACK_geqrfp(d)
{
    {
        LAPACK_RETURN_CHECK( dgeqrfp_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w, lwork,
                                            info ) )
    }
    {
        LAPACK_geqrf_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqrfp(c)
{
    {
        LAPACK_RETURN_CHECK( cgeqrfp_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w, lwork,
                                            info ) )
    }
    {
        LAPACK_geqrf_body(c)
    }
}
LAPACK_geqrfp(z)
{
    {
        LAPACK_RETURN_CHECK( zgeqrfp_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w, lwork,
                                            info ) )
    }
    {
        LAPACK_geqrf_body(z)
    }
}
#endif

#define LAPACK_geqr2p(prefix)                                            \
  int F77_ ## prefix ## geqr2p( int* m,                                  \
                                int* n,                                 \
                                PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                                PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,  \
                                PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,  \
                                int* info )

LAPACK_geqr2p(s)
{
    {
        LAPACK_RETURN_CHECK( sgeqr2p_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w,
                                            info ) )
    }
    {
        LAPACK_geqrf_body(s)
    }
}
LAPACK_geqr2p(d)
{
    {
        LAPACK_RETURN_CHECK( dgeqr2p_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w,
                                            info ) )
    }
    {
        LAPACK_geqrf_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqr2p(c)
{
    {
        LAPACK_RETURN_CHECK( cgeqr2p_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w,
                                            info ) )
    }
    {
        LAPACK_geqrf_body(c)
    }
}
LAPACK_geqr2p(z)
{
    {
        LAPACK_RETURN_CHECK( zgeqr2p_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w,
                                            info ) )
    }
    {
        LAPACK_geqrf_body(z)
    }
}
#endif

#endif
