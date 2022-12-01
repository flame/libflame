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
   GELQF computes an LQ factorization of a M-by-N matrix A: A = L * Q.
*/

#define LAPACK_gelqf(prefix)                                            \
  int F77_ ## prefix ## gelqf( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                               int* info )

#define LAPACK_gelqf_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Obj      A, t, T;                                         \
  int          min_m_n  = min( *m, *n );                        \
  FLA_Error    init_result;                                     \
                                                                \
  FLA_Init_safe( &init_result );                                \
                                                                \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );        \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );              \
                                                                \
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &t );    \
  FLA_Obj_attach_buffer( buff_t, 1, min_m_n, &t );              \
                                                                \
  FLA_Set( FLA_ZERO, t );                                       \
                                                                \
  FLA_LQ_UT_create_T( A, &T );                                  \
  FLA_LQ_UT( A, T );                                            \
  FLA_LQ_UT_recover_tau( T, t );                                \
  PREFIX2FLAME_INVERT_TAU(prefix,t);                            \
                                                                \
  FLA_Obj_free_without_buffer( &A );                            \
  FLA_Obj_free_without_buffer( &t );                            \
  FLA_Obj_free( &T );                                           \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \
  *info = 0;                                                    \
                                                                \
  return 0;

LAPACK_gelqf(s)
{
    {
        LAPACK_RETURN_CHECK( sgelqf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(s)
    }
}
LAPACK_gelqf(d)
{
    {
        LAPACK_RETURN_CHECK( dgelqf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gelqf(c)
{
    {
        LAPACK_RETURN_CHECK( cgelqf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(c)
    }
}
LAPACK_gelqf(z)
{
    {
        LAPACK_RETURN_CHECK( zgelqf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(z)
    }
}
#endif

#define LAPACK_gelq2(prefix)                                            \
  int F77_ ## prefix ## gelq2( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,   \
                               int* info )
LAPACK_gelq2(s)
{
    {
        LAPACK_RETURN_CHECK( sgelq2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(s)
    }
}
LAPACK_gelq2(d)
{
    {
        LAPACK_RETURN_CHECK( dgelq2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gelq2(c)
{
    {
        LAPACK_RETURN_CHECK( cgelq2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(c)
    }
}
LAPACK_gelq2(z)
{
    {
        LAPACK_RETURN_CHECK( zgelq2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gelqf_body(z)
    }
}
#endif


#endif
