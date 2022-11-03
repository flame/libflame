/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_YET_LAPACK2FLAME

#include "FLASH_lapack2flash_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  GEHRD reduces a general matrix A to upper Hessenberg form H by
  an unitary similarity transformation: Q**H * A * Q = H .
*/

#define LAPACK_gehrd(prefix)                                            \
  int F77_ ## prefix ## gehrd( int* m,                                  \
                               int* ilo,                                \
                               int* ihi,                                \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                               int* info )

#define LAPACK_gehrd_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  dim_t        m_t      = ( *m - 1 );                           \
  FLA_Obj      A, t, T;                                         \
  FLA_Error    init_result;                                     \
                                                                \
  FLA_Init_safe( &init_result );                                \
                                                                \
  FLA_Obj_create_without_buffer( datatype, *m, *m, &A );        \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );              \
                                                                \
  FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );        \
  if ( m_t > 0 ) FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );   \
                                                                \
  FLA_Hess_UT_create_T( A, &T );                                \
  FLA_Hess_UT( A, T );                                          \
  FLA_Hess_UT_recover_tau( T, t );                              \
  FLA_Obj_free( &T );                                           \
                                                                \
  PREFIX2FLAME_INVERT_TAU(prefix,t);                            \
  FLA_Obj_free_without_buffer( &t );                            \
                                                                \
  FLA_Obj_free_without_buffer( &A );                            \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \
  *info = 0;                                                    \
                                                                \
  return 0;


LAPACK_gehrd(s)
{
    {
        LAPACK_RETURN_CHECK( sgehrd_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gehrd_body(s)
    }
}
LAPACK_gehrd(d)
{
    {
        LAPACK_RETURN_CHECK( dgehrd_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gehrd_body(d)
    }
}
#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gehrd(c)
{
    {
        LAPACK_RETURN_CHECK( cgehrd_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gehrd_body(c)
    }
}
LAPACK_gehrd(z)
{
    {
        LAPACK_RETURN_CHECK( zgehrd_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gehrd_body(z)
    }
}
#endif

#define LAPACK_gehd2(prefix)                                            \
  int F77_ ## prefix ## gehd2( int* m,                                  \
                               int* ilo,                                \
                               int* ihi,                                \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,   \
                               int* info )

LAPACK_gehd2(s)
{
    {
        LAPACK_RETURN_CHECK( sgehd2_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gehrd_body(s)
    }
}
LAPACK_gehd2(d)
{
    {
        LAPACK_RETURN_CHECK( dgehd2_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gehrd_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gehd2(c)
{
    {
        LAPACK_RETURN_CHECK( cgehd2_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gehrd_body(c)
    }
}
LAPACK_gehd2(z)
{
    {
        LAPACK_RETURN_CHECK( zgehd2_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gehrd_body(z)
    }
}
#endif

#endif
