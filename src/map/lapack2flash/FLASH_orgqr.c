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
  SORGQR generates an M-by-N real matrix Q with orthonormal columns,
  which is defined as the first N columns of a product of K elementary
  reflectors of order M

  Q = H(1) H(2) . . . H(k)

  as returned by SGEQRF.
*/

#define LAPACK_orgqr(prefix, name)                                      \
  int F77_ ## prefix ## name ## qr( int* m,                             \
                                    int* n,                             \
                                    int* k,                             \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                                    int* ldim_A,                        \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, \
                                    int* lwork,                         \
                                    int* info)

#define LAPACK_orgqr_body(prefix)                                       \
  FLA_Datatype datatype   = PREFIX2FLAME_DATATYPE(prefix);              \
  FLA_Obj      A, AL, AR, t, T;                                         \
  FLA_Error    init_result;                                             \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  if ( *k > 0 && !( PREFIX2FLAME_IS_ZERO(prefix, buff_t) ) )        \
    {                                                                   \
      FLA_Obj_create_without_buffer( datatype, *k, 1, &t );             \
      FLA_Obj_attach_buffer( buff_t, 1, *k, &t );                       \
      PREFIX2FLAME_INVERT_TAU(prefix,t);                                \
                                                                        \
      FLA_Part_1x2( A, &AL, &AR, *k, FLA_LEFT );                        \
      FLA_QR_UT_create_T( AL, &T );                                     \
      FLA_Set( FLA_ZERO, T );                                           \
      FLA_Accum_T_UT( FLA_FORWARD, FLA_COLUMNWISE, AL, t, T );          \
      FLA_QR_UT_form_Q( AL, T, A );                                     \
                                                                        \
      PREFIX2FLAME_INVERT_TAU(prefix,t);                                \
      FLA_Obj_free_without_buffer( &t );                                \
      FLA_Obj_free( &T );                                               \
    }                                                                   \
  else                                                                  \
    {                                                                   \
      FLA_Set_to_identity( A );                                         \
    }                                                                   \
  FLA_Obj_free_without_buffer( &A );                                    \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;

LAPACK_orgqr(s, org)
{
    {
        LAPACK_RETURN_CHECK( sorgqr_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orgqr_body(s)
    }
}
LAPACK_orgqr(d, org)
{
    {
        LAPACK_RETURN_CHECK( dorgqr_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orgqr_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_orgqr(c, ung)
{
    {
        LAPACK_RETURN_CHECK( cungqr_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orgqr_body(c)
    }
}
LAPACK_orgqr(z, ung)
{
    {
        LAPACK_RETURN_CHECK( zungqr_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orgqr_body(z)
    }
}
#endif

#define LAPACK_org2r(prefix, name)                                      \
  int F77_ ## prefix ## name ## 2r( int* m,                                  \
                                    int* n,                             \
                                    int* k,                             \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                                    int* ldim_A,                        \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, \
                                    int* info)

LAPACK_org2r(s, org)
{
    {
        LAPACK_RETURN_CHECK( sorg2r_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_orgqr_body(s)
    }
}
LAPACK_org2r(d, org)
{
    {
        LAPACK_RETURN_CHECK( dorg2r_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_orgqr_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_org2r(c, ung)
{
    {
        LAPACK_RETURN_CHECK( cung2r_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_orgqr_body(c)
    }
}
LAPACK_org2r(z, ung)
{
    {
        LAPACK_RETURN_CHECK( zung2r_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_orgqr_body(z)
    }
}
#endif

#endif
