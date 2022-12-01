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
  SORGLQ generates an M-by-N real matrix Q with orthonormal rows,
  which is defined as the first M rows of a product of K elementary
  reflectors of order N

  Q = H(k) . . . H(2) H(1)

  as returned by SGELQF.
*/

#define LAPACK_orglq(prefix, name)                                      \
  int F77_ ## prefix ## name ## lq( int* m,                             \
                                    int* n,                             \
                                    int* k,                             \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                                    int* ldim_A,                        \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, \
                                    int* lwork,                         \
                                    int* info)

#define LAPACK_orglq_body(prefix)                                       \
  FLA_Datatype datatype   = PREFIX2FLAME_DATATYPE(prefix);              \
  FLA_Obj      A, AT, AB, t, T;                                         \
  FLA_Error    init_result;                                             \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  if ( *k > 0 && !( PREFIX2FLAME_IS_ZERO(prefix, buff_t) ) )            \
    {                                                                   \
      FLA_Obj_create_without_buffer( datatype, *k, 1, &t );             \
      FLA_Obj_attach_buffer( buff_t, 1, *k, &t );                       \
      PREFIX2FLAME_INVERT_TAU(prefix,t);                                \
                                                                        \
      FLA_Part_2x1( A, &AT,                                             \
                    &AB, *k, FLA_TOP );                                 \
      FLA_LQ_UT_create_T( AT, &T );                                     \
      FLA_Set( FLA_ZERO, T );                                           \
      FLA_Accum_T_UT( FLA_FORWARD, FLA_ROWWISE, AT, t, T );             \
      FLA_LQ_UT_form_Q( AT, T, A );                                     \
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


LAPACK_orglq(s, org)
{
    {
        LAPACK_RETURN_CHECK( sorglq_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orglq_body(s)
    }
}
LAPACK_orglq(d, org)
{
    {
        LAPACK_RETURN_CHECK( dorglq_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orglq_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_orglq(c, ung)
{
    {
        LAPACK_RETURN_CHECK( cunglq_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orglq_body(c)
    }
}
LAPACK_orglq(z, ung)
{
    {
        LAPACK_RETURN_CHECK( zunglq_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orglq_body(z)
    }
}
#endif

#define LAPACK_orgl2(prefix, name)                                      \
  int F77_ ## prefix ## name ## l2( int* m,                                  \
                                    int* n,                             \
                                    int* k,                             \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                                    int* ldim_A,                        \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, \
                                    int* info)

LAPACK_orgl2(s, org)
{
    {
        LAPACK_RETURN_CHECK( sorgl2_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_orglq_body(s)
    }
}
LAPACK_orgl2(d, org)
{
    {
        LAPACK_RETURN_CHECK( dorgl2_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_orglq_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_orgl2(c, ung)
{
    {
        LAPACK_RETURN_CHECK( cungl2_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_orglq_body(c)
    }
}
LAPACK_orgl2(z, ung)
{
    {
        LAPACK_RETURN_CHECK( zungl2_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_orglq_body(z)
    }
}
#endif

#endif
