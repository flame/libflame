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
  DORMLQ overwrites the general real M-by-N matrix C with
  SIDE = 'L' SIDE = 'R'
  TRANS = 'N': Q * C C * Q
  TRANS = 'T': Q**T * C C * Q**T

  where Q is a real orthogonal matrix defined as the product of k
  elementary reflectors

  Q = H(k) . . . H(2) H(1)

  as returned by (real)GELQF. Q is of order M if SIDE = 'L' and of order N
  if SIDE = 'R'.
*/

#define LAPACK_ormlq(prefix, name)                                      \
  int F77_ ## prefix ## name ## lq( char* side,                         \
                                    char* trans,                        \
                                    int* m,                             \
                                    int* n,                             \
                                    int* k,                             \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, int* ldim_B, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                                    int* info )

#define LAPACK_ormlq_body(prefix)                                       \
  FLA_Datatype datatype   = PREFIX2FLAME_DATATYPE(prefix);              \
  FLA_Side     side_fla;                                                \
  FLA_Trans    trans_fla;                                               \
  FLA_Error    init_result;                                             \
  dim_t        /*mq, */ nq;                                             \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Param_map_netlib_to_flame_side( side, &side_fla );                \
  FLA_Param_map_netlib_to_flame_trans( trans, &trans_fla );             \
                                                                        \
  if    ( side_fla == FLA_LEFT )      { /* mq = *n; */ nq = *m; }       \
  else /* side_fla == FLA_RIGHT ) */  { /* mq = *m; */ nq = *n; }       \
                                                                        \
  if ( *k > 0 && !( PREFIX2FLAME_IS_ZERO(prefix, buff_t) ) )            \
    {                                                                   \
      FLA_Obj      A, t, B, T, W;                                       \
      FLA_Obj_create_without_buffer( datatype, *k, nq, &A );            \
      FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                  \
                                                                        \
      FLA_Obj_create_without_buffer( datatype, *m, *n, &B );            \
      FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );                  \
                                                                        \
      FLA_Obj_create_without_buffer( datatype, *k, 1, &t );             \
      FLA_Obj_attach_buffer( buff_t, 1, *k, &t );                       \
      PREFIX2FLAME_INVERT_TAU(prefix,t);                                \
                                                                        \
      FLA_LQ_UT_create_T( A, &T );                                      \
      FLA_Set( FLA_ZERO, T );                                           \
      FLA_Apply_Q_UT_create_workspace_side( side_fla, T, B, &W);        \
                                                                        \
      FLA_Accum_T_UT( FLA_FORWARD, FLA_ROWWISE, A, t, T );              \
      FLA_Apply_Q_UT( side_fla, trans_fla, FLA_BACKWARD, FLA_ROWWISE,   \
                      A, T, W, B );                                     \
                                                                        \
      FLA_Obj_free( &W );                                               \
      FLA_Obj_free( &T );                                               \
                                                                        \
      PREFIX2FLAME_INVERT_TAU(prefix,t);                                \
      FLA_Obj_free_without_buffer( &t );                                \
      FLA_Obj_free_without_buffer( &B );                                \
      FLA_Obj_free_without_buffer( &A );                                \
    }                                                                   \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;


LAPACK_ormlq(s, orm)
{
    {
        LAPACK_RETURN_CHECK( sormlq_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_ormlq_body(s)
    }
}
LAPACK_ormlq(d, orm)
{
    {
        LAPACK_RETURN_CHECK( dormlq_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_ormlq_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_ormlq(c, unm)
{
    {
        LAPACK_RETURN_CHECK( cunmlq_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_ormlq_body(c)
    }
}
LAPACK_ormlq(z, unm)
{
    {
        LAPACK_RETURN_CHECK( zunmlq_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_ormlq_body(z)
    }
}
#endif

#define LAPACK_orml2(prefix, name)                                      \
  int F77_ ## prefix ## name ## l2( char* side,                         \
                                    char* trans,                        \
                                    int* m,                             \
                                    int* n,                             \
                                    int* k,                             \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, int* ldim_B, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, \
                                    int* info )

LAPACK_orml2(s, orm)
{
    {
        LAPACK_RETURN_CHECK( sorml2_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_ormlq_body(s)
    }
}
LAPACK_orml2(d, orm)
{
    {
        LAPACK_RETURN_CHECK( dorml2_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_ormlq_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_orml2(c, unm)
{
    {
        LAPACK_RETURN_CHECK( cunml2_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_ormlq_body(c)
    }
}
LAPACK_orml2(z, unm)
{
    {
        LAPACK_RETURN_CHECK( zunml2_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_ormlq_body(z)
    }
}
#endif

#endif
