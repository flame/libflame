/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  ORMQR overwrites the general real M-by-N matrix C with

  SIDE = 'L' SIDE = 'R'
  TRANS = 'N': Q * C C * Q
  TRANS = 'T': Q**T * C C * Q**T

  where Q is a real orthogonal matrix defined as the product of k
  elementary reflectors

  Q = H(1) H(2) . . . H(k)

  as returned by GEQRF. Q is of order M if SIDE = 'L' and of order N
  if SIDE = 'R'.
*/

extern int dormqr_fla(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal * c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
extern int sormqr_fla(char *side, char *trans, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real * c__, integer *ldc, real *work, integer *lwork, integer *info);

#define LAPACK_ormqr(prefix, name)                                      \
  int F77_ ## prefix ## name ## qr( char* side,                              \
                                    char* trans,                        \
                                    integer* m,                             \
                                    integer* n,                             \
                                    integer* k,                             \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, integer* ldim_B, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, integer* lwork, \
                                    integer* info )

#define LAPACK_ormqr_body(prefix)                                       \
  FLA_Datatype datatype   = PREFIX2FLAME_DATATYPE(prefix);              \
  FLA_Side     side_fla;                                                \
  FLA_Trans    trans_fla;                                               \
  FLA_Error    init_result;                                             \
  dim_t        mq /*, nq */;                                            \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Param_map_netlib_to_flame_side( side, &side_fla );                \
  FLA_Param_map_netlib_to_flame_trans( trans, &trans_fla );             \
                                                                        \
  if      ( side_fla == FLA_LEFT)     { mq = *m; /* nq = *n; */ }       \
  else /* ( side_fla == FLA_RIGHT) */ { mq = *n; /* nq = *m; */ }       \
                                                                        \
  if ( *k > 0 && !( PREFIX2FLAME_IS_ZERO(prefix, buff_t) ) )            \
    {                                                                   \
      FLA_Obj      A, t, B, T, W;                                       \
                                                                        \
      FLA_Obj_create_without_buffer( datatype, mq, *k, &A );            \
      FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                  \
                                                                        \
      FLA_Obj_create_without_buffer( datatype, *m, *n, &B );            \
      FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );                  \
                                                                        \
      FLA_Obj_create_without_buffer( datatype, *k, 1, &t );             \
      FLA_Obj_attach_buffer( buff_t, 1, *k, &t );                       \
      PREFIX2FLAME_INVERT_TAU(prefix,t);                                \
                                                                        \
      FLA_QR_UT_create_T( A, &T );                                      \
      FLA_Set( FLA_ZERO, T );                                           \
      FLA_Apply_Q_UT_create_workspace_side( side_fla, T, B, &W);        \
      FLA_Accum_T_UT( FLA_FORWARD, FLA_COLUMNWISE, A, t, T );           \
      FLA_Apply_Q_UT( side_fla, trans_fla, FLA_FORWARD, FLA_COLUMNWISE, \
                      A, T, W, B );                                     \
      FLA_Obj_free( &W );                                               \
      FLA_Obj_free( &T );                                               \
                                                                        \
      PREFIX2FLAME_INVERT_TAU(prefix,t);                                \
      FLA_Obj_free_without_buffer( &t );                                \
      FLA_Obj_free_without_buffer( &B );                                \
      FLA_Obj_free_without_buffer( &A );                                \
    }                                                                   \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \


LAPACK_ormqr(s, orm)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sormqr inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "", *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
#if !FLA_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( sormqr_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_ormqr_body(s)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        sormqr_fla(side, trans, m, n, k,
                   buff_A, ldim_A,
                   buff_t,
                   buff_B, ldim_B,
                   buff_w, lwork, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}
LAPACK_ormqr(d, orm)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dormqr inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "", *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
#if !FLA_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( dormqr_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormqr_body(d)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        dormqr_fla(side, trans, m, n, k,
                   buff_A, ldim_A,
                   buff_t,
                   buff_B, ldim_B,
                   buff_w, lwork, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_ormqr(c, unm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cunmqr inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "", *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1( cunmqr_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormqr_body(c)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_ormqr(z, unm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zunmqr inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "", *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1( zunmqr_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormqr_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#define LAPACK_orm2r(prefix, name)                                      \
  int F77_ ## prefix ## name ## 2r( char* side,                              \
                                    char* trans,                        \
                                    integer* m,                             \
                                    integer* n,                             \
                                    integer* k,                             \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, integer* ldim_B, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, \
                                    integer* info )

LAPACK_orm2r(s, orm)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sorm2r inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "", *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    extern int sorm2r_fla(char *side, char *trans, integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, real *work, integer *info);

#if !FLA_AMD_OPT 
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( sorm2r_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormqr_body(s)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        sorm2r_fla( side, trans,
                    m, n, k,
                    buff_A, ldim_A,
                    buff_t,
                    buff_B, ldim_B,
                    buff_w,
                    info );
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}

LAPACK_orm2r(d, orm)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dorm2r inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "", *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    extern int dorm2r_fla(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal * c__, integer *ldc, doublereal *work, integer *info);

#if !FLA_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( dorm2r_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormqr_body(d)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        dorm2r_fla( side, trans,
                    m, n, k,
                    buff_A, ldim_A,
                    buff_t,
                    buff_B, ldim_B,
                    buff_w,
                    info );
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;

    }
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_orm2r(c, unm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cunm2r inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "", *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1( cunm2r_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormqr_body(c)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_orm2r(z, unm)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zunm2r inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS ", ldc %" FLA_IS "", *side, *trans, *m, *n, *k, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1( zunm2r_check( side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_B, ldim_B,
                                           buff_w,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_ormqr_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#endif
