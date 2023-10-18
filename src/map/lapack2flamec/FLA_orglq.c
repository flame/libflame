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

extern int sorglq_fla(integer* m, integer* n, integer* k, real* a, integer* lda, real* tau, real* work, integer* lwork, integer* info);
extern int dorglq_fla(integer* m, integer* n, integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);
/*
  SORGLQ generates an M-by-N real matrix Q with orthonormal rows,
  which is defined as the first M rows of a product of K elementary
  reflectors of order N

  Q = H(k) . . . H(2) H(1)

  as returned by SGELQF.
*/

#define LAPACK_orglq(prefix, name)                                      \
  int F77_ ## prefix ## name ## lq( integer* m,                             \
                                    integer* n,                             \
                                    integer* k,                             \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                                    integer* ldim_A,                        \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, \
                                    integer* lwork,                         \
                                    integer* info)

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



LAPACK_orglq(s, org)
{

    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sorglq inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "", *m, *n, *k, *ldim_A);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( sorglq_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_orglq_body(s)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        sorglq_fla(m, n, k,
            buff_A, ldim_A,
            buff_t,
            buff_w, lwork,
            info);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}
LAPACK_orglq(d, org)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dorglq inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "", *m, *n, *k, *ldim_A);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( dorglq_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orglq_body(d)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        dorglq_fla(m, n, k,
            buff_A, ldim_A,
            buff_t,
            buff_w, lwork,
            info);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_orglq(c, ung)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cunglq inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "", *m, *n, *k, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( cunglq_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orglq_body(c)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_orglq(z, ung)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zunglq inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "", *m, *n, *k, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( zunglq_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orglq_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#define LAPACK_orgl2(prefix, name)                                      \
  int F77_ ## prefix ## name ## l2( integer* m,                                  \
                                    integer* n,                             \
                                    integer* k,                             \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                                    integer* ldim_A,                        \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, \
                                    integer* info)

LAPACK_orgl2(s, org)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sorgl2 inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "", *m, *n, *k, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( sorgl2_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orglq_body(s)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_orgl2(d, org)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dorgl2 inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "", *m, *n, *k, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( dorgl2_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orglq_body(d)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_orgl2(c, ung)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cungl2 inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "", *m, *n, *k, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( cungl2_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orglq_body(c)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_orgl2(z, ung)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zungl2 inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "", *m, *n, *k, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( zungl2_check( m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orglq_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#endif
