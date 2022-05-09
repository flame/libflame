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
   GELQF computes an LQ factorization of a M-by-N matrix A: A = L * Q.
*/

#define LAPACK_gelqf(prefix)                                            \
  int F77_ ## prefix ## gelqf( integer* m,                                  \
                               integer* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, integer* lwork, \
                               integer* info )

#define LAPACK_gelqf_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Obj      A, t, T;                                         \
  integer          min_m_n  = min( *m, *n );                        \
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


LAPACK_gelqf(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgelqf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( sgelqf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ), fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_gelqf_body(s)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_gelqf(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgelqf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dgelqf_check(m, n,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w, lwork,
                                              info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gelqf_body(d)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gelqf(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgelqf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cgelqf_check(m, n,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w, lwork,
                                              info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gelqf_body(c)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_gelqf(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgelqf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zgelqf_check(m, n,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w, lwork,
                                              info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gelqf_body(z)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#define LAPACK_gelq2(prefix)                                            \
  int F77_ ## prefix ## gelq2( integer* m,                                  \
                               integer* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,   \
                               integer* info )
LAPACK_gelq2(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgelq2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( sgelq2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ),fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gelqf_body(s)
      /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_gelq2(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgelq2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dgelq2_check(m, n,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w,
                                              info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gelqf_body(d)
      /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gelq2(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgelq2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cgelq2_check(m, n,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w,
                                              info),
                                 fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gelqf_body(c)
      /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_gelq2(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgelq2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( zgelq2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ),fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gelqf_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error=0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif


#endif
