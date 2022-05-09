/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_YET_LAPACK2FLAME

#include "FLA_lapack2flame_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  GEHRD reduces a general matrix A to upper Hessenberg form H by
  an unitary similarity transformation: Q**H * A * Q = H .
*/

#define LAPACK_gehrd(prefix)                                            \
  int F77_ ## prefix ## gehrd( integer* m,                                  \
                               integer* ilo,                                \
                               integer* ihi,                                \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, integer* lwork, \
                               integer* info )

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


LAPACK_gehrd(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgehrd inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS "", *m, *ilo, *ihi, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( sgehrd_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ), fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gehrd_body(s)
        fla_error=0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_gehrd(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgehrd inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS "", *m, *ilo, *ihi, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( dgehrd_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gehrd_body(d)
             /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gehrd(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgehrd inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS "", *m, *ilo, *ihi, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cgehrd_check(m,
                                              ilo, ihi,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w, lwork,
                                              info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gehrd_body(c)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_gehrd(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgehrd inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS "", *m, *ilo, *ihi, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zgehrd_check(m,
                                              ilo, ihi,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w, lwork,
                                              info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gehrd_body(z)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#define LAPACK_gehd2(prefix)                                            \
  int F77_ ## prefix ## gehd2( integer* m,                                  \
                               integer* ilo,                                \
                               integer* ihi,                                \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,   \
                               int* info )

LAPACK_gehd2(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgehd2 inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS "", *m, *ilo, *ihi, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( sgehd2_check( m,
                                           ilo, ihi,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gehrd_body(s)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_gehd2(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgehd2 inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS "", *m, *ilo, *ihi, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dgehd2_check(m,
                                              ilo, ihi,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w,
                                              info),fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gehrd_body(d)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gehd2(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgehd2 inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS "", *m, *ilo, *ihi, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cgehd2_check(m,
                                              ilo, ihi,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w,
                                              info),fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gehrd_body(c)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_gehd2(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgehd2 inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", lda %" FLA_IS "", *m, *ilo, *ihi, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zgehd2_check(m,
                                              ilo, ihi,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w,
                                              info),fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gehrd_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#endif
