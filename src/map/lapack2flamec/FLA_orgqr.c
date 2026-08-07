/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_prototypes.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_util_defs.h"

/*
  SORGQR generates an M-by-N real matrix Q with orthonormal columns,
  which is defined as the first N columns of a product of K elementary
  reflectors of order M

  Q = H(1) H(2) . . . H(k)

  as returned by SGEQRF.
*/

extern int lapack_dorgqr(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a, aocl_int64_t *lda,
                         doublereal *tau, doublereal *work, aocl_int64_t *lwork, aocl_int64_t *info);
extern void sorgqr_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, real *a, aocl_int64_t *lda, real *tau,
                       real *work, aocl_int64_t *lwork, aocl_int64_t *info);
extern void dorg2r_fla(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, doublereal *a, aocl_int64_t *lda,
                       doublereal *tau, doublereal *work, aocl_int64_t *info);


/** Generated wrapper function */
void sorgqr_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sorgqr(m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sorgqr(&m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dorgqr_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dorgqr(m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dorgqr(&m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sorg2r_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sorg2r(m, n, k, buff_A, ldim_A, buff_t, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sorg2r(&m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dorg2r_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dorg2r(m, n, k, buff_A, ldim_A, buff_t, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dorg2r(&m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

#define LAPACK_orgqr(prefix, name)                                                  \
    void aocl_lapack_##prefix##name##qr(                                                    \
        aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, \
        aocl_int64_t * ldim_A, PREFIX2LAPACK_TYPEDEF(prefix) * buff_t,                   \
        PREFIX2LAPACK_TYPEDEF(prefix) * buff_w, aocl_int64_t * lwork, aocl_int64_t * info)

#define LAPACK_orgqr_body(prefix)                              \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);     \
    FLA_Obj A, AL, AR, t, T;                                   \
    FLA_Error init_result;                                     \
                                                               \
    FLA_Init_safe(&init_result);                               \
                                                               \
    FLA_Obj_create_without_buffer(datatype, *m, *n, &A);       \
    FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);             \
                                                               \
    if(*k > 0 && !(PREFIX2FLAME_IS_ZERO(prefix, buff_t)))      \
    {                                                          \
        FLA_Obj_create_without_buffer(datatype, *k, 1, &t);    \
        FLA_Obj_attach_buffer(buff_t, 1, *k, &t);              \
        PREFIX2FLAME_INVERT_TAU(prefix, t);                    \
                                                               \
        FLA_Part_1x2(A, &AL, &AR, *k, FLA_LEFT);               \
        FLA_QR_UT_create_T(AL, &T);                            \
        FLA_Set(FLA_ZERO, T);                                  \
        FLA_Accum_T_UT(FLA_FORWARD, FLA_COLUMNWISE, AL, t, T); \
        FLA_QR_UT_form_Q(AL, T, A);                            \
                                                               \
        PREFIX2FLAME_INVERT_TAU(prefix, t);                    \
        FLA_Obj_free_without_buffer(&t);                       \
        FLA_Obj_free(&T);                                      \
    }                                                          \
    else                                                       \
    {                                                          \
        FLA_Set_to_identity(A);                                \
    }                                                          \
    FLA_Obj_free_without_buffer(&A);                           \
    FLA_Finalize_safe(init_result);                            \
                                                               \
    *info = 0;

LAPACK_orgqr(s, org)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sorgqr inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "",
                      *m, *n, *k, *ldim_A);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(sorgqr_check(m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orgqr_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
#else
    {
        sorgqr_fla(m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}
LAPACK_orgqr(d, org)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dorgqr inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "",
                      *m, *n, *k, *ldim_A);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(dorgqr_check(m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orgqr_body(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
#else
    {
        lapack_dorgqr(m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_orgqr(c, ung)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cungqr inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "",
                      *m, *n, *k, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cungqr_check(m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orgqr_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_orgqr(z, ung)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zungqr inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "",
                      *m, *n, *k, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zungqr_check(m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orgqr_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#define LAPACK_org2r(prefix, name)                                                        \
    void aocl_lapack_##prefix##name##2r(aocl_int64_t * m, aocl_int64_t * n, aocl_int64_t * k,                    \
                                PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, aocl_int64_t * ldim_A, \
                                PREFIX2LAPACK_TYPEDEF(prefix) * buff_t,                   \
                                PREFIX2LAPACK_TYPEDEF(prefix) * buff_w, aocl_int64_t * info)

LAPACK_org2r(s, org)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sorg2r inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "",
                      *m, *n, *k, *ldim_A);
    extern void sorg2r_fla(aocl_int64_t * m, aocl_int64_t * n, aocl_int64_t * k, real * a, aocl_int64_t * lda,
                           real * tau, real * work, aocl_int64_t * info);

#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(sorg2r_check(m, n, k, buff_A, ldim_A, buff_t, buff_w, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orgqr_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
#else
    {
        sorg2r_fla(m, n, k, buff_A, ldim_A, buff_t, buff_w, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}

LAPACK_org2r(d, org)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dorg2r inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "",
                      *m, *n, *k, *ldim_A);

#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(dorg2r_check(m, n, k, buff_A, ldim_A, buff_t, buff_w, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orgqr_body(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
#else
    {
        dorg2r_fla(m, n, k, buff_A, ldim_A, buff_t, buff_w, info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_org2r(c, ung)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cung2r inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "",
                      *m, *n, *k, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cung2r_check(m, n, k, buff_A, ldim_A, buff_t, buff_w, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orgqr_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_org2r(z, ung)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zung2r inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "",
                      *m, *n, *k, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zung2r_check(m, n, k, buff_A, ldim_A, buff_t, buff_w, info),
                                 fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_orgqr_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
#endif

#endif
