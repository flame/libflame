/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
/** Generated wrapper function */
void sgebrd_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_d, real *buff_e, real *buff_tu, real *buff_tv, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgebrd(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu, buff_tv, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgebrd(&m_64, &n_64, buff_A, &ldim_A_64, buff_d, buff_e, buff_tu, buff_tv, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgebrd_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_d, doublereal *buff_e, doublereal *buff_tu, doublereal *buff_tv, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgebrd(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu, buff_tv, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgebrd(&m_64, &n_64, buff_A, &ldim_A_64, buff_d, buff_e, buff_tu, buff_tv, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgebd2_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_d, real *buff_e, real *buff_tu, real *buff_tv, real *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgebd2(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu, buff_tv, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgebd2(&m_64, &n_64, buff_A, &ldim_A_64, buff_d, buff_e, buff_tu, buff_tv, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgebd2_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_d, doublereal *buff_e, doublereal *buff_tu, doublereal *buff_tv, doublereal *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgebd2(m, n, buff_A, ldim_A, buff_d, buff_e, buff_tu, buff_tv, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgebd2(&m_64, &n_64, buff_A, &ldim_A_64, buff_d, buff_e, buff_tu, buff_tv, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgelqf_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgelqf(m, n, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgelqf(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgelqf_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgelqf(m, n, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgelqf(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgelq2_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgelq2(m, n, buff_A, ldim_A, buff_t, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgelq2(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgelq2_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgelq2(m, n, buff_A, ldim_A, buff_t, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgelq2(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgeqp3_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, real *buff_t, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgeqp3(m, n, buff_A, ldim_A, buff_p, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgeqp3(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgeqp3_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, doublereal *buff_t, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgeqp3(m, n, buff_A, ldim_A, buff_p, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgeqp3(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgeqrf_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgeqrf(m, n, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgeqrf(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgeqrf_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgeqrf(m, n, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgeqrf(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgeqr2_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgeqr2(m, n, buff_A, ldim_A, buff_t, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgeqr2(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgeqr2_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgeqr2(m, n, buff_A, ldim_A, buff_t, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgeqr2(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgeqrfp_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgeqrfp(m, n, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgeqrfp(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgeqrfp_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgeqrfp(m, n, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgeqrfp(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgeqr2p_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgeqr2p(m, n, buff_A, ldim_A, buff_t, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgeqr2p(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgeqr2p_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgeqr2p(m, n, buff_A, ldim_A, buff_t, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgeqr2p(&m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgesdd_(char *jobz, aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_s, real *buff_U, aocl_int_t *ldim_U, real *buff_Vh, aocl_int_t *ldim_Vh, real *buff_w, aocl_int_t *lwork, aocl_int_t *buff_i, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgesdd(jobz, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U, buff_Vh, ldim_Vh, buff_w, lwork, buff_i, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_U_64 = *ldim_U;
    aocl_int64_t ldim_Vh_64 = *ldim_Vh;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgesdd(jobz, &m_64, &n_64, buff_A, &ldim_A_64, buff_s, buff_U, &ldim_U_64, buff_Vh, &ldim_Vh_64, buff_w, &lwork_64, buff_i, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgesdd_(char *jobz, aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_s, doublereal *buff_U, aocl_int_t *ldim_U, doublereal *buff_Vh, aocl_int_t *ldim_Vh, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *buff_i, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgesdd(jobz, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U, buff_Vh, ldim_Vh, buff_w, lwork, buff_i, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_U_64 = *ldim_U;
    aocl_int64_t ldim_Vh_64 = *ldim_Vh;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgesdd(jobz, &m_64, &n_64, buff_A, &ldim_A_64, buff_s, buff_U, &ldim_U_64, buff_Vh, &ldim_Vh_64, buff_w, &lwork_64, buff_i, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgesvd_(char *jobu, char *jobv, aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_s, real *buff_U, aocl_int_t *ldim_U, real *buff_Vh, aocl_int_t *ldim_Vh, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgesvd(jobu, jobv, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U, buff_Vh, ldim_Vh, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_U_64 = *ldim_U;
    aocl_int64_t ldim_Vh_64 = *ldim_Vh;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgesvd(jobu, jobv, &m_64, &n_64, buff_A, &ldim_A_64, buff_s, buff_U, &ldim_U_64, buff_Vh, &ldim_Vh_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgesvd_(char *jobu, char *jobv, aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_s, doublereal *buff_U, aocl_int_t *ldim_U, doublereal *buff_Vh, aocl_int_t *ldim_Vh, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgesvd(jobu, jobv, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U, buff_Vh, ldim_Vh, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_U_64 = *ldim_U;
    aocl_int64_t ldim_Vh_64 = *ldim_Vh;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgesvd(jobu, jobv, &m_64, &n_64, buff_A, &ldim_A_64, buff_s, buff_U, &ldim_U_64, buff_Vh, &ldim_Vh_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgetrf_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgetrf(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgetrf(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgetrf_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgetrf(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgetrf(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void cgetrf_(aocl_int_t *m, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgetrf(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgetrf(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zgetrf_(aocl_int_t *m, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zgetrf(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zgetrf(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgetf2_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgetf2(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgetf2(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgetf2_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgetf2(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgetf2(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void cgetf2_(aocl_int_t *m, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgetf2(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgetf2(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zgetf2_(aocl_int_t *m, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zgetf2(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zgetf2(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void ssygst_(aocl_int_t *itype, char *uplo, aocl_int_t *m, real *buff_A, aocl_int_t *ldim_A, real *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ssygst(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ssygst(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dsygst_(aocl_int_t *itype, char *uplo, aocl_int_t *m, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dsygst(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dsygst(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void chegst_(aocl_int_t *itype, char *uplo, aocl_int_t *m, scomplex *buff_A, aocl_int_t *ldim_A, scomplex *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chegst(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chegst(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zhegst_(aocl_int_t *itype, char *uplo, aocl_int_t *m, dcomplex *buff_A, aocl_int_t *ldim_A, dcomplex *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zhegst(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zhegst(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void ssygs2_(aocl_int_t *itype, char *uplo, aocl_int_t *m, real *buff_A, aocl_int_t *ldim_A, real *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ssygs2(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ssygs2(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dsygs2_(aocl_int_t *itype, char *uplo, aocl_int_t *m, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dsygs2(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dsygs2(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void chegs2_(aocl_int_t *itype, char *uplo, aocl_int_t *m, scomplex *buff_A, aocl_int_t *ldim_A, scomplex *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chegs2(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chegs2(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zhegs2_(aocl_int_t *itype, char *uplo, aocl_int_t *m, dcomplex *buff_A, aocl_int_t *ldim_A, dcomplex *buff_B, aocl_int_t *ldim_B, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zhegs2(itype, uplo, m, buff_A, ldim_A, buff_B, ldim_B, info);
#else
    aocl_int64_t itype_64 = *itype;
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zhegs2(&itype_64, uplo, &m_64, buff_A, &ldim_A_64, buff_B, &ldim_B_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void ssytrd_(char *uplo, aocl_int_t *m, real *buff_A, aocl_int_t *ldim_A, real *buff_d, real *buff_e, real *buff_t, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ssytrd(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ssytrd(uplo, &m_64, buff_A, &ldim_A_64, buff_d, buff_e, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dsytrd_(char *uplo, aocl_int_t *m, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_d, doublereal *buff_e, doublereal *buff_t, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dsytrd(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dsytrd(uplo, &m_64, buff_A, &ldim_A_64, buff_d, buff_e, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void ssytd2_(char *uplo, aocl_int_t *m, real *buff_A, aocl_int_t *ldim_A, real *buff_d, real *buff_e, real *buff_t, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ssytd2(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ssytd2(uplo, &m_64, buff_A, &ldim_A_64, buff_d, buff_e, buff_t, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dsytd2_(char *uplo, aocl_int_t *m, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_d, doublereal *buff_e, doublereal *buff_t, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dsytd2(uplo, m, buff_A, ldim_A, buff_d, buff_e, buff_t, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dsytd2(uplo, &m_64, buff_A, &ldim_A_64, buff_d, buff_e, buff_t, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void slauum_(char *uplo, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_slauum(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_slauum(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dlauum_(char *uplo, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dlauum(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dlauum(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void clauum_(char *uplo, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clauum(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_clauum(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zlauum_(char *uplo, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlauum(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zlauum(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void slauu2_(char *uplo, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_slauu2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_slauu2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dlauu2_(char *uplo, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dlauu2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dlauu2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void clauu2_(char *uplo, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clauu2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_clauu2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zlauu2_(char *uplo, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zlauu2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zlauu2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sorgbr_(char *vect, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sorgbr(vect, m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sorgbr(vect, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dorgbr_(char *vect, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dorgbr(vect, m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dorgbr(vect, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sorglq_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sorglq(m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sorglq(&m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dorglq_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dorglq(m, n, k, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dorglq(&m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sorgl2_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sorgl2(m, n, k, buff_A, ldim_A, buff_t, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sorgl2(&m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dorgl2_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dorgl2(m, n, k, buff_A, ldim_A, buff_t, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dorgl2(&m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

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

/** Generated wrapper function */
void sorgtr_(char *uplo, aocl_int_t *m, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sorgtr(uplo, m, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sorgtr(uplo, &m_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dorgtr_(char *uplo, aocl_int_t *m, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dorgtr(uplo, m, buff_A, ldim_A, buff_t, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dorgtr(uplo, &m_64, buff_A, &ldim_A_64, buff_t, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sormbr_(char *vect, char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_C, aocl_int_t *ldim_C, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sormbr(vect, side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_C, ldim_C, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_C_64 = *ldim_C;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sormbr(vect, side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_C, &ldim_C_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dormbr_(char *vect, char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_C, aocl_int_t *ldim_C, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dormbr(vect, side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_C, ldim_C, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_C_64 = *ldim_C;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dormbr(vect, side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_C, &ldim_C_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sormlq_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_B, aocl_int_t *ldim_B, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sormlq(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sormlq(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dormlq_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_B, aocl_int_t *ldim_B, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dormlq(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dormlq(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sorml2_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_B, aocl_int_t *ldim_B, real *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sorml2(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sorml2(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dorml2_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_B, aocl_int_t *ldim_B, doublereal *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dorml2(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dorml2(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sormqr_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_B, aocl_int_t *ldim_B, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sormqr(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sormqr(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dormqr_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_B, aocl_int_t *ldim_B, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dormqr(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dormqr(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sorm2r_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_B, aocl_int_t *ldim_B, real *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sorm2r(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sorm2r(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dorm2r_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_B, aocl_int_t *ldim_B, doublereal *buff_w, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dorm2r(side, trans, m, n, k, buff_A, ldim_A, buff_t, buff_B, ldim_B, buff_w, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_B_64 = *ldim_B;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dorm2r(side, trans, &m_64, &n_64, &k_64, buff_A, &ldim_A_64, buff_t, buff_B, &ldim_B_64, buff_w, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sormtr_(char *side, char *uplo, char *trans, aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, real *buff_t, real *buff_C, aocl_int_t *ldim_C, real *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sormtr(side, uplo, trans, m, n, buff_A, ldim_A, buff_t, buff_C, ldim_C, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_C_64 = *ldim_C;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sormtr(side, uplo, trans, &m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_C, &ldim_C_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dormtr_(char *side, char *uplo, char *trans, aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, doublereal *buff_t, doublereal *buff_C, aocl_int_t *ldim_C, doublereal *buff_w, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dormtr(side, uplo, trans, m, n, buff_A, ldim_A, buff_t, buff_C, ldim_C, buff_w, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t ldim_C_64 = *ldim_C;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dormtr(side, uplo, trans, &m_64, &n_64, buff_A, &ldim_A_64, buff_t, buff_C, &ldim_C_64, buff_w, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void spotrf_(char *uplo, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_spotrf(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_spotrf(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dpotrf_(char *uplo, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dpotrf(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dpotrf(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void cpotrf_(char *uplo, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cpotrf(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cpotrf(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zpotrf_(char *uplo, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zpotrf(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zpotrf(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void spotf2_(char *uplo, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_spotf2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_spotf2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dpotf2_(char *uplo, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dpotf2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dpotf2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void cpotf2_(char *uplo, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cpotf2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cpotf2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zpotf2_(char *uplo, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zpotf2(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zpotf2(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void spotri_(char *uplo, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_spotri(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_spotri(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dpotri_(char *uplo, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dpotri(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dpotri(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void cpotri_(char *uplo, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cpotri(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cpotri(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zpotri_(char *uplo, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zpotri(uplo, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zpotri(uplo, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sspffrt2_(real *ap, aocl_int_t *n, aocl_int_t *ncolm, real *work, real *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sspffrt2(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm = *ncolm;

    aocl_lapack_sspffrt2(ap, &n_64, &ncolm, work, work2);
#endif
}

/** Generated wrapper function */
void dspffrt2_(doublereal *ap, aocl_int_t *n, aocl_int_t *ncolm, doublereal *work, doublereal *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dspffrt2(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm = *ncolm;

    aocl_lapack_dspffrt2(ap, &n_64, &ncolm, work, work2);
#endif
}

/** Generated wrapper function */
void cspffrt2_(scomplex *ap, aocl_int_t *n, aocl_int_t *ncolm, scomplex *work, scomplex *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cspffrt2(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm = *ncolm;

    aocl_lapack_cspffrt2(ap, &n_64, &ncolm, work, work2);
#endif
}

/** Generated wrapper function */
void zspffrt2_(dcomplex *ap, aocl_int_t *n, aocl_int_t *ncolm, dcomplex *work, dcomplex *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zspffrt2(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm = *ncolm;

    aocl_lapack_zspffrt2(ap, &n_64, &ncolm, work, work2);
#endif
}

/** Generated wrapper function */
void sspffrtx_(real *ap, aocl_int_t *n, aocl_int_t *ncolm, real *work, real *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sspffrtx(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm = *ncolm;

    aocl_lapack_sspffrtx(ap, &n_64, &ncolm, work, work2);
#endif
}

/** Generated wrapper function */
void dspffrtx_(doublereal *ap, aocl_int_t *n, aocl_int_t *ncolm, doublereal *work, doublereal *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dspffrtx(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm = *ncolm;

    aocl_lapack_dspffrtx(ap, &n_64, &ncolm, work, work2);
#endif
}

/** Generated wrapper function */
void cspffrtx_(scomplex *ap, aocl_int_t *n, aocl_int_t *ncolm, scomplex *work, scomplex *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cspffrtx(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm = *ncolm;

    aocl_lapack_cspffrtx(ap, &n_64, &ncolm, work, work2);
#endif
}

/** Generated wrapper function */
void zspffrtx_(dcomplex *ap, aocl_int_t *n, aocl_int_t *ncolm, dcomplex *work, dcomplex *work2)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zspffrtx(ap, n, ncolm, work, work2);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ncolm = *ncolm;

    aocl_lapack_zspffrtx(ap, &n_64, &ncolm, work, work2);
#endif
}

/** Generated wrapper function */
void strtri_(char *uplo, char *diag, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_strtri(uplo, diag, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_strtri(uplo, diag, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dtrtri_(char *uplo, char *diag, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dtrtri(uplo, diag, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dtrtri(uplo, diag, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void ctrtri_(char *uplo, char *diag, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ctrtri(uplo, diag, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ctrtri(uplo, diag, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void ztrtri_(char *uplo, char *diag, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ztrtri(uplo, diag, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ztrtri(uplo, diag, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void strti2_(char *uplo, char *diag, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_strti2(uplo, diag, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_strti2(uplo, diag, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dtrti2_(char *uplo, char *diag, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dtrti2(uplo, diag, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dtrti2(uplo, diag, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void ctrti2_(char *uplo, char *diag, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ctrti2(uplo, diag, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ctrti2(uplo, diag, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void ztrti2_(char *uplo, char *diag, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_ztrti2(uplo, diag, n, buff_A, ldim_A, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_ztrti2(uplo, diag, &n_64, buff_A, &ldim_A_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

