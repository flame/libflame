/******************************************************************************
 * Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dgesvd_small6T_avx2_.c
 *  @brief DGESVD Small path (path 6T)
 *  without the LQ Factorization.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* SVD for small fat-matrices with LQ factorization
 * already computed
 */
void fla_dgesvd_small6T_avx2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                             doublereal *ql, aocl_int64_t *ldql, doublereal *s, doublereal *u,
                             aocl_int64_t *ldu, doublereal *vt, aocl_int64_t *ldvt,
                             doublereal *work, aocl_int64_t *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    aocl_int64_t iu, ie;
    aocl_int64_t itau, itauq, itaup;
    aocl_int64_t rlen, knt;
    aocl_int64_t c__1 = 1;

    doublereal *tau, *tauq, *taup;
    doublereal *e, *vtau, *avt;
    doublereal stau, d__1;
    doublereal c_one = 1.;

    /* indices for partitioning work buffer */
    iu = 1;
    itau = iu + *lda * *m;
    ie = itau + *m;
    itauq = ie + *m;
    itaup = itauq + *m;

    /* parameter adjustments */
    a -= (1 + *lda);
    u -= (1 + *ldu);
    vt -= (1 + *ldvt);
    ql -= (1 + *ldql);
    --s;
    --work;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* Upper Bidiagonalization */
    FLA_BIDIAGONALIZE_DSMALL(*m, *m, a, lda, tauq, taup, s, e);

    /* Generate Qr (from bidiag) in vt from work[iu] (a here) */
    FLA_GESVD_ZERO_MAT(*m, *n, vt, ldvt);
    FLA_LARF_VTAPPLY_SMALL_SQR(m, a, lda, taup, vt, ldvt, c_one);

    /* Generate Ql (from bidiag) in u from a */
    FLA_LARF_UAPPLY_SMALL_SQR(m, a, lda, tauq, u, ldu, taup, c_one);

    /* Compute Singular Values and Vetcors */
    lapack_dbdsqr_small("U", m, m, m, &s[1], &e[1], &vt[1 + *ldvt], ldvt, &u[1 + *ldu], ldu, info);

    /* Apply HH from LQ factorization (ql) on vt from right */

    tau = &work[itau - 1];
    vtau = tau + *m;
    avt = vtau + *n;
    FLA_GESVD_VTAPPLY_LQ_SMALL(m, n, vt, ldvt, ql, ldql, tau, vtau, avt);

    return;
}
#endif
