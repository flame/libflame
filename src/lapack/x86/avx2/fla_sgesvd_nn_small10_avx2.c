/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_sgesvd_nn_small10_avx2.c
 *  @brief SGESVD Small path (Path 10)
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

void fla_sgesvd_xx_small10_avx2(aocl_int64_t wntu, aocl_int64_t wntv, aocl_int64_t *m,
                                aocl_int64_t *n, aocl_int64_t *ncu, real *a, aocl_int64_t *lda,
                                real *s, real *u, aocl_int64_t *ldu, real *vt, aocl_int64_t *ldvt,
                                real *work, aocl_int64_t *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_SSMALL();

    real d__1;
    real *tau, *tauq, *taup;
    real *e;
    real stau;
    real c_one = 1.f;
    real cosu = 0.f, sinu = 0.f;

    aocl_int64_t ncvt, nru;
    aocl_int64_t c__1 = 1;

    aocl_int64_t ie;
    aocl_int64_t itauq, itaup;
    aocl_int64_t rlen, knt;

    /* indices for partitioning work buffer */
    ie = 1;
    itauq = ie + *n;
    itaup = itauq + *n;

    /* parameter adjustments */
    a -= (1 + *lda);
    u -= (1 + *ldu);
    vt -= (1 + *ldvt);
    --s;
    --work;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* Upper Bidiagonalization */
    if(*m == 2 && *n == 2)
    {
        /* 2x2 matrix Bi-Diag using Givens */
        FLA_BIDIAG_2X2_GIVENS_SSMALL(a, lda, s, e, cosu, sinu);
    }
    else
    {
        FLA_BIDIAGONALIZE_SSMALL(*m, *n, a, lda, tauq, taup, s, e);

        /* Generate Qr (from bidiag) in vt from work[iu] (a here) */
        if(wntv)
        {
            FLA_GESVD_ZERO_MAT(*n, *n, vt, ldvt);
            FLA_LARF_VTAPPLY_SMALL_SQR(n, a, lda, taup, vt, ldvt, c_one);
        }
        /* Generate Ql (from bidiag) in u from a */
        if(wntu)
        {
            FLA_GESVD_FORM_U_RECT_SMALL(m, n, ncu, a, lda, tauq, u, ldu);
        }
    }

    /* Compute final Singular Values/Vectors */
    ncvt = 0;
    nru = 0;
    if(wntv)
    {
        ncvt = *n;
    }
    if(wntu)
    {
        nru = *m;
    }
    if(*m == 2 && *n == 2)
    {
        /* 2 by 2 block, handle separately */
        FLA_GESVD_LASV2_2X2_SSMALL(s, e, sigmn, sigmx, sinr, cosr, sinl, cosl);
        /* Compute singular vectors, if desired */
        if(ncvt > 0)
        {
            FLA_COMPUTE_VT_S2X2(vt, ldvt, sigmx, sigmn, cosr, sinr);
        }
        if(nru > 0)
        {
            FLA_GESVD_U_FROM_2X2_GIVENS_SSMALL(u, ldu, cosl, sinl, cosu, sinu);
        }

        /* Normalize singular values and scale corresponding vectors for 2x2 case */
        FLA_NORMALIZE_SINGULAR_VALUE_AND_VECTORS_2X2(1, wntu);
        FLA_NORMALIZE_SINGULAR_VALUE_AND_VECTORS_2X2(2, wntu);
    }
    else
    {
        /* Compute Singular Values and Vectors */
        lapack_sbdsqr_small("U", n, &ncvt, &nru, &s[1], &e[1], &vt[1 + *ldvt], ldvt, &u[1 + *ldu],
                            ldu, info);
    }
    return;
}
#endif
