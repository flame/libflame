/******************************************************************************
 * Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dgesvd_xx_small10_avx2.c
 *  @brief DGESVD Small path (Path 10)
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

void fla_dgesvd_xx_small10_avx2(aocl_int64_t wntu, aocl_int64_t wntv, aocl_int64_t *m,
                                aocl_int64_t *n, aocl_int64_t *ncu, doublereal *a,
                                aocl_int64_t *lda, doublereal *s, doublereal *u, aocl_int64_t *ldu,
                                doublereal *vt, aocl_int64_t *ldvt, doublereal *work,
                                aocl_int64_t *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    doublereal d__1;
    doublereal *tau, *tauq, *taup;
    doublereal *e;
    doublereal stau;
    doublereal c_one = 1.;
    doublereal cosu = 0., sinu = 0.;

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
        FLA_BIDIAG_2X2_GIVENS_DSMALL(a, lda, s, e, cosu, sinu);
    }
    else
    {
        FLA_BIDIAGONALIZE_DSMALL(*m, *n, a, lda, tauq, taup, s, e);

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
        FLA_GESVD_LASV2_2X2_DSMALL(s, e, sigmn, sigmx, sinr, cosr, sinl, cosl);
        /* Compute singular vectors, if desired */
        if(ncvt > 0)
        {
            FLA_COMPUTE_VT_D2X2(vt, ldvt, sigmx, sigmn, cosr, sinr);
        }
        if(nru > 0)
        {
            FLA_GESVD_U_FROM_2X2_GIVENS_DSMALL(u, ldu, cosl, sinl, cosu, sinu);
        }

        /* Normalize singular values and scale corresponding vectors for 2x2 case */
        FLA_NORMALIZE_SINGULAR_VALUE_AND_VECTORS_2X2(1, wntu);
        FLA_NORMALIZE_SINGULAR_VALUE_AND_VECTORS_2X2(2, wntu);
    }
    else
    {
        /* Compute Singular Values and Vectors */
        lapack_dbdsqr_small("U", n, &ncvt, &nru, &s[1], &e[1], &vt[1 + *ldvt], ldvt,
                            &u[1 + *ldu], ldu, info);
    }
    return;
}
#endif
