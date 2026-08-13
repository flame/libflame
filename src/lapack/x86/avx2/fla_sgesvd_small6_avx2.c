/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_sgesvd_small6_avx2.c
 *  @brief SGESVD Small path (path 6)
 *  without the LQ Factorization.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"
#include "fla_lapack_x86_common.h"

#if FLA_ENABLE_AMD_OPT

/* SVD for small tall-matrices with QR factorization
 * already computed
 */
void fla_sgesvd_small6_avx2(aocl_int64_t wntus, aocl_int64_t wntvs, aocl_int64_t *m,
                            aocl_int64_t *n, real *a, aocl_int64_t *lda, real *qr,
                            aocl_int64_t *ldqr, real *s, real *u, aocl_int64_t *ldu, real *vt,
                            aocl_int64_t *ldvt, real *work, aocl_int64_t *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_SSMALL();

    aocl_int64_t ie;
    aocl_int64_t itau, itauq, itaup;
    aocl_int64_t rlen, knt;
    aocl_int64_t ni;
    aocl_int64_t tn;
    aocl_int64_t ncvt, nru;
    aocl_int64_t *ldau;
    aocl_int64_t c__1 = 1;

    real *tau, *tauq, *taup;
    real *e, *au;
    real stau, d__1;
    real dum[2];
    real c_zero = 0.f;
    real c_one = 1.f;

    /* indices for partitioning work buffer */
    ie = 1;
    itau = ie + *n;
    itauq = itau + *n;
    itaup = itauq + *n;

    /* parameter adjustments */
    a -= (1 + *lda);
    u -= (1 + *ldu);
    vt -= (1 + *ldvt);
    qr -= (1 + *ldqr);
    --s;
    --work;

    /* local variables initialization */
    v = &dum[0];
    ncvt = 0;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* QR Factorization */
    fla_sgeqrf_small(m, n, &a[1 + *lda], lda, &work[itau], &work[ie]);

    /* Upper Bidiagonalization */
    if(wntus)
    {
        nru = *n;
        au = u;
        ldau = ldu;
        /* Copy R to U */
        aocl_lapack_slacpy("U", n, n, &a[1 + *lda], lda, &au[1 + *ldau], ldau);
    }
    else
    {
        nru = 0;
        au = a;
        ldau = lda;
    }
    /* Set lower part of U to zero */
    tn = *n - 1;
    aocl_lapack_slaset("L", &tn, &tn, &c_zero, &c_zero, &au[2 + *ldau], ldau);

    FLA_BIDIAGONALIZE_SSMALL(*n, *n, au, ldau, tauq, taup, s, e);

    /* Form Vt' in vt from HH vectors in U (right bi-diagonalizing Q) */
    if(wntvs)
    {
        ncvt = *n;
        FLA_GESVD_ZERO_MAT(*n, *n, vt, ldvt);
        FLA_LARF_VTAPPLY_SMALL_SQR(n, au, ldau, taup, vt, ldvt, c_one);
    }

    /* Form U' in U (left bi-diagonalizing Q) */
    if(wntus)
    {
        FLA_LARF_UAPPLY_SMALL_SQR(n, au, ldau, tauq, u, ldu, taup, c_one);
    }

    /* Compute SVD for bi-diagonal matrix
     * (sbdsqr with no lwork)
     * */
    if(*n == 2)
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
            fla_srot_avx2(&nru, &u[1 + *ldu], &c__1, &u[1 + 2 * *ldu], &c__1, &cosl, &sinl);
        }

        /* Normalize singular values and scale corresponding vectors for 2x2 case */
        FLA_NORMALIZE_SINGULAR_VALUE_AND_VECTORS_2X2(1, wntus);
        FLA_NORMALIZE_SINGULAR_VALUE_AND_VECTORS_2X2(2, wntus);
    }
    else
    {
        /* Compute Singular Values and Vectors */
        lapack_sbdsqr_small("U", n, &ncvt, &nru, &s[1], &e[1], &vt[1 + *ldvt], ldvt, &u[1 + *ldu],
                            ldu, info);
    }

    /* Compute U by updating U' by applying from the left the Q from QR */
    if(wntus)
    {
        tau = &work[itau - 1];
        FLA_GESVD_UAPPLY_QR_SSMALL(m, n, u, ldu, qr, ldqr, tau);
    }

    return;
}
#endif
