/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_dgesvd_small6_avx2_.c
 *  @brief DGESVD Small path (path 6)
 *  without the LQ Factorization.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

double d_sign(doublereal *, doublereal *);

static integer c__0 = 0;
static integer c__1 = 1;

/* SVD for small fat-matrices with LQ factorization
 * already computed
 */
void fla_dgesvd_small6_avx2(integer *m, integer *n,
                            doublereal *a, integer *lda,
                            doublereal *qr, integer *ldqr,
                            doublereal *s,
                            doublereal *u, integer *ldu,
                            doublereal *vt, integer *ldvt,
                            doublereal *work,
                            integer *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    integer iu, ie, iwork;
    integer itau, itauq, itaup;
    integer i__1, rlen, knt;
    integer ni;

    doublereal *tau, *tauq, *taup;
    doublereal *e, *au;
    doublereal stau, d__1;
    doublereal dum[1];

    /* indices for partitioning work buffer */
    iu = 1;
    itau = iu + *n * *lda;
    ie = itau + *n;
    itauq = ie + *n;
    itaup = itauq + *n;
    iwork = itaup + *n;

    /* parameter adjustments */
    a -= (1 + *lda);
    u -= (1 + *ldu);
    vt -= (1 + *ldvt);
    qr -= (1 + *ldqr);
    v = &dum[-1];
    --s;
    --work;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* Upper Bidiagonalization */
    FLA_BIDIAGONALIZE_SMALL(*n, *n);

    for (i = 1; i <= *n; i++)
        for (j = 1; j <= *n; j++)
            vt[i + j * *ldvt] = 0.;
    /* Generate Qr (from bidiag) in vt from work[iu] (a here) */
    if (*n > 2)
    {
        /* iteration corresponding to (n - 2) HH[n-2] */
        stau = taup[*n - 2];
        d__1 = a[*n - 2 + *n * *lda];
        dtmp = - (stau * d__1); /* tau * v2 */

        vt[*n - 1 + (*n - 1) * *ldvt] = 1.0 - stau; /* 1 - tau */
        vt[*n + (*n - 1) * *ldvt] = dtmp; /* tau * v2 */
        vt[*n - 1 + *n * *ldvt] = dtmp; /* tau * v2 */
        vt[*n + *n * *ldvt] = 1.0 + (dtmp * d__1); /* 1 - tau * v2^2 */

        /* for HH vectors [n-3:1] */
        for (i = *n - 3; i >= 1; i--)
        {
            stau = - taup[i];

            /* Scale row i by -tau and dlarf for rest of the rows */
            for (j = i + 2; j <= *n; j++)
            {
                vt[i + 1 + j * *ldvt] = stau * a[i + j * *lda];

                /* GEMV part of the dlarf excluding zero first column */
                dtmp = 0.;
                for (k = i + 2; k <= *n; k++)
                {
                    dtmp = dtmp + vt[j + k * *ldvt] * a[i + k * *lda];
                }
                vt[j + (i + 1) * *ldvt] = stau * dtmp;
            }
            vt[i + 1 + (i + 1) * *ldvt] = 1.0 + stau;

            for (j = i + 2; j <= *n; j++)
            {
                for (k = i + 2; k <= *n; k++)
                {
                    vt[j + k * *ldvt] = vt[j + k * *ldvt] + a[i + k * *lda] *
                                        vt[j + (i + 1) * *ldvt];
                }
            }
        }
    }

    /* Generate Ql (from bidiag) in u from a */

    if (*n > 1)
    {
        /* iteration corresponding to (n - 1) HH(n-1) */
        stau = tauq[*n - 1];
        d__1 = a[*n + (*n - 1) * *lda];
        dtmp = - (stau * d__1);

        u[*n - 1 + (*n - 1) * *ldu] = 1.0 - stau; /* 1 - tau */
        u[*n + (*n - 1) * *ldu] = dtmp; /* tau * v2 */
        u[*n - 1 + *n * *ldu] = dtmp; /* tau * v2 */
        u[*n + *n * *ldu] = 1.0 + (dtmp * d__1); /* 1 - tau * v2^2 */
    }
    else
    {
        u[1 + *ldu] = 1.0;
    }

    /* for HH vectors [n-2:1] */
    for (i = *n - 2; i >= 1; i--)
    {
        stau = - tauq[i];

        /* scale col i by -tau and dlarf for rest of the columns */
        for (j = i + 1; j <= *n; j++)
        {
            u[j + i * *ldu] = stau * a[j + i * *lda];

            /* GEMV part of dlarf excluding zero first row */
            dtmp = 0;
            for (k = i + 1; k <= *n; k++)
            {
                dtmp = dtmp + u[k + j * *ldu] * a[k + i * *lda];
            }
            u[i + j * *ldu] = stau * dtmp;
        }
        u[i + i * *ldu] = 1.0 + stau;

        for (j = i + 1; j <= *n; j++)
        {
            for (k = i + 1; k <= *n; k++)
            {
                u[k + j * *ldu] = u[k + j * *ldu] + a[k + i * *lda] * u[i + j * *ldu];
            }
        }
    }
    vt[1 + *ldvt] = 1.0;

    lapack_dbdsqr("U", n, n, n, &c__0, &s[1], &e[1],
                  &vt[1 + *ldvt], ldvt,
                  &u[1 + *ldu], ldu,
                  dum, &c__1,
                  &work[iwork], info);

    /* Apply HH from QR factorization (qr) on vt from left */

    tau = &work[itau - 1];
    /* First Iteration corresponding to HH(n) */
    i = *n;
    for (j = 1; j <= *n; j++)
    {
        /* - u[i][j] * tau[i] */
        d__1 = - u[i + j * *ldu] * tau[i];

        /* u[n+1:m, j] = d__1 * u[n+1:m, j] */
        for (k = *n + 1; k <= *m; k++)
        {
            u[k + j * *ldu] = d__1 * qr[k + *n * *ldqr];
        }
    }
    /* u[m, 1:m] = u[m, 1:m] * (1 - tau) */
    d__1 = 1 - tau[i];
    for (j = 1; j <= *n; j++)
    {
        u[*n + j * *ldu] = u[*n + j * *ldu] * d__1;
    }

    /* Second Iteration onwards */
    beta = 0;
    for (i = *n - 1; i >= 1; i--)
    {
        ni = *n + i;

        au = &u[-i * *ldqr];
        v = &qr[i + i * *ldqr - 1];
        FLA_ELEM_REFLECTOR_APPLY_DLARGE(i, m, &ni, au, ldqr, tau);
    }

    return;
}
#endif
