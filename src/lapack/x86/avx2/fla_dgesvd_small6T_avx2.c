/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_dgesvd_small6T_avx2_.c
 *  @brief DGESVD Small path (path 6T)
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
void fla_dgesvd_small6T_avx2(integer *m, integer *n,
                             doublereal *a, integer *lda,
                             doublereal *ql, integer *ldql,
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

    doublereal *tau, *tauq, *taup;
    doublereal *e, *vtau, *avt;
    doublereal stau, d__1;
    doublereal dum[1];

    /* indices for partitioning work buffer */
    iu = 1;
    itau = iu + *m * *lda;
    ie = itau + *m;
    itauq = ie + *m;
    itaup = itauq + *m;
    iwork = itaup + *m;

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
    FLA_BIDIAGONALIZE_SMALL(*m, *m);

    for (i = 1; i <= *m; i++)
        for (j = 1; j <= *n; j++)
            vt[i + j * *ldvt] = 0.;
    /* Generate Qr (from bidiag) in vt from work[iu] (a here) */
    if (*m > 2)
    {
        /* iteration corresponding to (m - 2) HH[m-2] */
        stau = taup[*m - 2];
        d__1 = a[*m - 2 + *m * *lda];
        dtmp = - (stau * d__1); /* tau * v2 */

        vt[*m - 1 + (*m - 1) * *ldvt] = 1.0 - stau; /* 1 - tau */
        vt[*m + (*m - 1) * *ldvt] = dtmp; /* tau * v2 */
        vt[*m - 1 + *m * *ldvt] = dtmp; /* tau * v2 */
        vt[*m + *m * *ldvt] = 1.0 + (dtmp * d__1); /* 1 - tau * v2^2 */

        /* for HH vectors [m-3:1] */
        for (i = *m - 3; i >= 1; i--)
        {
            stau = - taup[i];

            /* Scale row i by -tau and dlarf for rest of the rows */
            for (j = i + 2; j <= *m; j++)
            {
                vt[i + 1 + j * *ldvt] = stau * a[i + j * *lda];

                /* GEMV part of the dlarf excluding zero first column */
                dtmp = 0.;
                for (k = i + 2; k <= *m; k++)
                {
                    dtmp = dtmp + vt[j + k * *ldvt] * a[i + k * *lda];
                }
                vt[j + (i + 1) * *ldvt] = stau * dtmp;
            }
            vt[i + 1 + (i + 1) * *ldvt] = 1.0 + stau;

            for (j = i + 2; j <= *m; j++)
            {
                for (k = i + 2; k <= *m; k++)
                {
                    vt[j + k * *ldvt] = vt[j + k * *ldvt] + a[i + k * *lda] *
                                        vt[j + (i + 1) * *ldvt];
                }
            }
        }
    }

    /* Generate Ql (from bidiag) in u from a */

    if (*m > 1)
    {
        /* iteration corresponding to (m - 1) HH(m-1) */
        stau = tauq[*m - 1];
        d__1 = a[*m + (*m - 1) * *lda];
        dtmp = - (stau * d__1);

        u[*m - 1 + (*m - 1) * *ldu] = 1.0 - stau; /* 1 - tau */
        u[*m + (*m - 1) * *ldu] = dtmp; /* tau * v2 */
        u[*m - 1 + *m * *ldu] = dtmp; /* tau * v2 */
        u[*m + *m * *ldu] = 1.0 + (dtmp * d__1); /* 1 - tau * v2^2 */
    }
    else
    {
        u[1 + *ldu] = 1.0;
    }

    /* for HH vectors [m-2:1] */
    for (i = *m - 2; i >= 1; i--)
    {
        stau = - tauq[i];

        /* scale col i by -tau and dlarf for rest of the columns */
        for (j = i + 1; j <= *m; j++)
        {
            u[j + i * *ldu] = stau * a[j + i * *lda];

            /* GEMV part of dlarf excluding zero first row */
            dtmp = 0;
            for (k = i + 1; k <= *m; k++)
            {
                dtmp = dtmp + u[k + j * *ldu] * a[k + i * *lda];
            }
            u[i + j * *ldu] = stau * dtmp;
        }
        u[i + i * *ldu] = 1.0 + stau;

        for (j = i + 1; j <= *m; j++)
        {
            for (k = i + 1; k <= *m; k++)
            {
                u[k + j * *ldu] = u[k + j * *ldu] + a[k + i * *lda] * u[i + j * *ldu];
            }
        }
    }
    vt[1 + *ldvt] = 1.0;

    lapack_dbdsqr("U", m, m, m, &c__0, &s[1], &e[1],
                  &vt[1 + *ldvt], ldvt,
                  &u[1 + *ldu], ldu,
                  dum, &c__1,
                  &work[iwork], info);

    tau = &work[itau - 1];
    vtau = tau + *m;
    avt = vtau + *n;

    /* Apply HH from LQ factorization (ql) on vt from right */

    /* First Iteration corresponding to HH(m) */
    i = *m;
    for (j = i + 1; j <= *n; j++)
    {
        /* - ql[i][j] * tau[i] */
        d__1 = - ql[i + j * *ldql] * tau[i];

        /* vt[1:m, j] = d__1 * vt[1:m, j] */
        for (k = 1; k <= *m; k++)
        {
            vt[k + j * *ldvt] = d__1 * vt[k + i * *ldvt];
        }
    }
    /* vt[m, 1:m] = vt[m, 1:m] * (1 - tau) */
    d__1 = 1 - tau[i];
    for (j = 1; j <= *m; j++)
    {
        vt[j + *m * *ldvt] = vt[j + *m * *ldvt] * d__1;
    }

    /* Second Iteration onwards */
    for (i = *m - 1; i >= 1; i--)
    {
        /* Scale HH vector by tau, store in vtau */
        vtau[1] = - tau[i];
        for (j = 2; j <= (*n - i + 1); j++)
        {
            vtau[j] = vtau[1] * ql[i + (j + i - 1) * *ldql];
        }

        /* avt = Vt * vtau (gemv) */
        for (j = 1; j <= *m; j++)
        {
            avt[j] = 0.;
        }
        for (j = 1; j <= (*n - i + 1); j++) /* for every column of Vt */
        {
            for (k = 1; k <= *m; k++) /* Scale the col and accumulate */
            {
                avt[k] = avt[k] + vtau[j] * vt[k + (j + i - 1) * *ldvt];
            }
        }

        /* Vt = Vt + avt * v' (ger) */
        for (k = 1; k <= *m; k++)
        {
            vt[k + i * *ldvt] = vt[k + i * *ldvt] + avt[k];
        }
        for (j = 2; j <= (*n - i + 1); j++)
        {
            for (k = 1; k <= *m; k++)
            {
                vt[k + (j + i - 1) * *ldvt] = vt[k + (j + i - 1) * *ldvt] + avt[k] * ql[i + (j + i - 1) * *ldql];
            }
        }
    }

    return;
}
#endif
