/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_GESVD_SMALL_AVX2_KERNEL_H
#define FLA_GESVD_SMALL_AVX2_KERNEL_H

/*! @file fla_gesvd_small_avx2_kernel.h
 *  @brief SVD kernels for small sizes (real s/d).
 *  */

#if FLA_ENABLE_AMD_OPT

#include "fla_geqrf_small_avx2_kernel.h"

#define FLA_GESVD_SMALL_CAT_(a, b) a##b
#define FLA_GESVD_SMALL_CAT(a, b) FLA_GESVD_SMALL_CAT_(a, b)

#define d_NRM2 aocl_blas_dnrm2
#define d_LAPY2 dlapy2_
#define d_SIGN d_sign
#define d_LASV2 dlasv2_
#define d_LARTG dlartg_

#define s_NRM2 aocl_blas_snrm2
#define s_LAPY2 slapy2_
#define s_SIGN r_sign
#define s_LASV2 slasv2_
#define s_LARTG slartg_

doublereal d_sign(doublereal *, doublereal *);
doublereal r_sign(real *, real *);

/* Generate elementary reflector to annihilate the elements to the right of
 * the super-diagonal of the current row */
#define FLA_LARF_GEN_GSMALL_ROW(P, i, m, n, iptr, ldia, tau)                         \
    /* Compute norm2 */                                                              \
    xnorm = FLA_GESVD_SMALL_CAT(P, _NRM2)(&rlen, &iptr[2 * *ldia], ldia);            \
    if(xnorm == 0.)                                                                  \
    {                                                                                \
        tau[i] = 0.;                                                                 \
        beta = iptr[*ldia];                                                          \
    }                                                                                \
    else                                                                             \
    {                                                                                \
        knt = 0;                                                                     \
        v = iptr;                                                                    \
        alpha = v[*ldia];                                                            \
        d__1 = FLA_GESVD_SMALL_CAT(P, _LAPY2)(&v[*ldia], &xnorm);                    \
        beta = -FLA_GESVD_SMALL_CAT(P, _SIGN)(&d__1, &alpha);                        \
        if(f2c_abs(beta) < safmin)                                                   \
        {                                                                            \
            for(knt = 0; f2c_abs(beta) < safmin && knt < 20; knt++)                  \
            {                                                                        \
                FLA_GESVD_SMALL_CAT(P, _SCAL)(&rlen, &rsafmin, &v[2 * *ldia], ldia); \
                beta *= rsafmin;                                                     \
                alpha *= rsafmin;                                                    \
            }                                                                        \
            /* New BETA is at most 1, at least SAFMIN */                             \
            xnorm = FLA_GESVD_SMALL_CAT(P, _NRM2)(&rlen, &v[2 * *ldia], ldia);       \
            d__1 = FLA_GESVD_SMALL_CAT(P, _LAPY2)(&alpha, &xnorm);                   \
            beta = -FLA_GESVD_SMALL_CAT(P, _SIGN)(&d__1, &alpha);                    \
        }                                                                            \
        tau[i] = (beta - alpha) / beta;                                              \
        d__1 = 1. / (alpha - beta);                                                  \
        FLA_GESVD_SMALL_CAT(P, _SCAL)(&rlen, &d__1, &v[2 * *ldia], ldia);            \
        for(j = 1; j <= knt; ++j)                                                    \
        {                                                                            \
            beta *= safmin;                                                          \
        }                                                                            \
    }

/* Apply the row reflector on A(i+1:nr,i+1:nc) from the right */
#define FLA_LARF_APPLY_SMALL_ROW(i, m, n, iptr, ldia, tau)                   \
    if(xnorm == 0.)                                                          \
    {                                                                        \
        tau[i] = 0.;                                                         \
    }                                                                        \
    else                                                                     \
    {                                                                        \
        /* for every row ac of A(i+1:nr,i+1:nc) */                           \
        ac = iptr;                                                           \
        v[*ldia] = 1;                                                        \
        for(j = 1; j <= slen; j++)                                           \
        {                                                                    \
            dtmp = 0;                                                        \
            /* w = (ac .* v) */                                              \
            for(k = 1; k <= rlen + 1; k++)                                   \
            {                                                                \
                dtmp = dtmp + ac[j + k * *ldia] * v[k * *ldia];              \
            }                                                                \
                                                                             \
            /* (ac .* v) * tau */                                            \
            dtmp = dtmp * tau[i];                                            \
                                                                             \
            /* ac = ac - ac * dtmp */                                        \
            for(k = 1; k <= rlen + 1; k++)                                   \
            {                                                                \
                ac[j + k * *ldia] = ac[j + k * *ldia] - v[k * *ldia] * dtmp; \
            }                                                                \
        }                                                                    \
        v[*ldia] = beta;                                                     \
    }

/* Form U from the left bidiagonalizing reflectors of a square matrix */
#define FLA_LARF_UAPPLY_SMALL_SQR(m, a, lda, tauq, u, ldu, twork, one)          \
    if(*m > 1)                                                                  \
    {                                                                           \
        /* iteration corresponding to (m - 1) HH(m-1) */                        \
        stau = tauq[*m - 1];                                                    \
        d__1 = a[*m + (*m - 1) * *lda];                                         \
        dtmp = -(stau * d__1);                                                  \
        u[*m - 1 + (*m - 1) * *ldu] = (one) - stau;                             \
        u[*m + (*m - 1) * *ldu] = dtmp;                                         \
        u[*m - 1 + *m * *ldu] = dtmp;                                           \
        u[*m + *m * *ldu] = (one) + (dtmp * d__1);                              \
    }                                                                           \
    else                                                                        \
    {                                                                           \
        u[1 + *ldu] = (one);                                                    \
    }                                                                           \
    for(i = *m - 2; i >= 1; i--)                                                \
    {                                                                           \
        stau = -tauq[i];                                                        \
        for(j = i + 1; j <= *m; j++)                                            \
        {                                                                       \
            twork[j] = a[j + i * *lda];                                         \
            dtmp = 0;                                                           \
            for(k = i + 1; k <= *m; k++)                                        \
            {                                                                   \
                dtmp = dtmp + u[k + j * *ldu] * a[k + i * *lda];                \
            }                                                                   \
            u[i + j * *ldu] = stau * dtmp;                                      \
        }                                                                       \
        u[i + i * *ldu] = (one) + stau;                                         \
        for(j = i + 1; j <= *m; j++)                                            \
        {                                                                       \
            for(k = i + 1; k <= *m; k++)                                        \
            {                                                                   \
                u[k + j * *ldu] = u[k + j * *ldu] + twork[k] * u[i + j * *ldu]; \
            }                                                                   \
        }                                                                       \
        for(j = i + 1; j <= *m; j++)                                            \
        {                                                                       \
            u[j + i * *ldu] = stau * a[j + i * *lda];                           \
        }                                                                       \
    }

/* Form Vt from the right bidiagonalizing reflectors of a square matrix */
#define FLA_LARF_VTAPPLY_SMALL_SQR(m, a, lda, taup, vt, ldvt, one)                       \
    if(*m > 2)                                                                           \
    {                                                                                    \
        /* iteration corresponding to (m - 2) HH[m-2] */                                 \
        stau = taup[*m - 2];                                                             \
        d__1 = a[*m - 2 + *m * *lda];                                                    \
        dtmp = -(stau * d__1);                                                           \
        vt[*m - 1 + (*m - 1) * *ldvt] = (one) - stau;                                    \
        vt[*m + (*m - 1) * *ldvt] = dtmp;                                                \
        vt[*m - 1 + *m * *ldvt] = dtmp;                                                  \
        vt[*m + *m * *ldvt] = (one) + (dtmp * d__1);                                     \
        for(i = *m - 3; i >= 1; i--)                                                     \
        {                                                                                \
            stau = -taup[i];                                                             \
            for(j = i + 2; j <= *m; j++)                                                 \
            {                                                                            \
                vt[i + 1 + j * *ldvt] = stau * a[i + j * *lda];                          \
                dtmp = 0.;                                                               \
                for(k = i + 2; k <= *m; k++)                                             \
                {                                                                        \
                    dtmp = dtmp + vt[j + k * *ldvt] * a[i + k * *lda];                   \
                }                                                                        \
                vt[j + (i + 1) * *ldvt] = stau * dtmp;                                   \
            }                                                                            \
            vt[i + 1 + (i + 1) * *ldvt] = (one) + stau;                                  \
            for(j = i + 2; j <= *m; j++)                                                 \
            {                                                                            \
                for(k = i + 2; k <= *m; k++)                                             \
                {                                                                        \
                    vt[j + k * *ldvt]                                                    \
                        = vt[j + k * *ldvt] + a[i + k * *lda] * vt[j + (i + 1) * *ldvt]; \
                }                                                                        \
            }                                                                            \
        }                                                                                \
    }                                                                                    \
    else                                                                                 \
    {                                                                                    \
        for(i = 1; i <= *m; i++)                                                         \
        {                                                                                \
            vt[i + i * *ldvt] = (one);                                                   \
        }                                                                                \
    }                                                                                    \
    vt[1 + *ldvt] = (one);

/* Apply a row reflector while forming Vt, for the lower-bidiagonal path */
#define FLA_LARF_VTAPPLY_SMALL_ROW(i, m, n, tau, sv, ldsv)              \
    /* for every row ac of A(i+1:nr,i+1:nc) */                          \
    v[*lda] = 1;                                                        \
    for(j = 1; j <= slen; j++)                                          \
    {                                                                   \
        dtmp = 0;                                                       \
        /* w = (ac .* v) */                                             \
        for(k = 1; k <= rlen + 1; k++)                                  \
        {                                                               \
            dtmp = dtmp + sv[j + k * *ldsv] * v[k * *lda];              \
        }                                                               \
                                                                        \
        /* (ac .* v) * tau */                                           \
        dtmp = dtmp * tau[i];                                           \
                                                                        \
        /* ac = ac - ac * dtmp */                                       \
        for(k = 1; k <= rlen + 1; k++)                                  \
        {                                                               \
            sv[j + k * *ldsv] = sv[j + k * *ldsv] - v[k * *lda] * dtmp; \
        }                                                               \
    }                                                                   \
    v[*lda] = beta;

/* Upper bidiagonalization of an nr x nc block */
#define FLA_BIDIAGONALIZE_GSMALL(P, nr, nc, ia, ldia, qtau, ptau, dv, ev)  \
    for(i = 1; i <= fla_min(nr, nc); i++)                                  \
    {                                                                      \
        slen = nr - i;                                                     \
        /* input address */                                                \
        FLA_GESVD_SMALL_CAT(P, _ELT) * iptr;                               \
        aocl_int64_t has_outliers = 0;                                     \
                                                                           \
        /* Annihilate elements in current column */                        \
        iptr = (FLA_GESVD_SMALL_CAT(P, _ELT) *)&ia[i + 1 + i * *ldia - 1]; \
        if(slen == 0)                                                      \
        {                                                                  \
            qtau[i] = 0.;                                                  \
            beta = 0.;                                                     \
        }                                                                  \
        else if(slen < 4)                                                  \
        {                                                                  \
            /* Generate elementary reflector to annihilate                 \
             * elements below diagonal A(i+1:nr,i) */                      \
            FLA_LARF_GEN_GSMALL_COL(P, i, &nr, &nc, qtau);                 \
            /* Apply the reflector on A(i:nr,i+1:nc) from the left */      \
            FLA_LARF_APPLY_GSMALL_COL(P, i, &nr, &nc, ia, ldia, qtau);     \
        }                                                                  \
        else                                                               \
        {                                                                  \
            /* Generate elementary reflector to annihilate                 \
             * elements below diagonal A(i+1:nr,i) */                      \
            FLA_LARF_GEN_GLARGE_COL(P, i, &nr, &nc, qtau);                 \
            /* Apply the reflector on A(i:nr,i+1:nc) from the left */      \
            FLA_LARF_APPLY_GLARGE_COL(P, i, &nr, &nc, ia, ldia, qtau);     \
        }                                                                  \
        dv[i] = *iptr;                                                     \
                                                                           \
        /* Annihilate elements in current row */                           \
        beta = 0.;                                                         \
        rlen = nc - i - 1;                                                 \
        tau = ptau;                                                        \
        if(rlen <= 0)                                                      \
        {                                                                  \
            tau[i] = 0.;                                                   \
        }                                                                  \
        else                                                               \
        {                                                                  \
            /* Generate elementary reflector to annihilate                 \
             * elements to the right of current row's                      \
             * super diagonalA(i,i+2:nr) */                                \
            FLA_LARF_GEN_GSMALL_ROW(P, i, &nr, &nc, iptr, ldia, tau);      \
            /* Apply the reflector on A(i+1:nr,i+1:nc) from the right */   \
            FLA_LARF_APPLY_SMALL_ROW(i, &nr, &nc, iptr, ldia, tau);        \
        }                                                                  \
        if(rlen >= 0)                                                      \
            ev[i] = iptr[*ldia];                                           \
    }

/* Vt for the 2x2 case, with the singular value signs folded in */
#define FLA_COMPUTE_VT_G2X2(P, vt, ldvt, s1, s2, cr, sr) \
    FLA_GESVD_SMALL_CAT(P, _ELT) scl1, scl2;             \
                                                         \
    scl1 = (s1 < 0.) ? -1. : 1.;                         \
    scl2 = (s2 < 0.) ? -1. : 1.;                         \
                                                         \
    vt[1 + *ldvt] = scl1 * cr;                           \
    vt[2 + *ldvt] = scl2 * -sr;                          \
                                                         \
    vt[1 + 2 * *ldvt] = scl1 * sr;                       \
    vt[2 + 2 * *ldvt] = scl2 * cr;

/* Macro to normalize singular value sign and
   scale corresponding singular vectors for 2x2 matrices */
#define FLA_NORMALIZE_SINGULAR_VALUE_AND_VECTORS_2X2(idx, wntu_var) \
    if(s[idx] == 0.0)                                               \
    {                                                               \
        s[idx] = 0.0; /* Avoid -ZERO */                             \
    }                                                               \
    else if(s[idx] < 0.0)                                           \
    {                                                               \
        s[idx] = -s[idx]; /* Make singular value positive */        \
        if(wntu_var && u != NULL)                                   \
        {                                                           \
            /* Negate corresponding left singular vector column */  \
            u[1 + (idx) * *ldu] = -u[1 + (idx) * *ldu];             \
            u[2 + (idx) * *ldu] = -u[2 + (idx) * *ldu];             \
        }                                                           \
    }

/* Macro to ensure all singular values are positive
   (values only, no vector adjustments) */
#define FLA_ENSURE_POSITIVE_SINGULAR_VALUES(n_vals)   \
    for(aocl_int64_t i__ = 1; i__ <= (n_vals); i__++) \
    {                                                 \
        if(s[i__] == 0.0)                             \
        {                                             \
            s[i__] = 0.0; /* Avoid -ZERO */           \
        }                                             \
        else if(s[i__] < 0.0)                         \
        {                                             \
            s[i__] = -s[i__]; /* Make positive */     \
        }                                             \
    }

#define FLA_GESVD_ZERO_MAT(nr, nc, x, ldx) \
    for(i = 1; i <= (nr); i++)             \
        for(j = 1; j <= (nc); j++)         \
            x[i + j * *ldx] = 0.;

#define FLA_BIDIAG_2X2_GIVENS_GSMALL(P, a, lda, s, e, cosu, sinu)      \
    FLA_GESVD_SMALL_CAT(P, _ELT) s0;                                   \
                                                                       \
    FLA_GESVD_SMALL_CAT(P, _LARTG)                                     \
    (&a[1 + *lda], &a[2 + *lda], &cosu, &sinu, &s0);                   \
    s[1] = s0;                                                         \
                                                                       \
    /* Update 2nd columns of A */                                      \
    dtmp = cosu * a[1 + 2 * *lda] + sinu * a[2 + 2 * *lda];            \
    a[2 + 2 * *lda] = cosu * a[2 + 2 * *lda] - sinu * a[1 + 2 * *lda]; \
    a[1 + 2 * *lda] = dtmp;                                            \
                                                                       \
    /* Update Singular values and vectors */                           \
    s[2] = a[2 + 2 * *lda];                                            \
    e[1] = a[1 + 2 * *lda];

#define FLA_GESVD_LASV2_2X2_GSMALL(P, s, e, sigmn, sigmx, sinr, cosr, sinl, cosl) \
    FLA_GESVD_SMALL_CAT(P, _ELT) sigmn, sigmx, sinr, cosr, sinl, cosl;            \
                                                                                  \
    FLA_GESVD_SMALL_CAT(P, _LASV2)                                                \
    (&s[1], &e[1], &s[2], &sigmn, &sigmx, &sinr, &cosr, &sinl, &cosl);            \
    s[1] = f2c_abs(sigmx);                                                        \
    s[2] = f2c_abs(sigmn);

#define FLA_GESVD_U_FROM_2X2_GIVENS_GSMALL(P, u, ldu, cosl, sinl, cosu, sinu) \
    FLA_GESVD_SMALL_CAT(P, _ELT) p0, p1, p2, p3;                              \
                                                                              \
    p0 = cosl * cosu;                                                         \
    p1 = sinl * sinu;                                                         \
    p2 = sinl * cosu;                                                         \
    p3 = cosl * sinu;                                                         \
                                                                              \
    u[1 + *ldu] = p0 - p1;                                                    \
    u[2 + *ldu] = p2 + p3;                                                    \
                                                                              \
    u[1 + 2 * *ldu] = -(p3 + p2);                                             \
    u[2 + 2 * *ldu] = p0 - p1;

#define FLA_GESVD_UAPPLY_QR_GSMALL(P, m, n, u, ldu, qr, ldqr, tau) \
    /* First Iteration corresponding to HH(n) */                   \
    i = *n;                                                        \
    for(j = 1; j <= *n; j++)                                       \
    {                                                              \
        /* - u[i][j] * tau[i] */                                   \
        d__1 = -u[i + j * *ldu] * tau[i];                          \
                                                                   \
        /* u[n+1:m, j] = d__1 * u[n+1:m, j] */                     \
        for(k = *n + 1; k <= *m; k++)                              \
        {                                                          \
            u[k + j * *ldu] = d__1 * qr[k + *n * *ldqr];           \
        }                                                          \
    }                                                              \
    /* u[m, 1:m] = u[m, 1:m] * (1 - tau) */                        \
    d__1 = 1 - tau[i];                                             \
    for(j = 1; j <= *n; j++)                                       \
    {                                                              \
        u[*n + j * *ldu] = u[*n + j * *ldu] * d__1;                \
    }                                                              \
                                                                   \
    /* Second Iteration onwards */                                 \
    beta = 0;                                                      \
    xnorm = 1.;                                                    \
    for(i = *n - 1; i >= 1; i--)                                   \
    {                                                              \
        /* incrementing n by i to compensate for decrement         \
         * by i done in FLA_LARF_APPLY_{D,S}LARGE_COL              \
         */                                                        \
        ni = *n + i;                                               \
                                                                   \
        au = &u[-i * *ldu];                                        \
        v = &qr[i + i * *ldqr - 1];                                \
        FLA_LARF_APPLY_GLARGE_COL(P, i, m, &ni, au, ldu, tau);     \
    }

#define FLA_GESVD_VTAPPLY_LQ_SMALL(m, n, vt, ldvt, ql, ldql, tau, vtau, avt)              \
    /* First Iteration corresponding to HH(m) */                                          \
    i = *m;                                                                               \
    for(j = i + 1; j <= *n; j++)                                                          \
    {                                                                                     \
        /* - ql[i][j] * tau[i] */                                                         \
        d__1 = -ql[i + j * *ldql] * tau[i];                                               \
                                                                                          \
        /* vt[1:m, j] = d__1 * vt[1:m, j] */                                              \
        for(k = 1; k <= *m; k++)                                                          \
        {                                                                                 \
            vt[k + j * *ldvt] = d__1 * vt[k + i * *ldvt];                                 \
        }                                                                                 \
    }                                                                                     \
    /* vt[m, 1:m] = vt[m, 1:m] * (1 - tau) */                                             \
    d__1 = 1 - tau[i];                                                                    \
    for(j = 1; j <= *m; j++)                                                              \
    {                                                                                     \
        vt[j + *m * *ldvt] = vt[j + *m * *ldvt] * d__1;                                   \
    }                                                                                     \
                                                                                          \
    /* Second Iteration onwards */                                                        \
    for(i = *m - 1; i >= 1; i--)                                                          \
    {                                                                                     \
        /* Scale HH vector by tau, store in vtau */                                       \
        vtau[1] = -tau[i];                                                                \
        for(j = 2; j <= (*n - i + 1); j++)                                                \
        {                                                                                 \
            vtau[j] = vtau[1] * ql[i + (j + i - 1) * *ldql];                              \
        }                                                                                 \
                                                                                          \
        /* avt = Vt * vtau (gemv) */                                                      \
        for(j = 1; j <= *m; j++)                                                          \
        {                                                                                 \
            avt[j] = 0.;                                                                  \
        }                                                                                 \
        for(j = 1; j <= (*n - i + 1); j++) /* for every column of Vt */                   \
        {                                                                                 \
            for(k = 1; k <= *m; k++) /* Scale the col and accumulate */                   \
            {                                                                             \
                avt[k] = avt[k] + vtau[j] * vt[k + (j + i - 1) * *ldvt];                  \
            }                                                                             \
        }                                                                                 \
                                                                                          \
        /* Vt = Vt + avt * v' (ger) */                                                    \
        for(k = 1; k <= *m; k++)                                                          \
        {                                                                                 \
            vt[k + i * *ldvt] = vt[k + i * *ldvt] + avt[k];                               \
        }                                                                                 \
        for(j = 2; j <= (*n - i + 1); j++)                                                \
        {                                                                                 \
            for(k = 1; k <= *m; k++)                                                      \
            {                                                                             \
                vt[k + (j + i - 1) * *ldvt]                                               \
                    = vt[k + (j + i - 1) * *ldvt] + avt[k] * ql[i + (j + i - 1) * *ldql]; \
            }                                                                             \
        }                                                                                 \
    }

#define FLA_GESVD_FORM_U_RECT_SMALL(m, n, ncu, a, lda, tauq, u, ldu)        \
    /* Initialize columns n to ncu of U to eye */                           \
    for(i = *n + 1; i <= *ncu; i++)                                         \
    {                                                                       \
        for(j = *n + 1; j <= *ncu; j++)                                     \
        {                                                                   \
            u[i + j * *ldu] = 0.;                                           \
        }                                                                   \
        u[i + i * *ldu] = 1.;                                               \
    }                                                                       \
    /* for all HH vectors from the end */                                   \
    for(i = *n; i >= 1; i--)                                                \
    {                                                                       \
        /* Update current column */                                         \
        stau = -tauq[i];                                                    \
        for(j = i + 1; j <= *m; j++)                                        \
        {                                                                   \
            u[j + i * *ldu] = stau * a[j + i * *lda];                       \
        }                                                                   \
        u[i + i * *ldu] = 1 + stau;                                         \
                                                                            \
        /* Update rest of the columns from (i + 1) to n */                  \
        for(k = i + 1; k <= *ncu; k++)                                      \
        {                                                                   \
            dtmp = 0.;                                                      \
            for(j = i + 1; j <= *m; j++)                                    \
            {                                                               \
                dtmp = dtmp + a[j + i * *lda] * u[j + k * *ldu];            \
            }                                                               \
            dtmp = stau * dtmp;                                             \
                                                                            \
            for(j = i + 1; j <= *m; j++)                                    \
            {                                                               \
                u[j + k * *ldu] = u[j + k * *ldu] + dtmp * a[j + i * *lda]; \
            }                                                               \
            u[i + k * *ldu] = dtmp;                                         \
        }                                                                   \
    }

#define FLA_GESVD_FORM_VT_LOWER_SMALL(m, n, a, lda, taup, vt, ldvt) \
    xnorm = 1.0;                                                    \
    for(i = *m; i >= 1; i--)                                        \
    {                                                               \
        /* Update current row */                                    \
        for(j = 1; j <= i - 1; j++)                                 \
        {                                                           \
            vt[i + j * *ldvt] = 0.;                                 \
        }                                                           \
        vt[i + i * *ldvt] = 1 - taup[i];                            \
        for(j = i + 1; j <= *n; j++)                                \
        {                                                           \
            vt[i + j * *ldvt] = -taup[i] * a[i + j * *lda];         \
        }                                                           \
                                                                    \
        /* Update rows below current row using row-apply */         \
        v = &a[i + i * *lda - *lda];                                \
        beta = v[*lda];                                             \
        slen = *m - i;                                              \
        rlen = *n - i;                                              \
        ta = &vt[i + i * *ldvt - *ldvt];                            \
        FLA_LARF_VTAPPLY_SMALL_ROW(i, m, n, taup, ta, ldvt);        \
    }

#define FLA_GESVD_FORM_U_LOWER_SMALL(m, a, lda, tauq, u, ldu)                              \
    u[1 + *ldu] = 1.0;                                                                     \
    if(*m > 2)                                                                             \
    {                                                                                      \
        /* iteration corresponding to (m - 2) HH(m-2) */                                   \
        i = *m - 2;                                                                        \
        stau = tauq[i];                                                                    \
        d__1 = a[*m + i * *lda];                                                           \
        dtmp = -(stau * d__1);                                                             \
                                                                                           \
        u[*m - 1 + (*m - 1) * *ldu] = 1.0 - stau; /* 1 - tau */                            \
        u[*m + (*m - 1) * *ldu] = dtmp; /* tau * v2 */                                     \
        u[*m - 1 + *m * *ldu] = dtmp; /* tau * v2 */                                       \
        u[*m + *m * *ldu] = 1.0 + (dtmp * d__1); /* 1 - tau * v2^2 */                      \
                                                                                           \
        u[*m + *ldu] = 0.;                                                                 \
        u[*m - 1 + *ldu] = 0.;                                                             \
        u[1 + *m * *ldu] = 0.;                                                             \
        u[1 + (*m - 1) * *ldu] = 0.;                                                       \
    }                                                                                      \
    else if(*m > 1)                                                                        \
    {                                                                                      \
        /* 2x2 case where where U is identity */                                           \
        u[1 + *ldu] = 1.0;                                                                 \
        u[2 + *ldu] = 0.;                                                                  \
        u[1 + 2 * *ldu] = 0.;                                                              \
        u[2 + 2 * *ldu] = 1.0;                                                             \
    }                                                                                      \
    /* for HH vectors [m-3:1] */                                                           \
    for(i = *m - 3; i >= 1; i--)                                                           \
    {                                                                                      \
        stau = -tauq[i];                                                                   \
        /* scale col (i + 1) by -tau and larf for rest of the columns */                   \
        for(j = i + 2; j <= *m; j++)                                                       \
        {                                                                                  \
            u[j + (i + 1) * *ldu] = stau * a[j + i * *lda];                                \
        }                                                                                  \
        /* Columns (i + 2) to m */                                                         \
        for(j = i + 2; j <= *m; j++)                                                       \
        {                                                                                  \
            /* GEMV part of larf excluding zero first row .                                \
               Store the dot product in u.                                                 \
            */                                                                             \
            dtmp = 0;                                                                      \
            for(k = i + 2; k <= *m; k++)                                                   \
            {                                                                              \
                dtmp = dtmp + u[k + j * *ldu] * a[k + i * *lda];                           \
            }                                                                              \
            u[i + 1 + j * *ldu] = stau * dtmp;                                             \
        }                                                                                  \
        u[i + 1 + (i + 1) * *ldu] = 1.0 + stau;                                            \
                                                                                           \
        for(j = i + 2; j <= *m; j++)                                                       \
        {                                                                                  \
            for(k = i + 2; k <= *m; k++)                                                   \
            {                                                                              \
                u[k + j * *ldu] = u[k + j * *ldu] + a[k + i * *lda] * u[i + 1 + j * *ldu]; \
            }                                                                              \
        }                                                                                  \
                                                                                           \
        /* Initialize 1st row/col elements */                                              \
        u[i + 1 + *ldu] = 0.;                                                              \
        u[1 + (i + 1) * *ldu] = 0.;                                                        \
    }

/* Double precision aliases */
#define FLA_LARF_GEN_DSMALL_ROW(i, m, n, iptr, ldia, tau) \
    FLA_LARF_GEN_GSMALL_ROW(d, i, m, n, iptr, ldia, tau)
#define FLA_BIDIAGONALIZE_DSMALL(nr, nc, ia, ldia, qtau, ptau, dv, ev) \
    FLA_BIDIAGONALIZE_GSMALL(d, nr, nc, ia, ldia, qtau, ptau, dv, ev)
#define FLA_COMPUTE_VT_D2X2(vt, ldvt, s1, s2, cr, sr) \
    FLA_COMPUTE_VT_G2X2(d, vt, ldvt, s1, s2, cr, sr)
#define FLA_BIDIAG_2X2_GIVENS_DSMALL(a, lda, s, e, cosu, sinu) \
    FLA_BIDIAG_2X2_GIVENS_GSMALL(d, a, lda, s, e, cosu, sinu)
#define FLA_GESVD_LASV2_2X2_DSMALL(s, e, sigmn, sigmx, sinr, cosr, sinl, cosl) \
    FLA_GESVD_LASV2_2X2_GSMALL(d, s, e, sigmn, sigmx, sinr, cosr, sinl, cosl)
#define FLA_GESVD_U_FROM_2X2_GIVENS_DSMALL(u, ldu, cosl, sinl, cosu, sinu) \
    FLA_GESVD_U_FROM_2X2_GIVENS_GSMALL(d, u, ldu, cosl, sinl, cosu, sinu)
#define FLA_GESVD_UAPPLY_QR_DSMALL(m, n, u, ldu, qr, ldqr, tau) \
    FLA_GESVD_UAPPLY_QR_GSMALL(d, m, n, u, ldu, qr, ldqr, tau)

/* Single precision aliases */
#define FLA_LARF_GEN_SSMALL_ROW(i, m, n, iptr, ldia, tau) \
    FLA_LARF_GEN_GSMALL_ROW(s, i, m, n, iptr, ldia, tau)
#define FLA_BIDIAGONALIZE_SSMALL(nr, nc, ia, ldia, qtau, ptau, dv, ev) \
    FLA_BIDIAGONALIZE_GSMALL(s, nr, nc, ia, ldia, qtau, ptau, dv, ev)
#define FLA_COMPUTE_VT_S2X2(vt, ldvt, s1, s2, cr, sr) \
    FLA_COMPUTE_VT_G2X2(s, vt, ldvt, s1, s2, cr, sr)
#define FLA_BIDIAG_2X2_GIVENS_SSMALL(a, lda, sv, e, cosu, sinu) \
    FLA_BIDIAG_2X2_GIVENS_GSMALL(s, a, lda, sv, e, cosu, sinu)
#define FLA_GESVD_LASV2_2X2_SSMALL(sv, e, sigmn, sigmx, sinr, cosr, sinl, cosl) \
    FLA_GESVD_LASV2_2X2_GSMALL(s, sv, e, sigmn, sigmx, sinr, cosr, sinl, cosl)
#define FLA_GESVD_U_FROM_2X2_GIVENS_SSMALL(u, ldu, cosl, sinl, cosu, sinu) \
    FLA_GESVD_U_FROM_2X2_GIVENS_GSMALL(s, u, ldu, cosl, sinl, cosu, sinu)
#define FLA_GESVD_UAPPLY_QR_SSMALL(m, n, u, ldu, qr, ldqr, tau) \
    FLA_GESVD_UAPPLY_QR_GSMALL(s, m, n, u, ldu, qr, ldqr, tau)

#endif /* FLA_ENABLE_AMD_OPT */
#endif /* FLA_GESVD_SMALL_AVX2_KERNEL_H */
