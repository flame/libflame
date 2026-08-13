/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_bdsqr_small_kernel.h
 *  @brief Body of the small-size bidiagonal QR iteration.
 *  */

#define BDSQR_FNAME FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_FNAME)
#define BDSQR_ELT FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_ELT)
#define BDSQR_LAMCH FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_LAMCH)
#define BDSQR_LARTG FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_LARTG)
#define BDSQR_LAS2 FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_LAS2)
#define BDSQR_LASV2 FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_LASV2)
#define BDSQR_SIGN FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_SIGN)
#define BDSQR_ROT FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_ROT)
#define BDSQR_SCAL FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_SCAL)
#define BDSQR_SWAP FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_SWAP)
#define BDSQR_DECL_EXTRA FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_DECL_EXTRA)
#define BDSQR_TOLMUL FLA_BDSQR_SMALL_CAT(BDSQR_PRE, _BDSQR_TOLMUL)

/* Subroutine */
int BDSQR_FNAME(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru, BDSQR_ELT *d__,
                BDSQR_ELT *e, BDSQR_ELT *vt, aocl_int64_t *ldvt, BDSQR_ELT *u, aocl_int64_t *ldu,
                aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    BDSQR_ELT d__1, d__2, d__3, d__4;
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal),
        BDSQR_SIGN(BDSQR_ELT *, BDSQR_ELT *);
    /* Local variables */
    aocl_int64_t iterdivn;
    BDSQR_ELT f, g, h__;
    aocl_int64_t i__, j, m;
    BDSQR_ELT r__;
    aocl_int64_t maxitdivn;
    BDSQR_ELT cs;
    aocl_int64_t ll;
    BDSQR_ELT sn, mu;
    aocl_int64_t tidx, lll;
    BDSQR_ELT eps, sll, tol, abse;
    aocl_int64_t idir;
    BDSQR_ELT abss;
    aocl_int64_t oldm;
    BDSQR_ELT cosl;
    aocl_int64_t isub, iter;
    BDSQR_ELT unfl, sinl, cosr, smin, smax, sinr;
#ifndef FLA_ENABLE_AOCL_BLAS
    void BDSQR_LAS2(BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *);
    extern logical lsame_(char *, char *, aocl_int64_t a, aocl_int64_t b);
#endif
    BDSQR_ELT oldcs;
    aocl_int64_t oldll;
    BDSQR_ELT shift, sigmn, oldsn;
    /* Subroutine */
    BDSQR_ELT sminl, sigmx;
    logical lower;
    extern /* Subroutine */
        void
        BDSQR_LASV2(BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *,
                    BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *);
    extern BDSQR_ELT BDSQR_LAMCH(char *);
    extern /* Subroutine */
        void
        BDSQR_LARTG(BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *, BDSQR_ELT *);
    BDSQR_ELT sminoa, thresh;
    BDSQR_ELT tolmul;
    BDSQR_DECL_EXTRA
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    /* Function Body */
    *info = 0;
    lower = lsame_(uplo, "L", 1, 1);
    if(*n == 1)
    {
        goto L160;
    }
    idir = 0;
    /* Get machine constants */
    eps = BDSQR_LAMCH("Epsilon");
    unfl = BDSQR_LAMCH("Safe minimum");

    /* If matrix lower bidiagonal, rotate to be upper bidiagonal */
    /* by applying Givens rotations on the left */
    if(lower)
    {
        i__1 = *n - 1;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            BDSQR_LARTG(&d__[i__], &e[i__], &cs, &sn, &r__);
            d__[i__] = r__;
            e[i__] = sn * d__[i__ + 1];
            d__[i__ + 1] = cs * d__[i__ + 1];
            /* Update singular vectors if desired */
            if(*nru > 0)
            {
                FLA_APPLY_GIVENS_GRVX(BDSQR_PRE, nru, u, ldu, i__, cs, sn);
            }
        }
    }
    /* Compute singular values to relative accuracy TOL */
    /* (By setting TOL to be negative, algorithm will compute */
    /* singular values to absolute accuracy ABS(TOL)*norm(input matrix)) */
    /* Computing MAX */
    /* Computing MIN */
    BDSQR_TOLMUL
    tol = tolmul * eps;
    /* Compute approximate maximum, minimum singular values */
    smax = 0.;
    sminl = 0.;
    if(tol >= 0.)
    {
        /* Relative accuracy desired */
        sminoa = f2c_abs(d__[1]);
        if(sminoa == 0.)
        {
            goto L50;
        }
        mu = sminoa;
        i__1 = *n;
        for(i__ = 2; i__ <= i__1; ++i__)
        {
            mu = (d__2 = d__[i__], f2c_abs(d__2))
                 * (mu / (mu + (d__1 = e[i__ - 1], f2c_abs(d__1))));
            sminoa = fla_min(sminoa, mu);
            if(sminoa == 0.)
            {
                goto L50;
            }
        }
    L50:
        sminoa /= sqrt((doublereal)(*n));
        /* Computing MAX */
        d__1 = tol * sminoa;
        d__2 = *n * (*n * unfl) * 6; // , expr subst
        thresh = fla_max(d__1, d__2);
    }
    else
    {
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Computing MAX */
            d__2 = smax;
            d__3 = (d__1 = d__[i__], f2c_abs(d__1)); // , expr subst
            smax = fla_max(d__2, d__3);
        }
        i__1 = *n - 1;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Computing MAX */
            d__2 = smax;
            d__3 = (d__1 = e[i__], f2c_abs(d__1)); // , expr subst
            smax = fla_max(d__2, d__3);
        }
        /* Absolute accuracy desired */
        /* Computing MAX */
        d__1 = f2c_abs(tol) * smax;
        d__2 = *n * (*n * unfl) * 6; // , expr subst
        thresh = fla_max(d__1, d__2);
    }
    /* Prepare for main iteration loop for the singular values */
    /* (MAXIT is the maximum number of passes through the inner */
    /* loop permitted before nonconvergence signalled.) */
    maxitdivn = *n * 6;
    iterdivn = 0;
    iter = -1;
    oldll = -1;
    oldm = -1;
    /* M points to last element of unconverged part of matrix */
    m = *n;
/* Begin main iteration loop */
L60: /* Check for convergence or exceeding iteration count */
    if(m <= 1)
    {
        goto L160;
    }
    if(iter >= *n)
    {
        iter -= *n;
        ++iterdivn;
        if(iterdivn >= maxitdivn)
        {
            goto L200;
        }
    }
    /* Find diagonal block of matrix to work on */
    if(tol < 0. && (d__1 = d__[m], f2c_abs(d__1)) <= thresh)
    {
        d__[m] = 0.;
    }
    smax = (d__1 = d__[m], f2c_abs(d__1));
    smin = smax;
    i__1 = m - 1;
    for(lll = 1; lll <= i__1; ++lll)
    {
        ll = m - lll;
        abss = (d__1 = d__[ll], f2c_abs(d__1));
        abse = (d__1 = e[ll], f2c_abs(d__1));
        if(tol < 0. && abss <= thresh)
        {
            d__[ll] = 0.;
        }
        if(abse <= thresh)
        {
            goto L80;
        }
        smin = fla_min(smin, abss);
        /* Computing MAX */
        d__1 = fla_max(smax, abss);
        smax = fla_max(d__1, abse);
        /* L70: */
    }
    ll = 0;
    goto L90;
L80:
    e[ll] = 0.;
    /* Matrix splits since E(LL) = 0 */
    if(ll == m - 1)
    {
        /* Convergence of bottom singular value, return to top of loop */
        --m;
        goto L60;
    }
L90:
    ++ll;
    /* E(LL) through E(M-1) are nonzero, E(LL-1) is zero */
    if(ll == m - 1)
    {
        /* 2 by 2 block, handle separately */
        BDSQR_LASV2(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr, &sinl, &cosl);
        d__[m - 1] = sigmx;
        e[m - 1] = 0.;
        d__[m] = sigmn;
        /* Compute singular vectors, if desired */
        if(*ncvt > 0)
        {
            BDSQR_ROT(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &cosr, &sinr);
        }
        if(*nru > 0)
        {
            BDSQR_ROT(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &c__1, &cosl,
                      &sinl);
        }
        m += -2;
        goto L60;
    }
    /* If working on new submatrix, choose shift direction */
    /* (from larger end diagonal element towards smaller) */
    if(ll > oldm || m < oldll)
    {
        if((d__1 = d__[ll], f2c_abs(d__1)) >= (d__2 = d__[m], f2c_abs(d__2)))
        {
            /* Chase bulge from top (big end) to bottom (small end) */
            idir = 1;
        }
        else
        {
            /* Chase bulge from bottom (big end) to top (small end) */
            idir = 2;
        }
    }
    /* Apply convergence tests */
    if(idir == 1)
    {
        /* Run convergence test in forward direction */
        /* First apply standard test to bottom of matrix */
        if((d__2 = e[m - 1], f2c_abs(d__2)) <= f2c_abs(tol) * (d__1 = d__[m], f2c_abs(d__1))
           || tol < 0. && (d__3 = e[m - 1], f2c_abs(d__3)) <= thresh)
        {
            e[m - 1] = 0.;
            goto L60;
        }
        if(tol >= 0.)
        {
            /* If relative accuracy desired, */
            /* apply convergence criterion forward */
            mu = (d__1 = d__[ll], f2c_abs(d__1));
            sminl = mu;
            i__1 = m - 1;
            for(lll = ll; lll <= i__1; ++lll)
            {
                if((d__1 = e[lll], f2c_abs(d__1)) <= tol * mu)
                {
                    e[lll] = 0.;
                    goto L60;
                }
                mu = (d__2 = d__[lll + 1], f2c_abs(d__2))
                     * (mu / (mu + (d__1 = e[lll], f2c_abs(d__1))));
                sminl = fla_min(sminl, mu);
            }
        }
    }
    else
    {
        /* Run convergence test in backward direction */
        /* First apply standard test to top of matrix */
        if((d__2 = e[ll], f2c_abs(d__2)) <= f2c_abs(tol) * (d__1 = d__[ll], f2c_abs(d__1))
           || tol < 0. && (d__3 = e[ll], f2c_abs(d__3)) <= thresh)
        {
            e[ll] = 0.;
            goto L60;
        }
        if(tol >= 0.)
        {
            /* If relative accuracy desired, */
            /* apply convergence criterion backward */
            mu = (d__1 = d__[m], f2c_abs(d__1));
            sminl = mu;
            i__1 = ll;
            for(lll = m - 1; lll >= i__1; --lll)
            {
                if((d__1 = e[lll], f2c_abs(d__1)) <= tol * mu)
                {
                    e[lll] = 0.;
                    goto L60;
                }
                mu = (d__2 = d__[lll], f2c_abs(d__2))
                     * (mu / (mu + (d__1 = e[lll], f2c_abs(d__1))));
                sminl = fla_min(sminl, mu);
            }
        }
    }
    oldll = ll;
    oldm = m;
    /* Compute shift. First, test if shifting would ruin relative */
    /* accuracy, and if so set the shift to zero. */
    /* Computing MAX */
    d__1 = eps;
    d__2 = tol * .01; // , expr subst
    if(tol >= 0. && *n * tol * (sminl / smax) <= fla_max(d__1, d__2))
    {
        /* Use a zero shift to avoid loss of relative accuracy */
        shift = 0.;
    }
    else
    {
        /* Compute the shift from 2-by-2 block at end of matrix */
        if(idir == 1)
        {
            sll = (d__1 = d__[ll], f2c_abs(d__1));
            BDSQR_LAS2(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
        }
        else
        {
            sll = (d__1 = d__[m], f2c_abs(d__1));
            BDSQR_LAS2(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
        }
        /* Test if shift negligible, and if so set to zero */
        if(sll > 0.)
        {
            /* Computing 2nd power */
            d__1 = shift / sll;
            if(d__1 * d__1 < eps)
            {
                shift = 0.;
            }
        }
    }
    /* Increment iteration count */
    iter = iter + m - ll;
    /* If SHIFT = 0, do simplified QR iteration */
    if(shift == 0.)
    {
        if(idir == 1)
        {
            /* Chase bulge from top to bottom */
            /* Save cosines and sines for later singular vector updates */
            cs = 1.;
            oldcs = 1.;
            oldsn = 0.;
            i__1 = m - 1;
            for(i__ = ll; i__ <= i__1; ++i__)
            {
                d__1 = d__[i__] * cs;
                BDSQR_LARTG(&d__1, &e[i__], &cs, &sn, &r__);
                if(i__ > ll)
                {
                    e[i__ - 1] = oldsn * r__;
                }
                d__1 = oldcs * r__;
                d__2 = d__[i__ + 1] * sn;
                BDSQR_LARTG(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
                if(*ncvt > 0)
                {
                    FLA_APPLY_GIVENS_GLVX(BDSQR_PRE, ncvt, vt, ldvt, i__, cs, sn);
                }
                if(*nru > 0)
                {
                    FLA_APPLY_GIVENS_GRVX(BDSQR_PRE, nru, u, ldu, i__, oldcs, oldsn);
                }
            }
            h__ = d__[m] * cs;
            d__[m] = h__ * oldcs;
            e[m - 1] = h__ * oldsn;
            /* Test convergence */
            if((d__1 = e[m - 1], f2c_abs(d__1)) <= thresh)
            {
                e[m - 1] = 0.;
            }
        }
        else
        {
            /* Chase bulge from bottom to top */
            /* Save cosines and sines for later singular vector updates */
            cs = 1.;
            oldcs = 1.;
            i__1 = ll + 1;
            for(i__ = m; i__ >= i__1; --i__)
            {
                d__1 = d__[i__] * cs;
                BDSQR_LARTG(&d__1, &e[i__ - 1], &cs, &sn, &r__);
                if(i__ < m)
                {
                    e[i__] = oldsn * r__;
                }
                d__1 = oldcs * r__;
                d__2 = d__[i__ - 1] * sn;
                BDSQR_LARTG(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
                if(*ncvt > 0)
                {
                    tidx = i__ - 1;
                    FLA_APPLY_GIVENS_GLVX(BDSQR_PRE, ncvt, vt, ldvt, tidx, oldcs, (-oldsn));
                }
                if(*nru > 0)
                {
                    tidx = i__ - 1;
                    FLA_APPLY_GIVENS_GRVX(BDSQR_PRE, nru, u, ldu, tidx, cs, (-sn));
                }
            }
            h__ = d__[ll] * cs;
            d__[ll] = h__ * oldcs;
            e[ll] = h__ * oldsn;
            /* Test convergence */
            if((d__1 = e[ll], f2c_abs(d__1)) <= thresh)
            {
                e[ll] = 0.;
            }
        }
    }
    else
    {
        /* Use nonzero shift */
        if(idir == 1)
        {
            /* Chase bulge from top to bottom */
            /* Save cosines and sines for later singular vector updates */
            f = ((d__1 = d__[ll], f2c_abs(d__1)) - shift)
                * (BDSQR_SIGN(&c_b49, &d__[ll]) + shift / d__[ll]);
            g = e[ll];
            i__1 = m - 1;
            for(i__ = ll; i__ <= i__1; ++i__)
            {
                BDSQR_LARTG(&f, &g, &cosr, &sinr, &r__);
                if(i__ > ll)
                {
                    e[i__ - 1] = r__;
                }
                f = cosr * d__[i__] + sinr * e[i__];
                e[i__] = cosr * e[i__] - sinr * d__[i__];
                g = sinr * d__[i__ + 1];
                d__[i__ + 1] = cosr * d__[i__ + 1];
                BDSQR_LARTG(&f, &g, &cosl, &sinl, &r__);
                d__[i__] = r__;
                f = cosl * e[i__] + sinl * d__[i__ + 1];
                d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
                if(i__ < m - 1)
                {
                    g = sinl * e[i__ + 1];
                    e[i__ + 1] = cosl * e[i__ + 1];
                }
                if(*ncvt > 0)
                {
                    FLA_APPLY_GIVENS_GLVX(BDSQR_PRE, ncvt, vt, ldvt, i__, cosr, sinr);
                }
                if(*nru > 0)
                {
                    FLA_APPLY_GIVENS_GRVX(BDSQR_PRE, nru, u, ldu, i__, cosl, sinl);
                }
            }
            e[m - 1] = f;
            /* Test convergence */
            if((d__1 = e[m - 1], f2c_abs(d__1)) <= thresh)
            {
                e[m - 1] = 0.;
            }
        }
        else
        {
            /* Chase bulge from bottom to top */
            /* Save cosines and sines for later singular vector updates */
            f = ((d__1 = d__[m], f2c_abs(d__1)) - shift)
                * (BDSQR_SIGN(&c_b49, &d__[m]) + shift / d__[m]);
            g = e[m - 1];
            i__1 = ll + 1;
            for(i__ = m; i__ >= i__1; --i__)
            {
                BDSQR_LARTG(&f, &g, &cosr, &sinr, &r__);
                if(i__ < m)
                {
                    e[i__] = r__;
                }
                f = cosr * d__[i__] + sinr * e[i__ - 1];
                e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
                g = sinr * d__[i__ - 1];
                d__[i__ - 1] = cosr * d__[i__ - 1];
                BDSQR_LARTG(&f, &g, &cosl, &sinl, &r__);
                d__[i__] = r__;
                f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
                d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
                if(i__ > ll + 1)
                {
                    g = sinl * e[i__ - 2];
                    e[i__ - 2] = cosl * e[i__ - 2];
                }
                if(*ncvt > 0)
                {
                    tidx = i__ - 1;
                    FLA_APPLY_GIVENS_GLVX(BDSQR_PRE, ncvt, vt, ldvt, tidx, cosl, (-sinl));
                }
                if(*nru > 0)
                {
                    tidx = i__ - 1;
                    FLA_APPLY_GIVENS_GRVX(BDSQR_PRE, nru, u, ldu, tidx, cosr, (-sinr));
                }
            }
            e[ll] = f;
            /* Test convergence */
            if((d__1 = e[ll], f2c_abs(d__1)) <= thresh)
            {
                e[ll] = 0.;
            }
        }
    }
    /* QR iteration finished, go back and check convergence */
    goto L60;
/* All singular values converged, so make them positive */
L160:
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(d__[i__] == 0.)
        {
            /* Avoid -ZERO */
            d__[i__] = 0.;
        }
        else if(d__[i__] < 0.)
        {
            d__[i__] = -d__[i__];
            /* Change sign of singular vectors, if desired */
            if(*ncvt > 0)
            {
                BDSQR_SCAL(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
            }
        }
    }
    /* Sort the singular values into decreasing order (insertion sort on */
    /* singular values, but only one transposition per singular vector) */
    i__1 = *n - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        /* Scan for smallest D(I) */
        isub = 1;
        smin = d__[1];
        i__2 = *n + 1 - i__;
        for(j = 2; j <= i__2; ++j)
        {
            if(d__[j] <= smin)
            {
                isub = j;
                smin = d__[j];
            }
        }
        if(isub != *n + 1 - i__)
        {
            /* Swap singular values and vectors */
            d__[isub] = d__[*n + 1 - i__];
            d__[*n + 1 - i__] = smin;
            if(*ncvt > 0)
            {
                BDSQR_SWAP(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + vt_dim1], ldvt);
            }
            if(*nru > 0)
            {
                BDSQR_SWAP(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * u_dim1 + 1],
                           &c__1);
            }
        }
    }
    goto L220;
/* Maximum number of iterations exceeded, failure to converge */
L200:
    *info = 0;
    i__1 = *n - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(e[i__] != 0.)
        {
            ++(*info);
        }
    }
L220:
    return 0;
    /* End of BDSQR */
}

#undef BDSQR_FNAME
#undef BDSQR_ELT
#undef BDSQR_LAMCH
#undef BDSQR_LARTG
#undef BDSQR_LAS2
#undef BDSQR_LASV2
#undef BDSQR_SIGN
#undef BDSQR_ROT
#undef BDSQR_SCAL
#undef BDSQR_SWAP
#undef BDSQR_DECL_EXTRA
#undef BDSQR_TOLMUL
