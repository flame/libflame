/* ../netlib/clar1v.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLAR1V computes the (scaled) r-th column of the inverse of the submatrix in rows b1 through bn of the tridiagonal matrix LDLT - λI. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAR1V + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clar1v. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clar1v. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clar1v. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD, */
/* PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA, */
/* R, ISUPPZ, NRMINV, RESID, RQCORR, WORK ) */
/* .. Scalar Arguments .. */
/* LOGICAL WANTNC */
/* INTEGER B1, BN, N, NEGCNT, R */
/* REAL GAPTOL, LAMBDA, MINGMA, NRMINV, PIVMIN, RESID, */
/* $ RQCORR, ZTZ */
/* .. */
/* .. Array Arguments .. */
/* INTEGER ISUPPZ( * ) */
/* REAL D( * ), L( * ), LD( * ), LLD( * ), */
/* $ WORK( * ) */
/* COMPLEX Z( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAR1V computes the (scaled) r-th column of the inverse of */
/* > the sumbmatrix in rows B1 through BN of the tridiagonal matrix */
/* > L D L**T - sigma I. When sigma is close to an eigenvalue, the */
/* > computed vector is an accurate eigenvector. Usually, r corresponds */
/* > to the index where the eigenvector is largest in magnitude. */
/* > The following steps accomplish this computation : */
/* > (a) Stationary qd transform, L D L**T - sigma I = L(+) D(+) L(+)**T, */
/* > (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T, */
/* > (c) Computation of the diagonal elements of the inverse of */
/* > L D L**T - sigma I by combining the above transforms, and choosing */
/* > r as the index where the diagonal of the inverse is (one of the) */
/* > largest in magnitude. */
/* > (d) Computation of the (scaled) r-th column of the inverse using the */
/* > twisted factorization obtained by combining the top part of the */
/* > the stationary and the bottom part of the progressive transform. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix L D L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] B1 */
/* > \verbatim */
/* > B1 is INTEGER */
/* > First index of the submatrix of L D L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] BN */
/* > \verbatim */
/* > BN is INTEGER */
/* > Last index of the submatrix of L D L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LAMBDA */
/* > \verbatim */
/* > LAMBDA is REAL */
/* > The shift. In order to compute an accurate eigenvector, */
/* > LAMBDA should be a good approximation to an eigenvalue */
/* > of L D L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is REAL array, dimension (N-1) */
/* > The (n-1) subdiagonal elements of the unit bidiagonal matrix */
/* > L, in elements 1 to N-1. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The n diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] LD */
/* > \verbatim */
/* > LD is REAL array, dimension (N-1) */
/* > The n-1 elements L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] LLD */
/* > \verbatim */
/* > LLD is REAL array, dimension (N-1) */
/* > The n-1 elements L(i)*L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* > PIVMIN is REAL */
/* > The minimum pivot in the Sturm sequence. */
/* > \endverbatim */
/* > */
/* > \param[in] GAPTOL */
/* > \verbatim */
/* > GAPTOL is REAL */
/* > Tolerance that indicates when eigenvector entries are negligible */
/* > w.r.t. their contribution to the residual. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (N) */
/* > On input, all entries of Z must be set to 0. */
/* > On output, Z contains the (scaled) r-th column of the */
/* > inverse. The scaling is such that Z(R) equals 1. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTNC */
/* > \verbatim */
/* > WANTNC is LOGICAL */
/* > Specifies whether NEGCNT has to be computed. */
/* > \endverbatim */
/* > */
/* > \param[out] NEGCNT */
/* > \verbatim */
/* > NEGCNT is INTEGER */
/* > If WANTNC is .TRUE. then NEGCNT = the number of pivots < pivmin */
/* > in the matrix factorization L D L**T, and NEGCNT = -1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] ZTZ */
/* > \verbatim */
/* > ZTZ is REAL */
/* > The square of the 2-norm of Z. */
/* > \endverbatim */
/* > */
/* > \param[out] MINGMA */
/* > \verbatim */
/* > MINGMA is REAL */
/* > The reciprocal of the largest (in magnitude) diagonal */
/* > element of the inverse of L D L**T - sigma I. */
/* > \endverbatim */
/* > */
/* > \param[in,out] R */
/* > \verbatim */
/* > R is INTEGER */
/* > The twist index for the twisted factorization used to */
/* > compute Z. */
/* > On input, 0 <= R <= N. If R is input as 0, R is set to */
/* > the index where (L D L**T - sigma I)^{
-1}
is largest */
/* > in magnitude. If 1 <= R <= N, R is unchanged. */
/* > On output, R contains the twist index used to compute Z. */
/* > Ideally, R designates the position of the maximum entry in the */
/* > eigenvector. */
/* > \endverbatim */
/* > */
/* > \param[out] ISUPPZ */
/* > \verbatim */
/* > ISUPPZ is INTEGER array, dimension (2) */
/* > The support of the vector in Z, i.e., the vector Z is */
/* > nonzero only in elements ISUPPZ(1) through ISUPPZ( 2 ). */
/* > \endverbatim */
/* > */
/* > \param[out] NRMINV */
/* > \verbatim */
/* > NRMINV is REAL */
/* > NRMINV = 1/SQRT( ZTZ ) */
/* > \endverbatim */
/* > */
/* > \param[out] RESID */
/* > \verbatim */
/* > RESID is REAL */
/* > The residual of the FP vector. */
/* > RESID = ABS( MINGMA )/SQRT( ZTZ ) */
/* > \endverbatim */
/* > */
/* > \param[out] RQCORR */
/* > \verbatim */
/* > RQCORR is REAL */
/* > The Rayleigh Quotient correction to LAMBDA. */
/* > RQCORR = MINGMA*TMP */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (4*N) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Beresford Parlett, University of California, Berkeley, USA \n */
/* > Jim Demmel, University of California, Berkeley, USA \n */
/* > Inderjit Dhillon, University of Texas, Austin, USA \n */
/* > Osni Marques, LBNL/NERSC, USA \n */
/* > Christof Voemel, University of California, Berkeley, USA */
/* ===================================================================== */
/* Subroutine */
int clar1v_(integer *n, integer *b1, integer *bn, real * lambda, real *d__, real *l, real *ld, real *lld, real *pivmin, real * gaptol, complex *z__, logical *wantnc, integer *negcnt, real *ztz, real *mingma, integer *r__, integer *isuppz, real *nrminv, real * resid, real *rqcorr, real *work)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1;
    complex q__1, q__2;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    integer i__;
    real s;
    integer r1, r2;
    real eps, tmp;
    integer neg1, neg2, indp, inds;
    real dplus;
    extern real slamch_(char *);
    integer indlpl, indumn;
    extern logical sisnan_(real *);
    real dminus;
    logical sawnan1, sawnan2;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --work;
    --isuppz;
    --z__;
    --lld;
    --ld;
    --l;
    --d__;
    /* Function Body */
    eps = slamch_("Precision");
    if (*r__ == 0)
    {
        r1 = *b1;
        r2 = *bn;
    }
    else
    {
        r1 = *r__;
        r2 = *r__;
    }
    /* Storage for LPLUS */
    indlpl = 0;
    /* Storage for UMINUS */
    indumn = *n;
    inds = (*n << 1) + 1;
    indp = *n * 3 + 1;
    if (*b1 == 1)
    {
        work[inds] = 0.f;
    }
    else
    {
        work[inds + *b1 - 1] = lld[*b1 - 1];
    }
    /* Compute the stationary transform (using the differential form) */
    /* until the index R2. */
    sawnan1 = FALSE_;
    neg1 = 0;
    s = work[inds + *b1 - 1] - *lambda;
    i__1 = r1 - 1;
    for (i__ = *b1;
            i__ <= i__1;
            ++i__)
    {
        dplus = d__[i__] + s;
        work[indlpl + i__] = ld[i__] / dplus;
        if (dplus < 0.f)
        {
            ++neg1;
        }
        work[inds + i__] = s * work[indlpl + i__] * l[i__];
        s = work[inds + i__] - *lambda;
        /* L50: */
    }
    sawnan1 = sisnan_(&s);
    if (sawnan1)
    {
        goto L60;
    }
    i__1 = r2 - 1;
    for (i__ = r1;
            i__ <= i__1;
            ++i__)
    {
        dplus = d__[i__] + s;
        work[indlpl + i__] = ld[i__] / dplus;
        work[inds + i__] = s * work[indlpl + i__] * l[i__];
        s = work[inds + i__] - *lambda;
        /* L51: */
    }
    sawnan1 = sisnan_(&s);
L60:
    if (sawnan1)
    {
        /* Runs a slower version of the above loop if a NaN is detected */
        neg1 = 0;
        s = work[inds + *b1 - 1] - *lambda;
        i__1 = r1 - 1;
        for (i__ = *b1;
                i__ <= i__1;
                ++i__)
        {
            dplus = d__[i__] + s;
            if (f2c_abs(dplus) < *pivmin)
            {
                dplus = -(*pivmin);
            }
            work[indlpl + i__] = ld[i__] / dplus;
            if (dplus < 0.f)
            {
                ++neg1;
            }
            work[inds + i__] = s * work[indlpl + i__] * l[i__];
            if (work[indlpl + i__] == 0.f)
            {
                work[inds + i__] = lld[i__];
            }
            s = work[inds + i__] - *lambda;
            /* L70: */
        }
        i__1 = r2 - 1;
        for (i__ = r1;
                i__ <= i__1;
                ++i__)
        {
            dplus = d__[i__] + s;
            if (f2c_abs(dplus) < *pivmin)
            {
                dplus = -(*pivmin);
            }
            work[indlpl + i__] = ld[i__] / dplus;
            work[inds + i__] = s * work[indlpl + i__] * l[i__];
            if (work[indlpl + i__] == 0.f)
            {
                work[inds + i__] = lld[i__];
            }
            s = work[inds + i__] - *lambda;
            /* L71: */
        }
    }
    /* Compute the progressive transform (using the differential form) */
    /* until the index R1 */
    sawnan2 = FALSE_;
    neg2 = 0;
    work[indp + *bn - 1] = d__[*bn] - *lambda;
    i__1 = r1;
    for (i__ = *bn - 1;
            i__ >= i__1;
            --i__)
    {
        dminus = lld[i__] + work[indp + i__];
        tmp = d__[i__] / dminus;
        if (dminus < 0.f)
        {
            ++neg2;
        }
        work[indumn + i__] = l[i__] * tmp;
        work[indp + i__ - 1] = work[indp + i__] * tmp - *lambda;
        /* L80: */
    }
    tmp = work[indp + r1 - 1];
    sawnan2 = sisnan_(&tmp);
    if (sawnan2)
    {
        /* Runs a slower version of the above loop if a NaN is detected */
        neg2 = 0;
        i__1 = r1;
        for (i__ = *bn - 1;
                i__ >= i__1;
                --i__)
        {
            dminus = lld[i__] + work[indp + i__];
            if (f2c_abs(dminus) < *pivmin)
            {
                dminus = -(*pivmin);
            }
            tmp = d__[i__] / dminus;
            if (dminus < 0.f)
            {
                ++neg2;
            }
            work[indumn + i__] = l[i__] * tmp;
            work[indp + i__ - 1] = work[indp + i__] * tmp - *lambda;
            if (tmp == 0.f)
            {
                work[indp + i__ - 1] = d__[i__] - *lambda;
            }
            /* L100: */
        }
    }
    /* Find the index (from R1 to R2) of the largest (in magnitude) */
    /* diagonal element of the inverse */
    *mingma = work[inds + r1 - 1] + work[indp + r1 - 1];
    if (*mingma < 0.f)
    {
        ++neg1;
    }
    if (*wantnc)
    {
        *negcnt = neg1 + neg2;
    }
    else
    {
        *negcnt = -1;
    }
    if (f2c_abs(*mingma) == 0.f)
    {
        *mingma = eps * work[inds + r1 - 1];
    }
    *r__ = r1;
    i__1 = r2 - 1;
    for (i__ = r1;
            i__ <= i__1;
            ++i__)
    {
        tmp = work[inds + i__] + work[indp + i__];
        if (tmp == 0.f)
        {
            tmp = eps * work[inds + i__];
        }
        if (f2c_abs(tmp) <= f2c_abs(*mingma))
        {
            *mingma = tmp;
            *r__ = i__ + 1;
        }
        /* L110: */
    }
    /* Compute the FP vector: solve N^T v = e_r */
    isuppz[1] = *b1;
    isuppz[2] = *bn;
    i__1 = *r__;
    z__[i__1].r = 1.f;
    z__[i__1].i = 0.f; // , expr subst
    *ztz = 1.f;
    /* Compute the FP vector upwards from R */
    if (! sawnan1 && ! sawnan2)
    {
        i__1 = *b1;
        for (i__ = *r__ - 1;
                i__ >= i__1;
                --i__)
        {
            i__2 = i__;
            i__3 = indlpl + i__;
            i__4 = i__ + 1;
            q__2.r = work[i__3] * z__[i__4].r;
            q__2.i = work[i__3] * z__[i__4] .i; // , expr subst
            q__1.r = -q__2.r;
            q__1.i = -q__2.i; // , expr subst
            z__[i__2].r = q__1.r;
            z__[i__2].i = q__1.i; // , expr subst
            if ((c_abs(&z__[i__]) + c_abs(&z__[i__ + 1])) * (r__1 = ld[i__], f2c_abs(r__1)) < *gaptol)
            {
                i__2 = i__;
                z__[i__2].r = 0.f;
                z__[i__2].i = 0.f; // , expr subst
                isuppz[1] = i__ + 1;
                goto L220;
            }
            i__2 = i__;
            i__3 = i__;
            q__1.r = z__[i__2].r * z__[i__3].r - z__[i__2].i * z__[i__3].i;
            q__1.i = z__[i__2].r * z__[i__3].i + z__[i__2].i * z__[ i__3].r; // , expr subst
            *ztz += q__1.r;
            /* L210: */
        }
L220:
        ;
    }
    else
    {
        /* Run slower loop if NaN occurred. */
        i__1 = *b1;
        for (i__ = *r__ - 1;
                i__ >= i__1;
                --i__)
        {
            i__2 = i__ + 1;
            if (z__[i__2].r == 0.f && z__[i__2].i == 0.f)
            {
                i__2 = i__;
                r__1 = -(ld[i__ + 1] / ld[i__]);
                i__3 = i__ + 2;
                q__1.r = r__1 * z__[i__3].r;
                q__1.i = r__1 * z__[i__3].i; // , expr subst
                z__[i__2].r = q__1.r;
                z__[i__2].i = q__1.i; // , expr subst
            }
            else
            {
                i__2 = i__;
                i__3 = indlpl + i__;
                i__4 = i__ + 1;
                q__2.r = work[i__3] * z__[i__4].r;
                q__2.i = work[i__3] * z__[ i__4].i; // , expr subst
                q__1.r = -q__2.r;
                q__1.i = -q__2.i; // , expr subst
                z__[i__2].r = q__1.r;
                z__[i__2].i = q__1.i; // , expr subst
            }
            if ((c_abs(&z__[i__]) + c_abs(&z__[i__ + 1])) * (r__1 = ld[i__], f2c_abs(r__1)) < *gaptol)
            {
                i__2 = i__;
                z__[i__2].r = 0.f;
                z__[i__2].i = 0.f; // , expr subst
                isuppz[1] = i__ + 1;
                goto L240;
            }
            i__2 = i__;
            i__3 = i__;
            q__1.r = z__[i__2].r * z__[i__3].r - z__[i__2].i * z__[i__3].i;
            q__1.i = z__[i__2].r * z__[i__3].i + z__[i__2].i * z__[ i__3].r; // , expr subst
            *ztz += q__1.r;
            /* L230: */
        }
L240:
        ;
    }
    /* Compute the FP vector downwards from R in blocks of size BLKSIZ */
    if (! sawnan1 && ! sawnan2)
    {
        i__1 = *bn - 1;
        for (i__ = *r__;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__ + 1;
            i__3 = indumn + i__;
            i__4 = i__;
            q__2.r = work[i__3] * z__[i__4].r;
            q__2.i = work[i__3] * z__[i__4] .i; // , expr subst
            q__1.r = -q__2.r;
            q__1.i = -q__2.i; // , expr subst
            z__[i__2].r = q__1.r;
            z__[i__2].i = q__1.i; // , expr subst
            if ((c_abs(&z__[i__]) + c_abs(&z__[i__ + 1])) * (r__1 = ld[i__], f2c_abs(r__1)) < *gaptol)
            {
                i__2 = i__ + 1;
                z__[i__2].r = 0.f;
                z__[i__2].i = 0.f; // , expr subst
                isuppz[2] = i__;
                goto L260;
            }
            i__2 = i__ + 1;
            i__3 = i__ + 1;
            q__1.r = z__[i__2].r * z__[i__3].r - z__[i__2].i * z__[i__3].i;
            q__1.i = z__[i__2].r * z__[i__3].i + z__[i__2].i * z__[ i__3].r; // , expr subst
            *ztz += q__1.r;
            /* L250: */
        }
L260:
        ;
    }
    else
    {
        /* Run slower loop if NaN occurred. */
        i__1 = *bn - 1;
        for (i__ = *r__;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__;
            if (z__[i__2].r == 0.f && z__[i__2].i == 0.f)
            {
                i__2 = i__ + 1;
                r__1 = -(ld[i__ - 1] / ld[i__]);
                i__3 = i__ - 1;
                q__1.r = r__1 * z__[i__3].r;
                q__1.i = r__1 * z__[i__3].i; // , expr subst
                z__[i__2].r = q__1.r;
                z__[i__2].i = q__1.i; // , expr subst
            }
            else
            {
                i__2 = i__ + 1;
                i__3 = indumn + i__;
                i__4 = i__;
                q__2.r = work[i__3] * z__[i__4].r;
                q__2.i = work[i__3] * z__[ i__4].i; // , expr subst
                q__1.r = -q__2.r;
                q__1.i = -q__2.i; // , expr subst
                z__[i__2].r = q__1.r;
                z__[i__2].i = q__1.i; // , expr subst
            }
            if ((c_abs(&z__[i__]) + c_abs(&z__[i__ + 1])) * (r__1 = ld[i__], f2c_abs(r__1)) < *gaptol)
            {
                i__2 = i__ + 1;
                z__[i__2].r = 0.f;
                z__[i__2].i = 0.f; // , expr subst
                isuppz[2] = i__;
                goto L280;
            }
            i__2 = i__ + 1;
            i__3 = i__ + 1;
            q__1.r = z__[i__2].r * z__[i__3].r - z__[i__2].i * z__[i__3].i;
            q__1.i = z__[i__2].r * z__[i__3].i + z__[i__2].i * z__[ i__3].r; // , expr subst
            *ztz += q__1.r;
            /* L270: */
        }
L280:
        ;
    }
    /* Compute quantities for convergence test */
    tmp = 1.f / *ztz;
    *nrminv = sqrt(tmp);
    *resid = f2c_abs(*mingma) * *nrminv;
    *rqcorr = *mingma * tmp;
    return 0;
    /* End of CLAR1V */
}
/* clar1v_ */
