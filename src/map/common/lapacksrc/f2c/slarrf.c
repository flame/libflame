/* ../netlib/slarrf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SLARRF finds a new relatively robust representation such that at least one of the eigenvalues i s relatively isolated. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLARRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLARRF( N, D, L, LD, CLSTRT, CLEND, */
/* W, WGAP, WERR, */
/* SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA, */
/* DPLUS, LPLUS, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER CLSTRT, CLEND, INFO, N */
/* REAL CLGAPL, CLGAPR, PIVMIN, SIGMA, SPDIAM */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), DPLUS( * ), L( * ), LD( * ), */
/* $ LPLUS( * ), W( * ), WGAP( * ), WERR( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given the initial representation L D L^T and its cluster of close */
/* > eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ... */
/* > W( CLEND ), SLARRF finds a new relatively robust representation */
/* > L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the */
/* > eigenvalues of L(+) D(+) L(+)^T is relatively isolated. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix (subblock, if the matrix splitted). */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The N diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is REAL array, dimension (N-1) */
/* > The (N-1) subdiagonal elements of the unit bidiagonal */
/* > matrix L. */
/* > \endverbatim */
/* > */
/* > \param[in] LD */
/* > \verbatim */
/* > LD is REAL array, dimension (N-1) */
/* > The (N-1) elements L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] CLSTRT */
/* > \verbatim */
/* > CLSTRT is INTEGER */
/* > The index of the first eigenvalue in the cluster. */
/* > \endverbatim */
/* > */
/* > \param[in] CLEND */
/* > \verbatim */
/* > CLEND is INTEGER */
/* > The index of the last eigenvalue in the cluster. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* > W is REAL array, dimension */
/* > dimension is >= (CLEND-CLSTRT+1) */
/* > The eigenvalue APPROXIMATIONS of L D L^T in ascending order. */
/* > W( CLSTRT ) through W( CLEND ) form the cluster of relatively */
/* > close eigenalues. */
/* > \endverbatim */
/* > */
/* > \param[in,out] WGAP */
/* > \verbatim */
/* > WGAP is REAL array, dimension */
/* > dimension is >= (CLEND-CLSTRT+1) */
/* > The separation from the right neighbor eigenvalue in W. */
/* > \endverbatim */
/* > */
/* > \param[in] WERR */
/* > \verbatim */
/* > WERR is REAL array, dimension */
/* > dimension is >= (CLEND-CLSTRT+1) */
/* > WERR contain the semiwidth of the uncertainty */
/* > interval of the corresponding eigenvalue APPROXIMATION in W */
/* > \endverbatim */
/* > */
/* > \param[in] SPDIAM */
/* > \verbatim */
/* > SPDIAM is REAL */
/* > estimate of the spectral diameter obtained from the */
/* > Gerschgorin intervals */
/* > \endverbatim */
/* > */
/* > \param[in] CLGAPL */
/* > \verbatim */
/* > CLGAPL is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] CLGAPR */
/* > \verbatim */
/* > CLGAPR is REAL */
/* > absolute gap on each end of the cluster. */
/* > Set by the calling routine to protect against shifts too close */
/* > to eigenvalues outside the cluster. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* > PIVMIN is REAL */
/* > The minimum pivot allowed in the Sturm sequence. */
/* > \endverbatim */
/* > */
/* > \param[out] SIGMA */
/* > \verbatim */
/* > SIGMA is REAL */
/* > The shift used to form L(+) D(+) L(+)^T. */
/* > \endverbatim */
/* > */
/* > \param[out] DPLUS */
/* > \verbatim */
/* > DPLUS is REAL array, dimension (N) */
/* > The N diagonal elements of the diagonal matrix D(+). */
/* > \endverbatim */
/* > */
/* > \param[out] LPLUS */
/* > \verbatim */
/* > LPLUS is REAL array, dimension (N-1) */
/* > The first (N-1) elements of LPLUS contain the subdiagonal */
/* > elements of the unit bidiagonal matrix L(+). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (2*N) */
/* > Workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > Signals processing OK (=0) or failure (=1) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERauxiliary */
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
int slarrf_(integer *n, real *d__, real *l, real *ld, integer *clstrt, integer *clend, real *w, real *wgap, real *werr, real *spdiam, real *clgapl, real *clgapr, real *pivmin, real *sigma, real *dplus, real *lplus, real *work, integer *info)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__;
    real s, bestshift, smlgrowth, eps, tmp, max1, max2, rrr1, rrr2, znm2, growthbound, fail, fact, oldp;
    integer indx;
    real prod;
    integer ktry;
    real fail2, avgap, ldmax, rdmax;
    integer shift;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *);
    logical dorrr1;
    real ldelta;
    extern real slamch_(char *);
    logical nofail;
    real mingap, lsigma, rdelta;
    logical forcer;
    real rsigma, clwdth;
    extern logical sisnan_(real *);
    logical sawnan1, sawnan2, tryrrr1;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --work;
    --lplus;
    --dplus;
    --werr;
    --wgap;
    --w;
    --ld;
    --l;
    --d__;
    /* Function Body */
    *info = 0;
    fact = 2.f;
    eps = slamch_("Precision");
    shift = 0;
    forcer = FALSE_;
    /* Note that we cannot guarantee that for any of the shifts tried, */
    /* the factorization has a small or even moderate element growth. */
    /* There could be Ritz values at both ends of the cluster and despite */
    /* backing off, there are examples where all factorizations tried */
    /* (in IEEE mode, allowing zero pivots & infinities) have INFINITE */
    /* element growth. */
    /* For this reason, we should use PIVMIN in this subroutine so that at */
    /* least the L D L^T factorization exists. It can be checked afterwards */
    /* whether the element growth caused bad residuals/orthogonality. */
    /* Decide whether the code should accept the best among all */
    /* representations despite large element growth or signal INFO=1 */
    nofail = TRUE_;
    /* Compute the average gap length of the cluster */
    clwdth = (r__1 = w[*clend] - w[*clstrt], f2c_abs(r__1)) + werr[*clend] + werr[ *clstrt];
    avgap = clwdth / (real) (*clend - *clstrt);
    mingap = min(*clgapl,*clgapr);
    /* Initial values for shifts to both ends of cluster */
    /* Computing MIN */
    r__1 = w[*clstrt];
    r__2 = w[*clend]; // , expr subst
    lsigma = min(r__1,r__2) - werr[*clstrt];
    /* Computing MAX */
    r__1 = w[*clstrt];
    r__2 = w[*clend]; // , expr subst
    rsigma = max(r__1,r__2) + werr[*clend];
    /* Use a small fudge to make sure that we really shift to the outside */
    lsigma -= f2c_abs(lsigma) * 2.f * eps;
    rsigma += f2c_abs(rsigma) * 2.f * eps;
    /* Compute upper bounds for how much to back off the initial shifts */
    ldmax = mingap * .25f + *pivmin * 2.f;
    rdmax = mingap * .25f + *pivmin * 2.f;
    /* Computing MAX */
    r__1 = avgap;
    r__2 = wgap[*clstrt]; // , expr subst
    ldelta = max(r__1,r__2) / fact;
    /* Computing MAX */
    r__1 = avgap;
    r__2 = wgap[*clend - 1]; // , expr subst
    rdelta = max(r__1,r__2) / fact;
    /* Initialize the record of the best representation found */
    s = slamch_("S");
    smlgrowth = 1.f / s;
    fail = (real) (*n - 1) * mingap / (*spdiam * eps);
    fail2 = (real) (*n - 1) * mingap / (*spdiam * sqrt(eps));
    bestshift = lsigma;
    /* while (KTRY <= KTRYMAX) */
    ktry = 0;
    growthbound = *spdiam * 8.f;
L5:
    sawnan1 = FALSE_;
    sawnan2 = FALSE_;
    /* Ensure that we do not back off too much of the initial shifts */
    ldelta = min(ldmax,ldelta);
    rdelta = min(rdmax,rdelta);
    /* Compute the element growth when shifting to both ends of the cluster */
    /* accept the shift if there is no element growth at one of the two ends */
    /* Left end */
    s = -lsigma;
    dplus[1] = d__[1] + s;
    if (f2c_abs(dplus[1]) < *pivmin)
    {
        dplus[1] = -(*pivmin);
        /* Need to set SAWNAN1 because refined RRR test should not be used */
        /* in this case */
        sawnan1 = TRUE_;
    }
    max1 = f2c_abs(dplus[1]);
    i__1 = *n - 1;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        lplus[i__] = ld[i__] / dplus[i__];
        s = s * lplus[i__] * l[i__] - lsigma;
        dplus[i__ + 1] = d__[i__ + 1] + s;
        if ((r__1 = dplus[i__ + 1], f2c_abs(r__1)) < *pivmin)
        {
            dplus[i__ + 1] = -(*pivmin);
            /* Need to set SAWNAN1 because refined RRR test should not be used */
            /* in this case */
            sawnan1 = TRUE_;
        }
        /* Computing MAX */
        r__2 = max1;
        r__3 = (r__1 = dplus[i__ + 1], f2c_abs(r__1)); // , expr subst
        max1 = max(r__2,r__3);
        /* L6: */
    }
    sawnan1 = sawnan1 || sisnan_(&max1);
    if (forcer || max1 <= growthbound && ! sawnan1)
    {
        *sigma = lsigma;
        shift = 1;
        goto L100;
    }
    /* Right end */
    s = -rsigma;
    work[1] = d__[1] + s;
    if (f2c_abs(work[1]) < *pivmin)
    {
        work[1] = -(*pivmin);
        /* Need to set SAWNAN2 because refined RRR test should not be used */
        /* in this case */
        sawnan2 = TRUE_;
    }
    max2 = f2c_abs(work[1]);
    i__1 = *n - 1;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        work[*n + i__] = ld[i__] / work[i__];
        s = s * work[*n + i__] * l[i__] - rsigma;
        work[i__ + 1] = d__[i__ + 1] + s;
        if ((r__1 = work[i__ + 1], f2c_abs(r__1)) < *pivmin)
        {
            work[i__ + 1] = -(*pivmin);
            /* Need to set SAWNAN2 because refined RRR test should not be used */
            /* in this case */
            sawnan2 = TRUE_;
        }
        /* Computing MAX */
        r__2 = max2;
        r__3 = (r__1 = work[i__ + 1], f2c_abs(r__1)); // , expr subst
        max2 = max(r__2,r__3);
        /* L7: */
    }
    sawnan2 = sawnan2 || sisnan_(&max2);
    if (forcer || max2 <= growthbound && ! sawnan2)
    {
        *sigma = rsigma;
        shift = 2;
        goto L100;
    }
    /* If we are at this point, both shifts led to too much element growth */
    /* Record the better of the two shifts (provided it didn't lead to NaN) */
    if (sawnan1 && sawnan2)
    {
        /* both MAX1 and MAX2 are NaN */
        goto L50;
    }
    else
    {
        if (! sawnan1)
        {
            indx = 1;
            if (max1 <= smlgrowth)
            {
                smlgrowth = max1;
                bestshift = lsigma;
            }
        }
        if (! sawnan2)
        {
            if (sawnan1 || max2 <= max1)
            {
                indx = 2;
            }
            if (max2 <= smlgrowth)
            {
                smlgrowth = max2;
                bestshift = rsigma;
            }
        }
    }
    /* If we are here, both the left and the right shift led to */
    /* element growth. If the element growth is moderate, then */
    /* we may still accept the representation, if it passes a */
    /* refined test for RRR. This test supposes that no NaN occurred. */
    /* Moreover, we use the refined RRR test only for isolated clusters. */
    if (clwdth < mingap / 128.f && min(max1,max2) < fail2 && ! sawnan1 && ! sawnan2)
    {
        dorrr1 = TRUE_;
    }
    else
    {
        dorrr1 = FALSE_;
    }
    tryrrr1 = TRUE_;
    if (tryrrr1 && dorrr1)
    {
        if (indx == 1)
        {
            tmp = (r__1 = dplus[*n], f2c_abs(r__1));
            znm2 = 1.f;
            prod = 1.f;
            oldp = 1.f;
            for (i__ = *n - 1;
                    i__ >= 1;
                    --i__)
            {
                if (prod <= eps)
                {
                    prod = dplus[i__ + 1] * work[*n + i__ + 1] / (dplus[i__] * work[*n + i__]) * oldp;
                }
                else
                {
                    prod *= (r__1 = work[*n + i__], f2c_abs(r__1));
                }
                oldp = prod;
                /* Computing 2nd power */
                r__1 = prod;
                znm2 += r__1 * r__1;
                /* Computing MAX */
                r__2 = tmp;
                r__3 = (r__1 = dplus[i__] * prod, f2c_abs(r__1)); // , expr subst
                tmp = max(r__2,r__3);
                /* L15: */
            }
            rrr1 = tmp / (*spdiam * sqrt(znm2));
            if (rrr1 <= 8.f)
            {
                *sigma = lsigma;
                shift = 1;
                goto L100;
            }
        }
        else if (indx == 2)
        {
            tmp = (r__1 = work[*n], f2c_abs(r__1));
            znm2 = 1.f;
            prod = 1.f;
            oldp = 1.f;
            for (i__ = *n - 1;
                    i__ >= 1;
                    --i__)
            {
                if (prod <= eps)
                {
                    prod = work[i__ + 1] * lplus[i__ + 1] / (work[i__] * lplus[i__]) * oldp;
                }
                else
                {
                    prod *= (r__1 = lplus[i__], f2c_abs(r__1));
                }
                oldp = prod;
                /* Computing 2nd power */
                r__1 = prod;
                znm2 += r__1 * r__1;
                /* Computing MAX */
                r__2 = tmp;
                r__3 = (r__1 = work[i__] * prod, f2c_abs(r__1)); // , expr subst
                tmp = max(r__2,r__3);
                /* L16: */
            }
            rrr2 = tmp / (*spdiam * sqrt(znm2));
            if (rrr2 <= 8.f)
            {
                *sigma = rsigma;
                shift = 2;
                goto L100;
            }
        }
    }
L50:
    if (ktry < 1)
    {
        /* If we are here, both shifts failed also the RRR test. */
        /* Back off to the outside */
        /* Computing MAX */
        r__1 = lsigma - ldelta;
        r__2 = lsigma - ldmax; // , expr subst
        lsigma = max(r__1,r__2);
        /* Computing MIN */
        r__1 = rsigma + rdelta;
        r__2 = rsigma + rdmax; // , expr subst
        rsigma = min(r__1,r__2);
        ldelta *= 2.f;
        rdelta *= 2.f;
        ++ktry;
        goto L5;
    }
    else
    {
        /* None of the representations investigated satisfied our */
        /* criteria. Take the best one we found. */
        if (smlgrowth < fail || nofail)
        {
            lsigma = bestshift;
            rsigma = bestshift;
            forcer = TRUE_;
            goto L5;
        }
        else
        {
            *info = 1;
            return 0;
        }
    }
L100:
    if (shift == 1)
    {
    }
    else if (shift == 2)
    {
        /* store new L and D back into DPLUS, LPLUS */
        scopy_(n, &work[1], &c__1, &dplus[1], &c__1);
        i__1 = *n - 1;
        scopy_(&i__1, &work[*n + 1], &c__1, &lplus[1], &c__1);
    }
    return 0;
    /* End of SLARRF */
}
/* slarrf_ */
