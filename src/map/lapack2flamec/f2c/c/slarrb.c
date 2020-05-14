/* ../netlib/slarrb.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLARRB provides limited bisection to locate eigenvalues for more accuracy. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLARRB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrb. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrb. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrb. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLARRB( N, D, LLD, IFIRST, ILAST, RTOL1, */
/* RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK, */
/* PIVMIN, SPDIAM, TWIST, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER IFIRST, ILAST, INFO, N, OFFSET, TWIST */
/* REAL PIVMIN, RTOL1, RTOL2, SPDIAM */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL D( * ), LLD( * ), W( * ), */
/* $ WERR( * ), WGAP( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Given the relatively robust representation(RRR) L D L^T, SLARRB */
/* > does "limited" bisection to refine the eigenvalues of L D L^T, */
/* > W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial */
/* > guesses for these eigenvalues are input in W, the corresponding estimate */
/* > of the error in these guesses and their gaps are input in WERR */
/* > and WGAP, respectively. During bisection, intervals */
/* > [left, right] are maintained by storing their mid-points and */
/* > semi-widths in the arrays W and WERR respectively. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The N diagonal elements of the diagonal matrix D. */
/* > \endverbatim */
/* > */
/* > \param[in] LLD */
/* > \verbatim */
/* > LLD is REAL array, dimension (N-1) */
/* > The (N-1) elements L(i)*L(i)*D(i). */
/* > \endverbatim */
/* > */
/* > \param[in] IFIRST */
/* > \verbatim */
/* > IFIRST is INTEGER */
/* > The index of the first eigenvalue to be computed. */
/* > \endverbatim */
/* > */
/* > \param[in] ILAST */
/* > \verbatim */
/* > ILAST is INTEGER */
/* > The index of the last eigenvalue to be computed. */
/* > \endverbatim */
/* > */
/* > \param[in] RTOL1 */
/* > \verbatim */
/* > RTOL1 is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] RTOL2 */
/* > \verbatim */
/* > RTOL2 is REAL */
/* > Tolerance for the convergence of the bisection intervals. */
/* > An interval [LEFT,RIGHT] has converged if */
/* > RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) ) */
/* > where GAP is the (estimated) distance to the nearest */
/* > eigenvalue. */
/* > \endverbatim */
/* > */
/* > \param[in] OFFSET */
/* > \verbatim */
/* > OFFSET is INTEGER */
/* > Offset for the arrays W, WGAP and WERR, i.e., the IFIRST-OFFSET */
/* > through ILAST-OFFSET elements of these arrays are to be used. */
/* > \endverbatim */
/* > */
/* > \param[in,out] W */
/* > \verbatim */
/* > W is REAL array, dimension (N) */
/* > On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are */
/* > estimates of the eigenvalues of L D L^T indexed IFIRST throug */
/* > ILAST. */
/* > On output, these estimates are refined. */
/* > \endverbatim */
/* > */
/* > \param[in,out] WGAP */
/* > \verbatim */
/* > WGAP is REAL array, dimension (N-1) */
/* > On input, the (estimated) gaps between consecutive */
/* > eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap between */
/* > eigenvalues I and I+1. Note that if IFIRST.EQ.ILAST */
/* > then WGAP(IFIRST-OFFSET) must be set to ZERO. */
/* > On output, these gaps are refined. */
/* > \endverbatim */
/* > */
/* > \param[in,out] WERR */
/* > \verbatim */
/* > WERR is REAL array, dimension (N) */
/* > On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are */
/* > the errors in the estimates of the corresponding elements in W. */
/* > On output, these errors are refined. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (2*N) */
/* > Workspace. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (2*N) */
/* > Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] PIVMIN */
/* > \verbatim */
/* > PIVMIN is REAL */
/* > The minimum pivot in the Sturm sequence. */
/* > \endverbatim */
/* > */
/* > \param[in] SPDIAM */
/* > \verbatim */
/* > SPDIAM is REAL */
/* > The spectral diameter of the matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] TWIST */
/* > \verbatim */
/* > TWIST is INTEGER */
/* > The twist index for the twisted factorization that is used */
/* > for the negcount. */
/* > TWIST = N: Compute negcount from L D L^T - LAMBDA I = L+ D+ L+^T */
/* > TWIST = 1: Compute negcount from L D L^T - LAMBDA I = U- D- U-^T */
/* > TWIST = R: Compute negcount from L D L^T - LAMBDA I = N(r) D(r) N(r) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > Error flag. */
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
int slarrb_(integer *n, real *d__, real *lld, integer * ifirst, integer *ilast, real *rtol1, real *rtol2, integer *offset, real *w, real *wgap, real *werr, real *work, integer *iwork, real * pivmin, real *spdiam, integer *twist, integer *info)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Builtin functions */
    double log(doublereal);
    /* Local variables */
    integer i__, k, r__, i1, ii, ip;
    real gap, mid, tmp, back, lgap, rgap, left;
    integer iter, nint, prev, next;
    real cvrgd, right, width;
    extern integer slaneg_(integer *, real *, real *, real *, real *, integer *);
    integer negcnt;
    real mnwdth;
    integer olnint, maxitr;
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
    --iwork;
    --work;
    --werr;
    --wgap;
    --w;
    --lld;
    --d__;
    /* Function Body */
    *info = 0;
    maxitr = (integer) ((log(*spdiam + *pivmin) - log(*pivmin)) / log(2.f)) + 2;
    mnwdth = *pivmin * 2.f;
    r__ = *twist;
    if (r__ < 1 || r__ > *n)
    {
        r__ = *n;
    }
    /* Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ]. */
    /* The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while */
    /* Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 ) */
    /* for an unconverged interval is set to the index of the next unconverged */
    /* interval, and is -1 or 0 for a converged interval. Thus a linked */
    /* list of unconverged intervals is set up. */
    i1 = *ifirst;
    /* The number of unconverged intervals */
    nint = 0;
    /* The last unconverged interval found */
    prev = 0;
    rgap = wgap[i1 - *offset];
    i__1 = *ilast;
    for (i__ = i1;
            i__ <= i__1;
            ++i__)
    {
        k = i__ << 1;
        ii = i__ - *offset;
        left = w[ii] - werr[ii];
        right = w[ii] + werr[ii];
        lgap = rgap;
        rgap = wgap[ii];
        gap = min(lgap,rgap);
        /* Make sure that [LEFT,RIGHT] contains the desired eigenvalue */
        /* Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT */
        /* Do while( NEGCNT(LEFT).GT.I-1 ) */
        back = werr[ii];
L20:
        negcnt = slaneg_(n, &d__[1], &lld[1], &left, pivmin, &r__);
        if (negcnt > i__ - 1)
        {
            left -= back;
            back *= 2.f;
            goto L20;
        }
        /* Do while( NEGCNT(RIGHT).LT.I ) */
        /* Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT */
        back = werr[ii];
L50:
        negcnt = slaneg_(n, &d__[1], &lld[1], &right, pivmin, &r__);
        if (negcnt < i__)
        {
            right += back;
            back *= 2.f;
            goto L50;
        }
        width = (r__1 = left - right, f2c_abs(r__1)) * .5f;
        /* Computing MAX */
        r__1 = f2c_abs(left);
        r__2 = f2c_abs(right); // , expr subst
        tmp = max(r__1,r__2);
        /* Computing MAX */
        r__1 = *rtol1 * gap;
        r__2 = *rtol2 * tmp; // , expr subst
        cvrgd = max(r__1,r__2);
        if (width <= cvrgd || width <= mnwdth)
        {
            /* This interval has already converged and does not need refinement. */
            /* (Note that the gaps might change through refining the */
            /* eigenvalues, however, they can only get bigger.) */
            /* Remove it from the list. */
            iwork[k - 1] = -1;
            /* Make sure that I1 always points to the first unconverged interval */
            if (i__ == i1 && i__ < *ilast)
            {
                i1 = i__ + 1;
            }
            if (prev >= i1 && i__ <= *ilast)
            {
                iwork[(prev << 1) - 1] = i__ + 1;
            }
        }
        else
        {
            /* unconverged interval found */
            prev = i__;
            ++nint;
            iwork[k - 1] = i__ + 1;
            iwork[k] = negcnt;
        }
        work[k - 1] = left;
        work[k] = right;
        /* L75: */
    }
    /* Do while( NINT.GT.0 ), i.e. there are still unconverged intervals */
    /* and while (ITER.LT.MAXITR) */
    iter = 0;
L80:
    prev = i1 - 1;
    i__ = i1;
    olnint = nint;
    i__1 = olnint;
    for (ip = 1;
            ip <= i__1;
            ++ip)
    {
        k = i__ << 1;
        ii = i__ - *offset;
        rgap = wgap[ii];
        lgap = rgap;
        if (ii > 1)
        {
            lgap = wgap[ii - 1];
        }
        gap = min(lgap,rgap);
        next = iwork[k - 1];
        left = work[k - 1];
        right = work[k];
        mid = (left + right) * .5f;
        /* semiwidth of interval */
        width = right - mid;
        /* Computing MAX */
        r__1 = f2c_abs(left);
        r__2 = f2c_abs(right); // , expr subst
        tmp = max(r__1,r__2);
        /* Computing MAX */
        r__1 = *rtol1 * gap;
        r__2 = *rtol2 * tmp; // , expr subst
        cvrgd = max(r__1,r__2);
        if (width <= cvrgd || width <= mnwdth || iter == maxitr)
        {
            /* reduce number of unconverged intervals */
            --nint;
            /* Mark interval as converged. */
            iwork[k - 1] = 0;
            if (i1 == i__)
            {
                i1 = next;
            }
            else
            {
                /* Prev holds the last unconverged interval previously examined */
                if (prev >= i1)
                {
                    iwork[(prev << 1) - 1] = next;
                }
            }
            i__ = next;
            goto L100;
        }
        prev = i__;
        /* Perform one bisection step */
        negcnt = slaneg_(n, &d__[1], &lld[1], &mid, pivmin, &r__);
        if (negcnt <= i__ - 1)
        {
            work[k - 1] = mid;
        }
        else
        {
            work[k] = mid;
        }
        i__ = next;
L100:
        ;
    }
    ++iter;
    /* do another loop if there are still unconverged intervals */
    /* However, in the last iteration, all intervals are accepted */
    /* since this is the best we can do. */
    if (nint > 0 && iter <= maxitr)
    {
        goto L80;
    }
    /* At this point, all the intervals have converged */
    i__1 = *ilast;
    for (i__ = *ifirst;
            i__ <= i__1;
            ++i__)
    {
        k = i__ << 1;
        ii = i__ - *offset;
        /* All intervals marked by '0' have been refined. */
        if (iwork[k - 1] == 0)
        {
            w[ii] = (work[k - 1] + work[k]) * .5f;
            werr[ii] = work[k] - w[ii];
        }
        /* L110: */
    }
    i__1 = *ilast;
    for (i__ = *ifirst + 1;
            i__ <= i__1;
            ++i__)
    {
        k = i__ << 1;
        ii = i__ - *offset;
        /* Computing MAX */
        r__1 = 0.f;
        r__2 = w[ii] - werr[ii] - w[ii - 1] - werr[ii - 1]; // , expr subst
        wgap[ii - 1] = max(r__1,r__2);
        /* L111: */
    }
    return 0;
    /* End of SLARRB */
}
/* slarrb_ */
