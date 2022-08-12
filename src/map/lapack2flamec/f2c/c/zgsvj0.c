/* ../netlib/v3.9.0/zgsvj0.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b27 = 1.;
/* > \brief <b> ZGSVJ0 pre-processor for the routine zgesvj. </b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGSVJ0 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgsvj0. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgsvj0. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgsvj0. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS, */
/* SFMIN, TOL, NSWEEP, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP */
/* DOUBLE PRECISION EPS, SFMIN, TOL */
/* CHARACTER*1 JOBV */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), D( N ), V( LDV, * ), WORK( LWORK ) */
/* DOUBLE PRECISION SVA( N ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGSVJ0 is called from ZGESVJ as a pre-processor and that is its main */
/* > purpose. It applies Jacobi rotations in the same way as ZGESVJ does, but */
/* > it does not check convergence (stopping criterion). Few tuning */
/* > parameters (marked by [TP]) are available for the implementer. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBV */
/* > \verbatim */
/* > JOBV is CHARACTER*1 */
/* > Specifies whether the output from this procedure is used */
/* > to compute the matrix V: */
/* > = 'V': the product of the Jacobi rotations is accumulated */
/* > by postmulyiplying the N-by-N array V. */
/* > (See the description of V.) */
/* > = 'A': the product of the Jacobi rotations is accumulated */
/* > by postmulyiplying the MV-by-N array V. */
/* > (See the descriptions of MV and V.) */
/* > = 'N': the Jacobi rotations are not accumulated. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the input matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the input matrix A. */
/* > M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, M-by-N matrix A, such that A*diag(D) represents */
/* > the input matrix. */
/* > On exit, */
/* > A_onexit * diag(D_onexit) represents the input matrix A*diag(D) */
/* > post-multiplied by a sequence of Jacobi rotations, where the */
/* > rotation threshold and the total number of sweeps are given in */
/* > TOL and NSWEEP, respectively. */
/* > (See the descriptions of D, TOL and NSWEEP.) */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is COMPLEX*16 array, dimension (N) */
/* > The array D accumulates the scaling factors from the complex scaled */
/* > Jacobi rotations. */
/* > On entry, A*diag(D) represents the input matrix. */
/* > On exit, A_onexit*diag(D_onexit) represents the input matrix */
/* > post-multiplied by a sequence of Jacobi rotations, where the */
/* > rotation threshold and the total number of sweeps are given in */
/* > TOL and NSWEEP, respectively. */
/* > (See the descriptions of A, TOL and NSWEEP.) */
/* > \endverbatim */
/* > */
/* > \param[in,out] SVA */
/* > \verbatim */
/* > SVA is DOUBLE PRECISION array, dimension (N) */
/* > On entry, SVA contains the Euclidean norms of the columns of */
/* > the matrix A*diag(D). */
/* > On exit, SVA contains the Euclidean norms of the columns of */
/* > the matrix A_onexit*diag(D_onexit). */
/* > \endverbatim */
/* > */
/* > \param[in] MV */
/* > \verbatim */
/* > MV is INTEGER */
/* > If JOBV = 'A', then MV rows of V are post-multipled by a */
/* > sequence of Jacobi rotations. */
/* > If JOBV = 'N', then MV is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* > V is COMPLEX*16 array, dimension (LDV,N) */
/* > If JOBV = 'V' then N rows of V are post-multipled by a */
/* > sequence of Jacobi rotations. */
/* > If JOBV = 'A' then MV rows of V are post-multipled by a */
/* > sequence of Jacobi rotations. */
/* > If JOBV = 'N', then V is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V, LDV >= 1. */
/* > If JOBV = 'V', LDV >= N. */
/* > If JOBV = 'A', LDV >= MV. */
/* > \endverbatim */
/* > */
/* > \param[in] EPS */
/* > \verbatim */
/* > EPS is DOUBLE PRECISION */
/* > EPS = DLAMCH('Epsilon') */
/* > \endverbatim */
/* > */
/* > \param[in] SFMIN */
/* > \verbatim */
/* > SFMIN is DOUBLE PRECISION */
/* > SFMIN = DLAMCH('Safe Minimum') */
/* > \endverbatim */
/* > */
/* > \param[in] TOL */
/* > \verbatim */
/* > TOL is DOUBLE PRECISION */
/* > TOL is the threshold for Jacobi rotations. For a pair */
/* > A(:,p), A(:,q) of pivot columns, the Jacobi rotation is */
/* > applied only if ABS(COS(angle(A(:,p),A(:,q)))) > TOL. */
/* > \endverbatim */
/* > */
/* > \param[in] NSWEEP */
/* > \verbatim */
/* > NSWEEP is INTEGER */
/* > NSWEEP is the number of sweeps of Jacobi rotations to be */
/* > performed. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > LWORK is the dimension of WORK. LWORK >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, then the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2016 */
/* > \ingroup complex16OTHERcomputational */
/* > */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > ZGSVJ0 is used just to enable ZGESVJ to call a simplified version of */
/* > itself to work on a submatrix of the original matrix. */
/* > */
/* > Contributor: */
/* ============= */
/* > */
/* > Zlatko Drmac (Zagreb, Croatia) */
/* > */
/* > \par Bugs, Examples and Comments: */
/* ============================ */
/* > */
/* > Please report all bugs and send interesting test examples and comments to */
/* > drmac@math.hr. Thank you. */
/* ===================================================================== */
/* Subroutine */
int zgsvj0_(char *jobv, integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *d__, doublereal *sva, integer *mv, doublecomplex *v, integer *ldv, doublereal *eps, doublereal *sfmin, doublereal *tol, integer *nsweep, doublecomplex * work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgsvj0 inputs: jobv %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", mv %" FLA_IS ", ldv %" FLA_IS ", nsweep %" FLA_IS "",*jobv, *m, *n, *lda, *mv, *ldv, *nsweep);
    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    double sqrt(doublereal), z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    double d_sign(doublereal *, doublereal *);
    /* Local variables */
    doublereal bigtheta;
    integer pskipped, i__, p, q;
    doublereal t, rootsfmin, cs, sn;
    integer ir1, jbc;
    doublereal big;
    integer kbl, igl, ibr, jgl, nbl, mvl;
    doublereal aapp;
    doublecomplex aapq;
    doublereal aaqq;
    integer ierr;
    doublecomplex ompq;
    extern /* Subroutine */
    int zrot_(integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, doublecomplex *);
    doublereal aapp0, aapq1, temp1, apoaq, aqoap;
    extern logical lsame_(char *, char *);
    doublereal theta, small;
    logical applv, rsvec;
    extern /* Double Complex */
    VOID zdotc_f2c_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    logical rotok;
    extern /* Subroutine */
    int zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    integer ijblsk, swband, blskip;
    doublereal mxaapq;
    extern /* Subroutine */
    int zlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublecomplex *, integer *, integer *);
    doublereal thsign, mxsinj;
    extern /* Subroutine */
    int zlassq_(integer *, doublecomplex *, integer *, doublereal *, doublereal *);
    integer emptsw, notrot, iswrot, lkahead;
    doublereal rootbig, rooteps;
    integer rowskip;
    doublereal roottol;
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2016 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* from BLAS */
    /* from LAPACK */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --sva;
    --d__;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --work;
    /* Function Body */
    applv = lsame_(jobv, "A");
    rsvec = lsame_(jobv, "V");
    if (! (rsvec || applv || lsame_(jobv, "N")))
    {
        *info = -1;
    }
    else if (*m < 0)
    {
        *info = -2;
    }
    else if (*n < 0 || *n > *m)
    {
        *info = -3;
    }
    else if (*lda < *m)
    {
        *info = -5;
    }
    else if ((rsvec || applv) && *mv < 0)
    {
        *info = -8;
    }
    else if (rsvec && *ldv < *n || applv && *ldv < *mv)
    {
        *info = -10;
    }
    else if (*tol <= *eps)
    {
        *info = -13;
    }
    else if (*nsweep < 0)
    {
        *info = -14;
    }
    else if (*lwork < *m)
    {
        *info = -16;
    }
    else
    {
        *info = 0;
    }
    /* #:( */
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGSVJ0", &i__1);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (rsvec)
    {
        mvl = *n;
    }
    else if (applv)
    {
        mvl = *mv;
    }
    rsvec = rsvec || applv;
    rooteps = sqrt(*eps);
    rootsfmin = sqrt(*sfmin);
    small = *sfmin / *eps;
    big = 1. / *sfmin;
    rootbig = 1. / rootsfmin;
    bigtheta = 1. / rooteps;
    roottol = sqrt(*tol);
    /* .. Row-cyclic Jacobi SVD algorithm with column pivoting .. */
    emptsw = *n * (*n - 1) / 2;
    notrot = 0;
    /* .. Row-cyclic pivot strategy with de Rijk's pivoting .. */
    swband = 0;
    /* [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective */
    /* if ZGESVJ is used as a computational routine in the preconditioned */
    /* Jacobi SVD algorithm ZGEJSV. For sweeps i=1:SWBAND the procedure */
    /* works on pivots inside a band-like region around the diagonal. */
    /* The boundaries are determined dynamically, based on the number of */
    /* pivots above a threshold. */
    kbl = min(8,*n);
    /* [TP] KBL is a tuning parameter that defines the tile size in the */
    /* tiling of the p-q loops of pivot pairs. In general, an optimal */
    /* value of KBL depends on the matrix dimensions and on the */
    /* parameters of the computer's memory. */
    nbl = *n / kbl;
    if (nbl * kbl != *n)
    {
        ++nbl;
    }
    /* Computing 2nd power */
    i__1 = kbl;
    blskip = i__1 * i__1;
    /* [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL. */
    rowskip = min(5,kbl);
    /* [TP] ROWSKIP is a tuning parameter. */
    lkahead = 1;
    /* [TP] LKAHEAD is a tuning parameter. */
    /* Quasi block transformations, using the lower (upper) triangular */
    /* structure of the input matrix. The quasi-block-cycling usually */
    /* invokes cubic convergence. Big part of this cycle is done inside */
    /* canonical subspaces of dimensions less than M. */
    /* .. Row-cyclic pivot strategy with de Rijk's pivoting .. */
    i__1 = *nsweep;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        /* .. go go go ... */
        mxaapq = 0.;
        mxsinj = 0.;
        iswrot = 0;
        notrot = 0;
        pskipped = 0;
        /* Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs */
        /* 1 <= p < q <= N. This is the first step toward a blocked implementation */
        /* of the rotations. New implementation, based on block transformations, */
        /* is under development. */
        i__2 = nbl;
        for (ibr = 1;
                ibr <= i__2;
                ++ibr)
        {
            igl = (ibr - 1) * kbl + 1;
            /* Computing MIN */
            i__4 = lkahead;
            i__5 = nbl - ibr; // , expr subst
            i__3 = min(i__4,i__5);
            for (ir1 = 0;
                    ir1 <= i__3;
                    ++ir1)
            {
                igl += ir1 * kbl;
                /* Computing MIN */
                i__5 = igl + kbl - 1;
                i__6 = *n - 1; // , expr subst
                i__4 = min(i__5,i__6);
                for (p = igl;
                        p <= i__4;
                        ++p)
                {
                    /* .. de Rijk's pivoting */
                    i__5 = *n - p + 1;
                    q = idamax_(&i__5, &sva[p], &c__1) + p - 1;
                    if (p != q)
                    {
                        zswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
                        if (rsvec)
                        {
                            zswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &c__1);
                        }
                        temp1 = sva[p];
                        sva[p] = sva[q];
                        sva[q] = temp1;
                        i__5 = p;
                        aapq.r = d__[i__5].r;
                        aapq.i = d__[i__5].i; // , expr subst
                        i__5 = p;
                        i__6 = q;
                        d__[i__5].r = d__[i__6].r;
                        d__[i__5].i = d__[i__6].i; // , expr subst
                        i__5 = q;
                        d__[i__5].r = aapq.r;
                        d__[i__5].i = aapq.i; // , expr subst
                    }
                    if (ir1 == 0)
                    {
                        /* Column norms are periodically updated by explicit */
                        /* norm computation. */
                        /* Caveat: */
                        /* Unfortunately, some BLAS implementations compute SNCRM2(M,A(1,p),1) */
                        /* as SQRT(S=ZDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to */
                        /* overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to */
                        /* underflow for ||A(:,p)||_2 < SQRT(underflow_threshold). */
                        /* Hence, DZNRM2 cannot be trusted, not even in the case when */
                        /* the true norm is far from the under(over)flow boundaries. */
                        /* If properly implemented DZNRM2 is available, the IF-THEN-ELSE-END IF */
                        /* below should be replaced with "AAPP = DZNRM2( M, A(1,p), 1 )". */
                        if (sva[p] < rootbig && sva[p] > rootsfmin)
                        {
                            sva[p] = dznrm2_(m, &a[p * a_dim1 + 1], &c__1);
                        }
                        else
                        {
                            temp1 = 0.;
                            aapp = 1.;
                            zlassq_(m, &a[p * a_dim1 + 1], &c__1, &temp1, & aapp);
                            sva[p] = temp1 * sqrt(aapp);
                        }
                        aapp = sva[p];
                    }
                    else
                    {
                        aapp = sva[p];
                    }
                    if (aapp > 0.)
                    {
                        pskipped = 0;
                        /* Computing MIN */
                        i__6 = igl + kbl - 1;
                        i__5 = min(i__6,*n);
                        for (q = p + 1;
                                q <= i__5;
                                ++q)
                        {
                            aaqq = sva[q];
                            if (aaqq > 0.)
                            {
                                aapp0 = aapp;
                                if (aaqq >= 1.)
                                {
                                    rotok = small * aapp <= aaqq;
                                    if (aapp < big / aaqq)
                                    {
                                        zdotc_f2c_(&z__3, m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], & c__1);
                                        z__2.r = z__3.r / aaqq;
                                        z__2.i = z__3.i / aaqq; // , expr subst
                                        z__1.r = z__2.r / aapp;
                                        z__1.i = z__2.i / aapp; // , expr subst
                                        aapq.r = z__1.r;
                                        aapq.i = z__1.i; // , expr subst
                                    }
                                    else
                                    {
                                        zcopy_(m, &a[p * a_dim1 + 1], &c__1, & work[1], &c__1);
                                        zlascl_("G", &c__0, &c__0, &aapp, & c_b27, m, &c__1, &work[1], lda, &ierr);
                                        zdotc_f2c_(&z__2, m, &work[1], &c__1, &a[ q * a_dim1 + 1], &c__1);
                                        z__1.r = z__2.r / aaqq;
                                        z__1.i = z__2.i / aaqq; // , expr subst
                                        aapq.r = z__1.r;
                                        aapq.i = z__1.i; // , expr subst
                                    }
                                }
                                else
                                {
                                    rotok = aapp <= aaqq / small;
                                    if (aapp > small / aaqq)
                                    {
                                        zdotc_f2c_(&z__3, m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], & c__1);
                                        z__2.r = z__3.r / aapp;
                                        z__2.i = z__3.i / aapp; // , expr subst
                                        z__1.r = z__2.r / aaqq;
                                        z__1.i = z__2.i / aaqq; // , expr subst
                                        aapq.r = z__1.r;
                                        aapq.i = z__1.i; // , expr subst
                                    }
                                    else
                                    {
                                        zcopy_(m, &a[q * a_dim1 + 1], &c__1, & work[1], &c__1);
                                        zlascl_("G", &c__0, &c__0, &aaqq, & c_b27, m, &c__1, &work[1], lda, &ierr);
                                        zdotc_f2c_(&z__2, m, &a[p * a_dim1 + 1], & c__1, &work[1], &c__1);
                                        z__1.r = z__2.r / aapp;
                                        z__1.i = z__2.i / aapp; // , expr subst
                                        aapq.r = z__1.r;
                                        aapq.i = z__1.i; // , expr subst
                                    }
                                }
                                /* AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q) */
                                aapq1 = -z_abs(&aapq);
                                /* Computing MAX */
                                d__1 = mxaapq;
                                d__2 = -aapq1; // , expr subst
                                mxaapq = max(d__1,d__2);
                                /* TO rotate or NOT to rotate, THAT is the question ... */
                                if (f2c_dabs(aapq1) > *tol)
                                {
                                    d__1 = z_abs(&aapq);
                                    z__1.r = aapq.r / d__1;
                                    z__1.i = aapq.i / d__1; // , expr subst
                                    ompq.r = z__1.r;
                                    ompq.i = z__1.i; // , expr subst
                                    /* .. rotate */
                                    /* [RTD] ROTATED = ROTATED + ONE */
                                    if (ir1 == 0)
                                    {
                                        notrot = 0;
                                        pskipped = 0;
                                        ++iswrot;
                                    }
                                    if (rotok)
                                    {
                                        aqoap = aaqq / aapp;
                                        apoaq = aapp / aaqq;
                                        theta = (d__1 = aqoap - apoaq, f2c_dabs( d__1)) * -.5 / aapq1;
                                        if (f2c_dabs(theta) > bigtheta)
                                        {
                                            t = .5 / theta;
                                            cs = 1.;
                                            d_cnjg(&z__2, &ompq);
                                            z__1.r = t * z__2.r;
                                            z__1.i = t * z__2.i; // , expr subst
                                            zrot_(m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], &c__1, &cs, &z__1);
                                            if (rsvec)
                                            {
                                                d_cnjg(&z__2, &ompq);
                                                z__1.r = t * z__2.r;
                                                z__1.i = t * z__2.i; // , expr subst
                                                zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &c__1, &cs, &z__1);
                                            }
                                            /* Computing MAX */
                                            d__1 = 0.;
                                            d__2 = t * apoaq * aapq1 + 1.; // , expr subst
                                            sva[q] = aaqq * sqrt((max(d__1, d__2)));
                                            /* Computing MAX */
                                            d__1 = 0.;
                                            d__2 = 1. - t * aqoap * aapq1; // , expr subst
                                            aapp *= sqrt((max(d__1,d__2)));
                                            /* Computing MAX */
                                            d__1 = mxsinj;
                                            d__2 = f2c_dabs(t); // , expr subst
                                            mxsinj = max(d__1,d__2);
                                        }
                                        else
                                        {
                                            /* .. choose correct signum for THETA and rotate */
                                            thsign = -d_sign(&c_b27, &aapq1);
                                            t = 1. / (theta + thsign * sqrt( theta * theta + 1.));
                                            cs = sqrt(1. / (t * t + 1.));
                                            sn = t * cs;
                                            /* Computing MAX */
                                            d__1 = mxsinj;
                                            d__2 = f2c_dabs(sn); // , expr subst
                                            mxsinj = max(d__1,d__2);
                                            /* Computing MAX */
                                            d__1 = 0.;
                                            d__2 = t * apoaq * aapq1 + 1.; // , expr subst
                                            sva[q] = aaqq * sqrt((max(d__1, d__2)));
                                            /* Computing MAX */
                                            d__1 = 0.;
                                            d__2 = 1. - t * aqoap * aapq1; // , expr subst
                                            aapp *= sqrt((max(d__1,d__2)));
                                            d_cnjg(&z__2, &ompq);
                                            z__1.r = sn * z__2.r;
                                            z__1.i = sn * z__2.i; // , expr subst
                                            zrot_(m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], &c__1, &cs, &z__1);
                                            if (rsvec)
                                            {
                                                d_cnjg(&z__2, &ompq);
                                                z__1.r = sn * z__2.r;
                                                z__1.i = sn * z__2.i; // , expr subst
                                                zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &c__1, &cs, &z__1);
                                            }
                                        }
                                        i__6 = p;
                                        i__7 = q;
                                        z__2.r = -d__[i__7].r;
                                        z__2.i = -d__[ i__7].i; // , expr subst
                                        z__1.r = z__2.r * ompq.r - z__2.i * ompq.i;
                                        z__1.i = z__2.r * ompq.i + z__2.i * ompq.r; // , expr subst
                                        d__[i__6].r = z__1.r;
                                        d__[i__6].i = z__1.i; // , expr subst
                                    }
                                    else
                                    {
                                        /* .. have to use modified Gram-Schmidt like transformation */
                                        zcopy_(m, &a[p * a_dim1 + 1], &c__1, & work[1], &c__1);
                                        zlascl_("G", &c__0, &c__0, &aapp, & c_b27, m, &c__1, &work[1], lda, &ierr);
                                        zlascl_("G", &c__0, &c__0, &aaqq, & c_b27, m, &c__1, &a[q * a_dim1 + 1], lda, &ierr);
                                        z__1.r = -aapq.r;
                                        z__1.i = -aapq.i; // , expr subst
                                        zaxpy_(m, &z__1, &work[1], &c__1, &a[ q * a_dim1 + 1], &c__1);
                                        zlascl_("G", &c__0, &c__0, &c_b27, & aaqq, m, &c__1, &a[q * a_dim1 + 1], lda, &ierr);
                                        /* Computing MAX */
                                        d__1 = 0.;
                                        d__2 = 1. - aapq1 * aapq1; // , expr subst
                                        sva[q] = aaqq * sqrt((max(d__1,d__2))) ;
                                        mxsinj = max(mxsinj,*sfmin);
                                    }
                                    /* END IF ROTOK THEN ... ELSE */
                                    /* In the case of cancellation in updating SVA(q), SVA(p) */
                                    /* recompute SVA(q), SVA(p). */
                                    /* Computing 2nd power */
                                    d__1 = sva[q] / aaqq;
                                    if (d__1 * d__1 <= rooteps)
                                    {
                                        if (aaqq < rootbig && aaqq > rootsfmin)
                                        {
                                            sva[q] = dznrm2_(m, &a[q * a_dim1 + 1], &c__1);
                                        }
                                        else
                                        {
                                            t = 0.;
                                            aaqq = 1.;
                                            zlassq_(m, &a[q * a_dim1 + 1], & c__1, &t, &aaqq);
                                            sva[q] = t * sqrt(aaqq);
                                        }
                                    }
                                    if (aapp / aapp0 <= rooteps)
                                    {
                                        if (aapp < rootbig && aapp > rootsfmin)
                                        {
                                            aapp = dznrm2_(m, &a[p * a_dim1 + 1], &c__1);
                                        }
                                        else
                                        {
                                            t = 0.;
                                            aapp = 1.;
                                            zlassq_(m, &a[p * a_dim1 + 1], & c__1, &t, &aapp);
                                            aapp = t * sqrt(aapp);
                                        }
                                        sva[p] = aapp;
                                    }
                                }
                                else
                                {
                                    /* A(:,p) and A(:,q) already numerically orthogonal */
                                    if (ir1 == 0)
                                    {
                                        ++notrot;
                                    }
                                    /* [RTD] SKIPPED = SKIPPED + 1 */
                                    ++pskipped;
                                }
                            }
                            else
                            {
                                /* A(:,q) is zero column */
                                if (ir1 == 0)
                                {
                                    ++notrot;
                                }
                                ++pskipped;
                            }
                            if (i__ <= swband && pskipped > rowskip)
                            {
                                if (ir1 == 0)
                                {
                                    aapp = -aapp;
                                }
                                notrot = 0;
                                goto L2103;
                            }
                            /* L2002: */
                        }
                        /* END q-LOOP */
L2103: /* bailed out of q-loop */
                        sva[p] = aapp;
                    }
                    else
                    {
                        sva[p] = aapp;
                        if (ir1 == 0 && aapp == 0.)
                        {
                            /* Computing MIN */
                            i__5 = igl + kbl - 1;
                            notrot = notrot + min(i__5,*n) - p;
                        }
                    }
                    /* L2001: */
                }
                /* end of the p-loop */
                /* end of doing the block ( ibr, ibr ) */
                /* L1002: */
            }
            /* end of ir1-loop */
            /* ... go to the off diagonal blocks */
            igl = (ibr - 1) * kbl + 1;
            i__3 = nbl;
            for (jbc = ibr + 1;
                    jbc <= i__3;
                    ++jbc)
            {
                jgl = (jbc - 1) * kbl + 1;
                /* doing the block at ( ibr, jbc ) */
                ijblsk = 0;
                /* Computing MIN */
                i__5 = igl + kbl - 1;
                i__4 = min(i__5,*n);
                for (p = igl;
                        p <= i__4;
                        ++p)
                {
                    aapp = sva[p];
                    if (aapp > 0.)
                    {
                        pskipped = 0;
                        /* Computing MIN */
                        i__6 = jgl + kbl - 1;
                        i__5 = min(i__6,*n);
                        for (q = jgl;
                                q <= i__5;
                                ++q)
                        {
                            aaqq = sva[q];
                            if (aaqq > 0.)
                            {
                                aapp0 = aapp;
                                /* .. M x 2 Jacobi SVD .. */
                                /* Safe Gram matrix computation */
                                if (aaqq >= 1.)
                                {
                                    if (aapp >= aaqq)
                                    {
                                        rotok = small * aapp <= aaqq;
                                    }
                                    else
                                    {
                                        rotok = small * aaqq <= aapp;
                                    }
                                    if (aapp < big / aaqq)
                                    {
                                        zdotc_f2c_(&z__3, m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], & c__1);
                                        z__2.r = z__3.r / aaqq;
                                        z__2.i = z__3.i / aaqq; // , expr subst
                                        z__1.r = z__2.r / aapp;
                                        z__1.i = z__2.i / aapp; // , expr subst
                                        aapq.r = z__1.r;
                                        aapq.i = z__1.i; // , expr subst
                                    }
                                    else
                                    {
                                        zcopy_(m, &a[p * a_dim1 + 1], &c__1, & work[1], &c__1);
                                        zlascl_("G", &c__0, &c__0, &aapp, & c_b27, m, &c__1, &work[1], lda, &ierr);
                                        zdotc_f2c_(&z__2, m, &work[1], &c__1, &a[ q * a_dim1 + 1], &c__1);
                                        z__1.r = z__2.r / aaqq;
                                        z__1.i = z__2.i / aaqq; // , expr subst
                                        aapq.r = z__1.r;
                                        aapq.i = z__1.i; // , expr subst
                                    }
                                }
                                else
                                {
                                    if (aapp >= aaqq)
                                    {
                                        rotok = aapp <= aaqq / small;
                                    }
                                    else
                                    {
                                        rotok = aaqq <= aapp / small;
                                    }
                                    if (aapp > small / aaqq)
                                    {
                                        zdotc_f2c_(&z__3, m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], & c__1);
                                        d__1 = max(aaqq,aapp);
                                        z__2.r = z__3.r / d__1;
                                        z__2.i = z__3.i / d__1; // , expr subst
                                        d__2 = min(aaqq,aapp);
                                        z__1.r = z__2.r / d__2;
                                        z__1.i = z__2.i / d__2; // , expr subst
                                        aapq.r = z__1.r;
                                        aapq.i = z__1.i; // , expr subst
                                    }
                                    else
                                    {
                                        zcopy_(m, &a[q * a_dim1 + 1], &c__1, & work[1], &c__1);
                                        zlascl_("G", &c__0, &c__0, &aaqq, & c_b27, m, &c__1, &work[1], lda, &ierr);
                                        zdotc_f2c_(&z__2, m, &a[p * a_dim1 + 1], & c__1, &work[1], &c__1);
                                        z__1.r = z__2.r / aapp;
                                        z__1.i = z__2.i / aapp; // , expr subst
                                        aapq.r = z__1.r;
                                        aapq.i = z__1.i; // , expr subst
                                    }
                                }
                                /* AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q) */
                                aapq1 = -z_abs(&aapq);
                                /* Computing MAX */
                                d__1 = mxaapq;
                                d__2 = -aapq1; // , expr subst
                                mxaapq = max(d__1,d__2);
                                /* TO rotate or NOT to rotate, THAT is the question ... */
                                if (f2c_dabs(aapq1) > *tol)
                                {
                                    d__1 = z_abs(&aapq);
                                    z__1.r = aapq.r / d__1;
                                    z__1.i = aapq.i / d__1; // , expr subst
                                    ompq.r = z__1.r;
                                    ompq.i = z__1.i; // , expr subst
                                    notrot = 0;
                                    /* [RTD] ROTATED = ROTATED + 1 */
                                    pskipped = 0;
                                    ++iswrot;
                                    if (rotok)
                                    {
                                        aqoap = aaqq / aapp;
                                        apoaq = aapp / aaqq;
                                        theta = (d__1 = aqoap - apoaq, f2c_dabs( d__1)) * -.5 / aapq1;
                                        if (aaqq > aapp0)
                                        {
                                            theta = -theta;
                                        }
                                        if (f2c_dabs(theta) > bigtheta)
                                        {
                                            t = .5 / theta;
                                            cs = 1.;
                                            d_cnjg(&z__2, &ompq);
                                            z__1.r = t * z__2.r;
                                            z__1.i = t * z__2.i; // , expr subst
                                            zrot_(m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], &c__1, &cs, &z__1);
                                            if (rsvec)
                                            {
                                                d_cnjg(&z__2, &ompq);
                                                z__1.r = t * z__2.r;
                                                z__1.i = t * z__2.i; // , expr subst
                                                zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &c__1, &cs, &z__1);
                                            }
                                            /* Computing MAX */
                                            d__1 = 0.;
                                            d__2 = t * apoaq * aapq1 + 1.; // , expr subst
                                            sva[q] = aaqq * sqrt((max(d__1, d__2)));
                                            /* Computing MAX */
                                            d__1 = 0.;
                                            d__2 = 1. - t * aqoap * aapq1; // , expr subst
                                            aapp *= sqrt((max(d__1,d__2)));
                                            /* Computing MAX */
                                            d__1 = mxsinj;
                                            d__2 = f2c_dabs(t); // , expr subst
                                            mxsinj = max(d__1,d__2);
                                        }
                                        else
                                        {
                                            /* .. choose correct signum for THETA and rotate */
                                            thsign = -d_sign(&c_b27, &aapq1);
                                            if (aaqq > aapp0)
                                            {
                                                thsign = -thsign;
                                            }
                                            t = 1. / (theta + thsign * sqrt( theta * theta + 1.));
                                            cs = sqrt(1. / (t * t + 1.));
                                            sn = t * cs;
                                            /* Computing MAX */
                                            d__1 = mxsinj;
                                            d__2 = f2c_dabs(sn); // , expr subst
                                            mxsinj = max(d__1,d__2);
                                            /* Computing MAX */
                                            d__1 = 0.;
                                            d__2 = t * apoaq * aapq1 + 1.; // , expr subst
                                            sva[q] = aaqq * sqrt((max(d__1, d__2)));
                                            /* Computing MAX */
                                            d__1 = 0.;
                                            d__2 = 1. - t * aqoap * aapq1; // , expr subst
                                            aapp *= sqrt((max(d__1,d__2)));
                                            d_cnjg(&z__2, &ompq);
                                            z__1.r = sn * z__2.r;
                                            z__1.i = sn * z__2.i; // , expr subst
                                            zrot_(m, &a[p * a_dim1 + 1], & c__1, &a[q * a_dim1 + 1], &c__1, &cs, &z__1);
                                            if (rsvec)
                                            {
                                                d_cnjg(&z__2, &ompq);
                                                z__1.r = sn * z__2.r;
                                                z__1.i = sn * z__2.i; // , expr subst
                                                zrot_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], &c__1, &cs, &z__1);
                                            }
                                        }
                                        i__6 = p;
                                        i__7 = q;
                                        z__2.r = -d__[i__7].r;
                                        z__2.i = -d__[ i__7].i; // , expr subst
                                        z__1.r = z__2.r * ompq.r - z__2.i * ompq.i;
                                        z__1.i = z__2.r * ompq.i + z__2.i * ompq.r; // , expr subst
                                        d__[i__6].r = z__1.r;
                                        d__[i__6].i = z__1.i; // , expr subst
                                    }
                                    else
                                    {
                                        /* .. have to use modified Gram-Schmidt like transformation */
                                        if (aapp > aaqq)
                                        {
                                            zcopy_(m, &a[p * a_dim1 + 1], & c__1, &work[1], &c__1);
                                            zlascl_("G", &c__0, &c__0, &aapp, &c_b27, m, &c__1, &work[1], lda, &ierr);
                                            zlascl_("G", &c__0, &c__0, &aaqq, &c_b27, m, &c__1, &a[q * a_dim1 + 1], lda, &ierr);
                                            z__1.r = -aapq.r;
                                            z__1.i = -aapq.i; // , expr subst
                                            zaxpy_(m, &z__1, &work[1], &c__1, &a[q * a_dim1 + 1], &c__1) ;
                                            zlascl_("G", &c__0, &c__0, &c_b27, &aaqq, m, &c__1, &a[q * a_dim1 + 1], lda, &ierr);
                                            /* Computing MAX */
                                            d__1 = 0.;
                                            d__2 = 1. - aapq1 * aapq1; // , expr subst
                                            sva[q] = aaqq * sqrt((max(d__1, d__2)));
                                            mxsinj = max(mxsinj,*sfmin);
                                        }
                                        else
                                        {
                                            zcopy_(m, &a[q * a_dim1 + 1], & c__1, &work[1], &c__1);
                                            zlascl_("G", &c__0, &c__0, &aaqq, &c_b27, m, &c__1, &work[1], lda, &ierr);
                                            zlascl_("G", &c__0, &c__0, &aapp, &c_b27, m, &c__1, &a[p * a_dim1 + 1], lda, &ierr);
                                            d_cnjg(&z__2, &aapq);
                                            z__1.r = -z__2.r;
                                            z__1.i = -z__2.i; // , expr subst
                                            zaxpy_(m, &z__1, &work[1], &c__1, &a[p * a_dim1 + 1], &c__1) ;
                                            zlascl_("G", &c__0, &c__0, &c_b27, &aapp, m, &c__1, &a[p * a_dim1 + 1], lda, &ierr);
                                            /* Computing MAX */
                                            d__1 = 0.;
                                            d__2 = 1. - aapq1 * aapq1; // , expr subst
                                            sva[p] = aapp * sqrt((max(d__1, d__2)));
                                            mxsinj = max(mxsinj,*sfmin);
                                        }
                                    }
                                    /* END IF ROTOK THEN ... ELSE */
                                    /* In the case of cancellation in updating SVA(q), SVA(p) */
                                    /* .. recompute SVA(q), SVA(p) */
                                    /* Computing 2nd power */
                                    d__1 = sva[q] / aaqq;
                                    if (d__1 * d__1 <= rooteps)
                                    {
                                        if (aaqq < rootbig && aaqq > rootsfmin)
                                        {
                                            sva[q] = dznrm2_(m, &a[q * a_dim1 + 1], &c__1);
                                        }
                                        else
                                        {
                                            t = 0.;
                                            aaqq = 1.;
                                            zlassq_(m, &a[q * a_dim1 + 1], & c__1, &t, &aaqq);
                                            sva[q] = t * sqrt(aaqq);
                                        }
                                    }
                                    /* Computing 2nd power */
                                    d__1 = aapp / aapp0;
                                    if (d__1 * d__1 <= rooteps)
                                    {
                                        if (aapp < rootbig && aapp > rootsfmin)
                                        {
                                            aapp = dznrm2_(m, &a[p * a_dim1 + 1], &c__1);
                                        }
                                        else
                                        {
                                            t = 0.;
                                            aapp = 1.;
                                            zlassq_(m, &a[p * a_dim1 + 1], & c__1, &t, &aapp);
                                            aapp = t * sqrt(aapp);
                                        }
                                        sva[p] = aapp;
                                    }
                                    /* end of OK rotation */
                                }
                                else
                                {
                                    ++notrot;
                                    /* [RTD] SKIPPED = SKIPPED + 1 */
                                    ++pskipped;
                                    ++ijblsk;
                                }
                            }
                            else
                            {
                                ++notrot;
                                ++pskipped;
                                ++ijblsk;
                            }
                            if (i__ <= swband && ijblsk >= blskip)
                            {
                                sva[p] = aapp;
                                notrot = 0;
                                goto L2011;
                            }
                            if (i__ <= swband && pskipped > rowskip)
                            {
                                aapp = -aapp;
                                notrot = 0;
                                goto L2203;
                            }
                            /* L2200: */
                        }
                        /* end of the q-loop */
L2203:
                        sva[p] = aapp;
                    }
                    else
                    {
                        if (aapp == 0.)
                        {
                            /* Computing MIN */
                            i__5 = jgl + kbl - 1;
                            notrot = notrot + min(i__5,*n) - jgl + 1;
                        }
                        if (aapp < 0.)
                        {
                            notrot = 0;
                        }
                    }
                    /* L2100: */
                }
                /* end of the p-loop */
                /* L2010: */
            }
            /* end of the jbc-loop */
L2011: /* 2011 bailed out of the jbc-loop */
            /* Computing MIN */
            i__4 = igl + kbl - 1;
            i__3 = min(i__4,*n);
            for (p = igl;
                    p <= i__3;
                    ++p)
            {
                sva[p] = (d__1 = sva[p], f2c_dabs(d__1));
                /* L2012: */
            }
            /* ** */
            /* L2000: */
        }
        /* 2000 :: end of the ibr-loop */
        /* .. update SVA(N) */
        if (sva[*n] < rootbig && sva[*n] > rootsfmin)
        {
            sva[*n] = dznrm2_(m, &a[*n * a_dim1 + 1], &c__1);
        }
        else
        {
            t = 0.;
            aapp = 1.;
            zlassq_(m, &a[*n * a_dim1 + 1], &c__1, &t, &aapp);
            sva[*n] = t * sqrt(aapp);
        }
        /* Additional steering devices */
        if (i__ < swband && (mxaapq <= roottol || iswrot <= *n))
        {
            swband = i__;
        }
        if (i__ > swband + 1 && mxaapq < sqrt((doublereal) (*n)) * *tol && ( doublereal) (*n) * mxaapq * mxsinj < *tol)
        {
            goto L1994;
        }
        if (notrot >= emptsw)
        {
            goto L1994;
        }
        /* L1993: */
    }
    /* end i=1:NSWEEP loop */
    /* #:( Reaching this point means that the procedure has not converged. */
    *info = *nsweep - 1;
    goto L1995;
L1994: /* #:) Reaching this point means numerical convergence after the i-th */
    /* sweep. */
    *info = 0;
    /* #:) INFO = 0 confirms successful iterations. */
L1995: /* Sort the vector SVA() of column norms. */
    i__1 = *n - 1;
    for (p = 1;
            p <= i__1;
            ++p)
    {
        i__2 = *n - p + 1;
        q = idamax_(&i__2, &sva[p], &c__1) + p - 1;
        if (p != q)
        {
            temp1 = sva[p];
            sva[p] = sva[q];
            sva[q] = temp1;
            i__2 = p;
            aapq.r = d__[i__2].r;
            aapq.i = d__[i__2].i; // , expr subst
            i__2 = p;
            i__3 = q;
            d__[i__2].r = d__[i__3].r;
            d__[i__2].i = d__[i__3].i; // , expr subst
            i__2 = q;
            d__[i__2].r = aapq.r;
            d__[i__2].i = aapq.i; // , expr subst
            zswap_(m, &a[p * a_dim1 + 1], &c__1, &a[q * a_dim1 + 1], &c__1);
            if (rsvec)
            {
                zswap_(&mvl, &v[p * v_dim1 + 1], &c__1, &v[q * v_dim1 + 1], & c__1);
            }
        }
        /* L5991: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* .. */
    /* .. END OF ZGSVJ0 */
    /* .. */
}
/* zgsvj0_ */

