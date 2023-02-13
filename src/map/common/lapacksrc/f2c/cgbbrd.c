/* ../netlib/cgbbrd.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    0.f,0.f
}
;
static complex c_b2 =
{
    1.f,0.f
}
;
static integer c__1 = 1;
/* > \brief \b CGBBRD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGBBRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbbrd. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbbrd. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbbrd. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q, */
/* LDQ, PT, LDPT, C, LDC, WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER VECT */
/* INTEGER INFO, KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), E( * ), RWORK( * ) */
/* COMPLEX AB( LDAB, * ), C( LDC, * ), PT( LDPT, * ), */
/* $ Q( LDQ, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBBRD reduces a complex general m-by-n band matrix A to real upper */
/* > bidiagonal form B by a unitary transformation: Q**H * A * P = B. */
/* > */
/* > The routine computes B, and optionally forms Q or P**H, or computes */
/* > Q**H*C for a given matrix C. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] VECT */
/* > \verbatim */
/* > VECT is CHARACTER*1 */
/* > Specifies whether or not the matrices Q and P**H are to be */
/* > formed. */
/* > = 'N': do not form Q or P**H;
*/
/* > = 'Q': form Q only;
*/
/* > = 'P': form P**H only;
*/
/* > = 'B': form both. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NCC */
/* > \verbatim */
/* > NCC is INTEGER */
/* > The number of columns of the matrix C. NCC >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* > KL is INTEGER */
/* > The number of subdiagonals of the matrix A. KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* > KU is INTEGER */
/* > The number of superdiagonals of the matrix A. KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB,N) */
/* > On entry, the m-by-n band matrix A, stored in rows 1 to */
/* > KL+KU+1. The j-th column of A is stored in the j-th column of */
/* > the array AB as follows: */
/* > AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl). */
/* > On exit, A is overwritten by values generated during the */
/* > reduction. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array A. LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is REAL array, dimension (min(M,N)) */
/* > The diagonal elements of the bidiagonal matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is REAL array, dimension (min(M,N)-1) */
/* > The superdiagonal elements of the bidiagonal matrix B. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* > Q is COMPLEX array, dimension (LDQ,M) */
/* > If VECT = 'Q' or 'B', the m-by-m unitary matrix Q. */
/* > If VECT = 'N' or 'P', the array Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. */
/* > LDQ >= max(1,M) if VECT = 'Q' or 'B';
LDQ >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] PT */
/* > \verbatim */
/* > PT is COMPLEX array, dimension (LDPT,N) */
/* > If VECT = 'P' or 'B', the n-by-n unitary matrix P'. */
/* > If VECT = 'N' or 'Q', the array PT is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDPT */
/* > \verbatim */
/* > LDPT is INTEGER */
/* > The leading dimension of the array PT. */
/* > LDPT >= max(1,N) if VECT = 'P' or 'B';
LDPT >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is COMPLEX array, dimension (LDC,NCC) */
/* > On entry, an m-by-ncc matrix C. */
/* > On exit, C is overwritten by Q**H*C. */
/* > C is not referenced if NCC = 0. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. */
/* > LDC >= max(1,M) if NCC > 0;
LDC >= 1 if NCC = 0. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (max(M,N)) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (max(M,N)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexGBcomputational */
/* ===================================================================== */
/* Subroutine */
int cgbbrd_(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku, complex *ab, integer *ldab, real *d__, real *e, complex *q, integer *ldq, complex *pt, integer *ldpt, complex *c__, integer *ldc, complex *work, real *rwork, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, c_dim1, c_offset, pt_dim1, pt_offset, q_dim1, q_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    double c_abs(complex *);
    /* Local variables */
    integer i__, j, l;
    complex t;
    integer j1, j2, kb;
    complex ra, rb;
    real rc;
    integer kk, ml, nr, mu;
    complex rs;
    integer kb1, ml0, mu0, klm, kun, nrt, klu1, inca;
    real abst;
    extern /* Subroutine */
    int crot_(integer *, complex *, integer *, complex *, integer *, real *, complex *), cscal_(integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *);
    logical wantb, wantc;
    integer minmn;
    logical wantq;
    extern /* Subroutine */
    int claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *), clartg_(complex *, complex *, real *, complex *, complex *), xerbla_(char *, integer *), clargv_(integer *, complex *, integer *, complex *, integer *, real *, integer *), clartv_(integer *, complex *, integer *, complex *, integer *, real *, complex *, integer *);
    logical wantpt;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --d__;
    --e;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    pt_dim1 = *ldpt;
    pt_offset = 1 + pt_dim1;
    pt -= pt_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    --rwork;
    /* Function Body */
    wantb = lsame_(vect, "B");
    wantq = lsame_(vect, "Q") || wantb;
    wantpt = lsame_(vect, "P") || wantb;
    wantc = *ncc > 0;
    klu1 = *kl + *ku + 1;
    *info = 0;
    if (! wantq && ! wantpt && ! lsame_(vect, "N"))
    {
        *info = -1;
    }
    else if (*m < 0)
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*ncc < 0)
    {
        *info = -4;
    }
    else if (*kl < 0)
    {
        *info = -5;
    }
    else if (*ku < 0)
    {
        *info = -6;
    }
    else if (*ldab < klu1)
    {
        *info = -8;
    }
    else if (*ldq < 1 || wantq && *ldq < max(1,*m))
    {
        *info = -12;
    }
    else if (*ldpt < 1 || wantpt && *ldpt < max(1,*n))
    {
        *info = -14;
    }
    else if (*ldc < 1 || wantc && *ldc < max(1,*m))
    {
        *info = -16;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGBBRD", &i__1);
        return 0;
    }
    /* Initialize Q and P**H to the unit matrix, if needed */
    if (wantq)
    {
        claset_("Full", m, m, &c_b1, &c_b2, &q[q_offset], ldq);
    }
    if (wantpt)
    {
        claset_("Full", n, n, &c_b1, &c_b2, &pt[pt_offset], ldpt);
    }
    /* Quick return if possible. */
    if (*m == 0 || *n == 0)
    {
        return 0;
    }
    minmn = min(*m,*n);
    if (*kl + *ku > 1)
    {
        /* Reduce to upper bidiagonal form if KU > 0;
        if KU = 0, reduce */
        /* first to lower bidiagonal form and then transform to upper */
        /* bidiagonal */
        if (*ku > 0)
        {
            ml0 = 1;
            mu0 = 2;
        }
        else
        {
            ml0 = 2;
            mu0 = 1;
        }
        /* Wherever possible, plane rotations are generated and applied in */
        /* vector operations of length NR over the index set J1:J2:KLU1. */
        /* The complex sines of the plane rotations are stored in WORK, */
        /* and the real cosines in RWORK. */
        /* Computing MIN */
        i__1 = *m - 1;
        klm = min(i__1,*kl);
        /* Computing MIN */
        i__1 = *n - 1;
        kun = min(i__1,*ku);
        kb = klm + kun;
        kb1 = kb + 1;
        inca = kb1 * *ldab;
        nr = 0;
        j1 = klm + 2;
        j2 = 1 - kun;
        i__1 = minmn;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            /* Reduce i-th column and i-th row of matrix to bidiagonal form */
            ml = klm + 1;
            mu = kun + 1;
            i__2 = kb;
            for (kk = 1;
                    kk <= i__2;
                    ++kk)
            {
                j1 += kb;
                j2 += kb;
                /* generate plane rotations to annihilate nonzero elements */
                /* which have been created below the band */
                if (nr > 0)
                {
                    clargv_(&nr, &ab[klu1 + (j1 - klm - 1) * ab_dim1], &inca, &work[j1], &kb1, &rwork[j1], &kb1);
                }
                /* apply plane rotations from the left */
                i__3 = kb;
                for (l = 1;
                        l <= i__3;
                        ++l)
                {
                    if (j2 - klm + l - 1 > *n)
                    {
                        nrt = nr - 1;
                    }
                    else
                    {
                        nrt = nr;
                    }
                    if (nrt > 0)
                    {
                        clartv_(&nrt, &ab[klu1 - l + (j1 - klm + l - 1) * ab_dim1], &inca, &ab[klu1 - l + 1 + (j1 - klm + l - 1) * ab_dim1], &inca, &rwork[j1], &work[ j1], &kb1);
                    }
                    /* L10: */
                }
                if (ml > ml0)
                {
                    if (ml <= *m - i__ + 1)
                    {
                        /* generate plane rotation to annihilate a(i+ml-1,i) */
                        /* within the band, and apply rotation from the left */
                        clartg_(&ab[*ku + ml - 1 + i__ * ab_dim1], &ab[*ku + ml + i__ * ab_dim1], &rwork[i__ + ml - 1], & work[i__ + ml - 1], &ra);
                        i__3 = *ku + ml - 1 + i__ * ab_dim1;
                        ab[i__3].r = ra.r;
                        ab[i__3].i = ra.i; // , expr subst
                        if (i__ < *n)
                        {
                            /* Computing MIN */
                            i__4 = *ku + ml - 2;
                            i__5 = *n - i__; // , expr subst
                            i__3 = min(i__4,i__5);
                            i__6 = *ldab - 1;
                            i__7 = *ldab - 1;
                            crot_(&i__3, &ab[*ku + ml - 2 + (i__ + 1) * ab_dim1], &i__6, &ab[*ku + ml - 1 + (i__ + 1) * ab_dim1], &i__7, &rwork[i__ + ml - 1], &work[i__ + ml - 1]);
                        }
                    }
                    ++nr;
                    j1 -= kb1;
                }
                if (wantq)
                {
                    /* accumulate product of plane rotations in Q */
                    i__3 = j2;
                    i__4 = kb1;
                    for (j = j1;
                            i__4 < 0 ? j >= i__3 : j <= i__3;
                            j += i__4)
                    {
                        r_cnjg(&q__1, &work[j]);
                        crot_(m, &q[(j - 1) * q_dim1 + 1], &c__1, &q[j * q_dim1 + 1], &c__1, &rwork[j], &q__1);
                        /* L20: */
                    }
                }
                if (wantc)
                {
                    /* apply plane rotations to C */
                    i__4 = j2;
                    i__3 = kb1;
                    for (j = j1;
                            i__3 < 0 ? j >= i__4 : j <= i__4;
                            j += i__3)
                    {
                        crot_(ncc, &c__[j - 1 + c_dim1], ldc, &c__[j + c_dim1] , ldc, &rwork[j], &work[j]);
                        /* L30: */
                    }
                }
                if (j2 + kun > *n)
                {
                    /* adjust J2 to keep within the bounds of the matrix */
                    --nr;
                    j2 -= kb1;
                }
                i__3 = j2;
                i__4 = kb1;
                for (j = j1;
                        i__4 < 0 ? j >= i__3 : j <= i__3;
                        j += i__4)
                {
                    /* create nonzero element a(j-1,j+ku) above the band */
                    /* and store it in WORK(n+1:2*n) */
                    i__5 = j + kun;
                    i__6 = j;
                    i__7 = (j + kun) * ab_dim1 + 1;
                    q__1.r = work[i__6].r * ab[i__7].r - work[i__6].i * ab[ i__7].i;
                    q__1.i = work[i__6].r * ab[i__7].i + work[i__6].i * ab[i__7].r; // , expr subst
                    work[i__5].r = q__1.r;
                    work[i__5].i = q__1.i; // , expr subst
                    i__5 = (j + kun) * ab_dim1 + 1;
                    i__6 = j;
                    i__7 = (j + kun) * ab_dim1 + 1;
                    q__1.r = rwork[i__6] * ab[i__7].r;
                    q__1.i = rwork[i__6] * ab[i__7].i; // , expr subst
                    ab[i__5].r = q__1.r;
                    ab[i__5].i = q__1.i; // , expr subst
                    /* L40: */
                }
                /* generate plane rotations to annihilate nonzero elements */
                /* which have been generated above the band */
                if (nr > 0)
                {
                    clargv_(&nr, &ab[(j1 + kun - 1) * ab_dim1 + 1], &inca, & work[j1 + kun], &kb1, &rwork[j1 + kun], &kb1);
                }
                /* apply plane rotations from the right */
                i__4 = kb;
                for (l = 1;
                        l <= i__4;
                        ++l)
                {
                    if (j2 + l - 1 > *m)
                    {
                        nrt = nr - 1;
                    }
                    else
                    {
                        nrt = nr;
                    }
                    if (nrt > 0)
                    {
                        clartv_(&nrt, &ab[l + 1 + (j1 + kun - 1) * ab_dim1], & inca, &ab[l + (j1 + kun) * ab_dim1], &inca, & rwork[j1 + kun], &work[j1 + kun], &kb1);
                    }
                    /* L50: */
                }
                if (ml == ml0 && mu > mu0)
                {
                    if (mu <= *n - i__ + 1)
                    {
                        /* generate plane rotation to annihilate a(i,i+mu-1) */
                        /* within the band, and apply rotation from the right */
                        clartg_(&ab[*ku - mu + 3 + (i__ + mu - 2) * ab_dim1], &ab[*ku - mu + 2 + (i__ + mu - 1) * ab_dim1], &rwork[i__ + mu - 1], &work[i__ + mu - 1], & ra);
                        i__4 = *ku - mu + 3 + (i__ + mu - 2) * ab_dim1;
                        ab[i__4].r = ra.r;
                        ab[i__4].i = ra.i; // , expr subst
                        /* Computing MIN */
                        i__3 = *kl + mu - 2;
                        i__5 = *m - i__; // , expr subst
                        i__4 = min(i__3,i__5);
                        crot_(&i__4, &ab[*ku - mu + 4 + (i__ + mu - 2) * ab_dim1], &c__1, &ab[*ku - mu + 3 + (i__ + mu - 1) * ab_dim1], &c__1, &rwork[i__ + mu - 1], &work[i__ + mu - 1]);
                    }
                    ++nr;
                    j1 -= kb1;
                }
                if (wantpt)
                {
                    /* accumulate product of plane rotations in P**H */
                    i__4 = j2;
                    i__3 = kb1;
                    for (j = j1;
                            i__3 < 0 ? j >= i__4 : j <= i__4;
                            j += i__3)
                    {
                        r_cnjg(&q__1, &work[j + kun]);
                        crot_(n, &pt[j + kun - 1 + pt_dim1], ldpt, &pt[j + kun + pt_dim1], ldpt, &rwork[j + kun], &q__1);
                        /* L60: */
                    }
                }
                if (j2 + kb > *m)
                {
                    /* adjust J2 to keep within the bounds of the matrix */
                    --nr;
                    j2 -= kb1;
                }
                i__3 = j2;
                i__4 = kb1;
                for (j = j1;
                        i__4 < 0 ? j >= i__3 : j <= i__3;
                        j += i__4)
                {
                    /* create nonzero element a(j+kl+ku,j+ku-1) below the */
                    /* band and store it in WORK(1:n) */
                    i__5 = j + kb;
                    i__6 = j + kun;
                    i__7 = klu1 + (j + kun) * ab_dim1;
                    q__1.r = work[i__6].r * ab[i__7].r - work[i__6].i * ab[ i__7].i;
                    q__1.i = work[i__6].r * ab[i__7].i + work[i__6].i * ab[i__7].r; // , expr subst
                    work[i__5].r = q__1.r;
                    work[i__5].i = q__1.i; // , expr subst
                    i__5 = klu1 + (j + kun) * ab_dim1;
                    i__6 = j + kun;
                    i__7 = klu1 + (j + kun) * ab_dim1;
                    q__1.r = rwork[i__6] * ab[i__7].r;
                    q__1.i = rwork[i__6] * ab[i__7].i; // , expr subst
                    ab[i__5].r = q__1.r;
                    ab[i__5].i = q__1.i; // , expr subst
                    /* L70: */
                }
                if (ml > ml0)
                {
                    --ml;
                }
                else
                {
                    --mu;
                }
                /* L80: */
            }
            /* L90: */
        }
    }
    if (*ku == 0 && *kl > 0)
    {
        /* A has been reduced to complex lower bidiagonal form */
        /* Transform lower bidiagonal form to upper bidiagonal by applying */
        /* plane rotations from the left, overwriting superdiagonal */
        /* elements on subdiagonal elements */
        /* Computing MIN */
        i__2 = *m - 1;
        i__1 = min(i__2,*n);
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            clartg_(&ab[i__ * ab_dim1 + 1], &ab[i__ * ab_dim1 + 2], &rc, &rs, &ra);
            i__2 = i__ * ab_dim1 + 1;
            ab[i__2].r = ra.r;
            ab[i__2].i = ra.i; // , expr subst
            if (i__ < *n)
            {
                i__2 = i__ * ab_dim1 + 2;
                i__4 = (i__ + 1) * ab_dim1 + 1;
                q__1.r = rs.r * ab[i__4].r - rs.i * ab[i__4].i;
                q__1.i = rs.r * ab[i__4].i + rs.i * ab[i__4].r; // , expr subst
                ab[i__2].r = q__1.r;
                ab[i__2].i = q__1.i; // , expr subst
                i__2 = (i__ + 1) * ab_dim1 + 1;
                i__4 = (i__ + 1) * ab_dim1 + 1;
                q__1.r = rc * ab[i__4].r;
                q__1.i = rc * ab[i__4].i; // , expr subst
                ab[i__2].r = q__1.r;
                ab[i__2].i = q__1.i; // , expr subst
            }
            if (wantq)
            {
                r_cnjg(&q__1, &rs);
                crot_(m, &q[i__ * q_dim1 + 1], &c__1, &q[(i__ + 1) * q_dim1 + 1], &c__1, &rc, &q__1);
            }
            if (wantc)
            {
                crot_(ncc, &c__[i__ + c_dim1], ldc, &c__[i__ + 1 + c_dim1], ldc, &rc, &rs);
            }
            /* L100: */
        }
    }
    else
    {
        /* A has been reduced to complex upper bidiagonal form or is */
        /* diagonal */
        if (*ku > 0 && *m < *n)
        {
            /* Annihilate a(m,m+1) by applying plane rotations from the */
            /* right */
            i__1 = *ku + (*m + 1) * ab_dim1;
            rb.r = ab[i__1].r;
            rb.i = ab[i__1].i; // , expr subst
            for (i__ = *m;
                    i__ >= 1;
                    --i__)
            {
                clartg_(&ab[*ku + 1 + i__ * ab_dim1], &rb, &rc, &rs, &ra);
                i__1 = *ku + 1 + i__ * ab_dim1;
                ab[i__1].r = ra.r;
                ab[i__1].i = ra.i; // , expr subst
                if (i__ > 1)
                {
                    r_cnjg(&q__3, &rs);
                    q__2.r = -q__3.r;
                    q__2.i = -q__3.i; // , expr subst
                    i__1 = *ku + i__ * ab_dim1;
                    q__1.r = q__2.r * ab[i__1].r - q__2.i * ab[i__1].i;
                    q__1.i = q__2.r * ab[i__1].i + q__2.i * ab[i__1] .r; // , expr subst
                    rb.r = q__1.r;
                    rb.i = q__1.i; // , expr subst
                    i__1 = *ku + i__ * ab_dim1;
                    i__2 = *ku + i__ * ab_dim1;
                    q__1.r = rc * ab[i__2].r;
                    q__1.i = rc * ab[i__2].i; // , expr subst
                    ab[i__1].r = q__1.r;
                    ab[i__1].i = q__1.i; // , expr subst
                }
                if (wantpt)
                {
                    r_cnjg(&q__1, &rs);
                    crot_(n, &pt[i__ + pt_dim1], ldpt, &pt[*m + 1 + pt_dim1], ldpt, &rc, &q__1);
                }
                /* L110: */
            }
        }
    }
    /* Make diagonal and superdiagonal elements real, storing them in D */
    /* and E */
    i__1 = *ku + 1 + ab_dim1;
    t.r = ab[i__1].r;
    t.i = ab[i__1].i; // , expr subst
    i__1 = minmn;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        abst = c_abs(&t);
        d__[i__] = abst;
        if (abst != 0.f)
        {
            q__1.r = t.r / abst;
            q__1.i = t.i / abst; // , expr subst
            t.r = q__1.r;
            t.i = q__1.i; // , expr subst
        }
        else
        {
            t.r = 1.f;
            t.i = 0.f; // , expr subst
        }
        if (wantq)
        {
            cscal_(m, &t, &q[i__ * q_dim1 + 1], &c__1);
        }
        if (wantc)
        {
            r_cnjg(&q__1, &t);
            cscal_(ncc, &q__1, &c__[i__ + c_dim1], ldc);
        }
        if (i__ < minmn)
        {
            if (*ku == 0 && *kl == 0)
            {
                e[i__] = 0.f;
                i__2 = (i__ + 1) * ab_dim1 + 1;
                t.r = ab[i__2].r;
                t.i = ab[i__2].i; // , expr subst
            }
            else
            {
                if (*ku == 0)
                {
                    i__2 = i__ * ab_dim1 + 2;
                    r_cnjg(&q__2, &t);
                    q__1.r = ab[i__2].r * q__2.r - ab[i__2].i * q__2.i;
                    q__1.i = ab[i__2].r * q__2.i + ab[i__2].i * q__2.r; // , expr subst
                    t.r = q__1.r;
                    t.i = q__1.i; // , expr subst
                }
                else
                {
                    i__2 = *ku + (i__ + 1) * ab_dim1;
                    r_cnjg(&q__2, &t);
                    q__1.r = ab[i__2].r * q__2.r - ab[i__2].i * q__2.i;
                    q__1.i = ab[i__2].r * q__2.i + ab[i__2].i * q__2.r; // , expr subst
                    t.r = q__1.r;
                    t.i = q__1.i; // , expr subst
                }
                abst = c_abs(&t);
                e[i__] = abst;
                if (abst != 0.f)
                {
                    q__1.r = t.r / abst;
                    q__1.i = t.i / abst; // , expr subst
                    t.r = q__1.r;
                    t.i = q__1.i; // , expr subst
                }
                else
                {
                    t.r = 1.f;
                    t.i = 0.f; // , expr subst
                }
                if (wantpt)
                {
                    cscal_(n, &t, &pt[i__ + 1 + pt_dim1], ldpt);
                }
                i__2 = *ku + 1 + (i__ + 1) * ab_dim1;
                r_cnjg(&q__2, &t);
                q__1.r = ab[i__2].r * q__2.r - ab[i__2].i * q__2.i;
                q__1.i = ab[i__2].r * q__2.i + ab[i__2].i * q__2.r; // , expr subst
                t.r = q__1.r;
                t.i = q__1.i; // , expr subst
            }
        }
        /* L120: */
    }
    return 0;
    /* End of CGBBRD */
}
/* cgbbrd_ */
