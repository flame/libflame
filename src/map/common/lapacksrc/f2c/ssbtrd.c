/* ../netlib/ssbtrd.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b9 = 0.f;
static real c_b10 = 1.f;
static integer c__1 = 1;
/* > \brief \b SSBTRD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSBTRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbtrd. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbtrd. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbtrd. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ, */
/* WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO, VECT */
/* INTEGER INFO, KD, LDAB, LDQ, N */
/* .. */
/* .. Array Arguments .. */
/* REAL AB( LDAB, * ), D( * ), E( * ), Q( LDQ, * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSBTRD reduces a real symmetric band matrix A to symmetric */
/* > tridiagonal form T by an orthogonal similarity transformation: */
/* > Q**T * A * Q = T. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] VECT */
/* > \verbatim */
/* > VECT is CHARACTER*1 */
/* > = 'N': do not form Q;
*/
/* > = 'V': form Q;
*/
/* > = 'U': update a matrix X, by forming X*Q. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
*/
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is REAL array, dimension (LDAB,N) */
/* > On entry, the upper or lower triangle of the symmetric band */
/* > matrix A, stored in the first KD+1 rows of the array. The */
/* > j-th column of A is stored in the j-th column of the array AB */
/* > as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*/
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd). */
/* > On exit, the diagonal elements of AB are overwritten by the */
/* > diagonal elements of the tridiagonal matrix T;
if KD > 0, the */
/* > elements on the first superdiagonal (if UPLO = 'U') or the */
/* > first subdiagonal (if UPLO = 'L') are overwritten by the */
/* > off-diagonal elements of T;
the rest of AB is overwritten by */
/* > values generated during the reduction. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > The off-diagonal elements of the tridiagonal matrix T: */
/* > E(i) = T(i,i+1) if UPLO = 'U';
E(i) = T(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (LDQ,N) */
/* > On entry, if VECT = 'U', then Q must contain an N-by-N */
/* > matrix X;
if VECT = 'N' or 'V', then Q need not be set. */
/* > */
/* > On exit: */
/* > if VECT = 'V', Q contains the N-by-N orthogonal matrix Q;
*/
/* > if VECT = 'U', Q contains the product X*Q;
*/
/* > if VECT = 'N', the array Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. */
/* > LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Modified by Linda Kaufman, Bell Labs. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int ssbtrd_(char *vect, char *uplo, integer *n, integer *kd, real *ab, integer *ldab, real *d__, real *e, real *q, integer *ldq, real *work, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, q_dim1, q_offset, i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
    integer i__, j, k, l, i2, j1, j2, nq, nr, kd1, ibl, iqb, kdn, jin, nrt, kdm1, inca, jend, lend, jinc, incx, last;
    real temp;
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *);
    integer j1end, j1inc, iqend;
    extern logical lsame_(char *, char *);
    logical initq, wantq, upper;
    extern /* Subroutine */
    int slar2v_(integer *, real *, real *, real *, integer *, real *, real *, integer *);
    integer iqaend;
    extern /* Subroutine */
    int xerbla_(char *, integer *), slaset_( char *, integer *, integer *, real *, real *, real *, integer *), slartg_(real *, real *, real *, real *, real *), slargv_( integer *, real *, integer *, real *, integer *, real *, integer * ), slartv_(integer *, real *, integer *, real *, integer *, real * , real *, integer *);
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
    --work;
    /* Function Body */
    initq = lsame_(vect, "V");
    wantq = initq || lsame_(vect, "U");
    upper = lsame_(uplo, "U");
    kd1 = *kd + 1;
    kdm1 = *kd - 1;
    incx = *ldab - 1;
    iqend = 1;
    *info = 0;
    if (! wantq && ! lsame_(vect, "N"))
    {
        *info = -1;
    }
    else if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*kd < 0)
    {
        *info = -4;
    }
    else if (*ldab < kd1)
    {
        *info = -6;
    }
    else if (*ldq < max(1,*n) && wantq)
    {
        *info = -10;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSBTRD", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Initialize Q to the unit matrix, if needed */
    if (initq)
    {
        slaset_("Full", n, n, &c_b9, &c_b10, &q[q_offset], ldq);
    }
    /* Wherever possible, plane rotations are generated and applied in */
    /* vector operations of length NR over the index set J1:J2:KD1. */
    /* The cosines and sines of the plane rotations are stored in the */
    /* arrays D and WORK. */
    inca = kd1 * *ldab;
    /* Computing MIN */
    i__1 = *n - 1;
    kdn = min(i__1,*kd);
    if (upper)
    {
        if (*kd > 1)
        {
            /* Reduce to tridiagonal form, working with upper triangle */
            nr = 0;
            j1 = kdn + 2;
            j2 = 1;
            i__1 = *n - 2;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                /* Reduce i-th row of matrix to tridiagonal form */
                for (k = kdn + 1;
                        k >= 2;
                        --k)
                {
                    j1 += kdn;
                    j2 += kdn;
                    if (nr > 0)
                    {
                        /* generate plane rotations to annihilate nonzero */
                        /* elements which have been created outside the band */
                        slargv_(&nr, &ab[(j1 - 1) * ab_dim1 + 1], &inca, & work[j1], &kd1, &d__[j1], &kd1);
                        /* apply rotations from the right */
                        /* Dependent on the the number of diagonals either */
                        /* SLARTV or SROT is used */
                        if (nr >= (*kd << 1) - 1)
                        {
                            i__2 = *kd - 1;
                            for (l = 1;
                                    l <= i__2;
                                    ++l)
                            {
                                slartv_(&nr, &ab[l + 1 + (j1 - 1) * ab_dim1], &inca, &ab[l + j1 * ab_dim1], &inca, & d__[j1], &work[j1], &kd1);
                                /* L10: */
                            }
                        }
                        else
                        {
                            jend = j1 + (nr - 1) * kd1;
                            i__2 = jend;
                            i__3 = kd1;
                            for (jinc = j1;
                                    i__3 < 0 ? jinc >= i__2 : jinc <= i__2;
                                    jinc += i__3)
                            {
                                srot_(&kdm1, &ab[(jinc - 1) * ab_dim1 + 2], & c__1, &ab[jinc * ab_dim1 + 1], &c__1, &d__[jinc], &work[jinc]);
                                /* L20: */
                            }
                        }
                    }
                    if (k > 2)
                    {
                        if (k <= *n - i__ + 1)
                        {
                            /* generate plane rotation to annihilate a(i,i+k-1) */
                            /* within the band */
                            slartg_(&ab[*kd - k + 3 + (i__ + k - 2) * ab_dim1] , &ab[*kd - k + 2 + (i__ + k - 1) * ab_dim1], &d__[i__ + k - 1], &work[i__ + k - 1], &temp);
                            ab[*kd - k + 3 + (i__ + k - 2) * ab_dim1] = temp;
                            /* apply rotation from the right */
                            i__3 = k - 3;
                            srot_(&i__3, &ab[*kd - k + 4 + (i__ + k - 2) * ab_dim1], &c__1, &ab[*kd - k + 3 + (i__ + k - 1) * ab_dim1], &c__1, &d__[i__ + k - 1], &work[i__ + k - 1]);
                        }
                        ++nr;
                        j1 = j1 - kdn - 1;
                    }
                    /* apply plane rotations from both sides to diagonal */
                    /* blocks */
                    if (nr > 0)
                    {
                        slar2v_(&nr, &ab[kd1 + (j1 - 1) * ab_dim1], &ab[kd1 + j1 * ab_dim1], &ab[*kd + j1 * ab_dim1], &inca, &d__[j1], &work[j1], &kd1);
                    }
                    /* apply plane rotations from the left */
                    if (nr > 0)
                    {
                        if ((*kd << 1) - 1 < nr)
                        {
                            /* Dependent on the the number of diagonals either */
                            /* SLARTV or SROT is used */
                            i__3 = *kd - 1;
                            for (l = 1;
                                    l <= i__3;
                                    ++l)
                            {
                                if (j2 + l > *n)
                                {
                                    nrt = nr - 1;
                                }
                                else
                                {
                                    nrt = nr;
                                }
                                if (nrt > 0)
                                {
                                    slartv_(&nrt, &ab[*kd - l + (j1 + l) * ab_dim1], &inca, &ab[*kd - l + 1 + (j1 + l) * ab_dim1], &inca, & d__[j1], &work[j1], &kd1);
                                }
                                /* L30: */
                            }
                        }
                        else
                        {
                            j1end = j1 + kd1 * (nr - 2);
                            if (j1end >= j1)
                            {
                                i__3 = j1end;
                                i__2 = kd1;
                                for (jin = j1;
                                        i__2 < 0 ? jin >= i__3 : jin <= i__3;
                                        jin += i__2)
                                {
                                    i__4 = *kd - 1;
                                    srot_(&i__4, &ab[*kd - 1 + (jin + 1) * ab_dim1], &incx, &ab[*kd + (jin + 1) * ab_dim1], &incx, &d__[jin], & work[jin]);
                                    /* L40: */
                                }
                            }
                            /* Computing MIN */
                            i__2 = kdm1;
                            i__3 = *n - j2; // , expr subst
                            lend = min(i__2,i__3);
                            last = j1end + kd1;
                            if (lend > 0)
                            {
                                srot_(&lend, &ab[*kd - 1 + (last + 1) * ab_dim1], &incx, &ab[*kd + (last + 1) * ab_dim1], &incx, &d__[last], &work[ last]);
                            }
                        }
                    }
                    if (wantq)
                    {
                        /* accumulate product of plane rotations in Q */
                        if (initq)
                        {
                            /* take advantage of the fact that Q was */
                            /* initially the Identity matrix */
                            iqend = max(iqend,j2);
                            /* Computing MAX */
                            i__2 = 0;
                            i__3 = k - 3; // , expr subst
                            i2 = max(i__2,i__3);
                            iqaend = i__ * *kd + 1;
                            if (k == 2)
                            {
                                iqaend += *kd;
                            }
                            iqaend = min(iqaend,iqend);
                            i__2 = j2;
                            i__3 = kd1;
                            for (j = j1;
                                    i__3 < 0 ? j >= i__2 : j <= i__2;
                                    j += i__3)
                            {
                                ibl = i__ - i2 / kdm1;
                                ++i2;
                                /* Computing MAX */
                                i__4 = 1;
                                i__5 = j - ibl; // , expr subst
                                iqb = max(i__4,i__5);
                                nq = iqaend + 1 - iqb;
                                /* Computing MIN */
                                i__4 = iqaend + *kd;
                                iqaend = min(i__4,iqend);
                                srot_(&nq, &q[iqb + (j - 1) * q_dim1], &c__1, &q[iqb + j * q_dim1], &c__1, &d__[j], &work[j]);
                                /* L50: */
                            }
                        }
                        else
                        {
                            i__3 = j2;
                            i__2 = kd1;
                            for (j = j1;
                                    i__2 < 0 ? j >= i__3 : j <= i__3;
                                    j += i__2)
                            {
                                srot_(n, &q[(j - 1) * q_dim1 + 1], &c__1, &q[ j * q_dim1 + 1], &c__1, &d__[j], & work[j]);
                                /* L60: */
                            }
                        }
                    }
                    if (j2 + kdn > *n)
                    {
                        /* adjust J2 to keep within the bounds of the matrix */
                        --nr;
                        j2 = j2 - kdn - 1;
                    }
                    i__2 = j2;
                    i__3 = kd1;
                    for (j = j1;
                            i__3 < 0 ? j >= i__2 : j <= i__2;
                            j += i__3)
                    {
                        /* create nonzero element a(j-1,j+kd) outside the band */
                        /* and store it in WORK */
                        work[j + *kd] = work[j] * ab[(j + *kd) * ab_dim1 + 1];
                        ab[(j + *kd) * ab_dim1 + 1] = d__[j] * ab[(j + *kd) * ab_dim1 + 1];
                        /* L70: */
                    }
                    /* L80: */
                }
                /* L90: */
            }
        }
        if (*kd > 0)
        {
            /* copy off-diagonal elements to E */
            i__1 = *n - 1;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                e[i__] = ab[*kd + (i__ + 1) * ab_dim1];
                /* L100: */
            }
        }
        else
        {
            /* set E to zero if original matrix was diagonal */
            i__1 = *n - 1;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                e[i__] = 0.f;
                /* L110: */
            }
        }
        /* copy diagonal elements to D */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            d__[i__] = ab[kd1 + i__ * ab_dim1];
            /* L120: */
        }
    }
    else
    {
        if (*kd > 1)
        {
            /* Reduce to tridiagonal form, working with lower triangle */
            nr = 0;
            j1 = kdn + 2;
            j2 = 1;
            i__1 = *n - 2;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                /* Reduce i-th column of matrix to tridiagonal form */
                for (k = kdn + 1;
                        k >= 2;
                        --k)
                {
                    j1 += kdn;
                    j2 += kdn;
                    if (nr > 0)
                    {
                        /* generate plane rotations to annihilate nonzero */
                        /* elements which have been created outside the band */
                        slargv_(&nr, &ab[kd1 + (j1 - kd1) * ab_dim1], &inca, & work[j1], &kd1, &d__[j1], &kd1);
                        /* apply plane rotations from one side */
                        /* Dependent on the the number of diagonals either */
                        /* SLARTV or SROT is used */
                        if (nr > (*kd << 1) - 1)
                        {
                            i__3 = *kd - 1;
                            for (l = 1;
                                    l <= i__3;
                                    ++l)
                            {
                                slartv_(&nr, &ab[kd1 - l + (j1 - kd1 + l) * ab_dim1], &inca, &ab[kd1 - l + 1 + ( j1 - kd1 + l) * ab_dim1], &inca, &d__[ j1], &work[j1], &kd1);
                                /* L130: */
                            }
                        }
                        else
                        {
                            jend = j1 + kd1 * (nr - 1);
                            i__3 = jend;
                            i__2 = kd1;
                            for (jinc = j1;
                                    i__2 < 0 ? jinc >= i__3 : jinc <= i__3;
                                    jinc += i__2)
                            {
                                srot_(&kdm1, &ab[*kd + (jinc - *kd) * ab_dim1] , &incx, &ab[kd1 + (jinc - *kd) * ab_dim1], &incx, &d__[jinc], &work[ jinc]);
                                /* L140: */
                            }
                        }
                    }
                    if (k > 2)
                    {
                        if (k <= *n - i__ + 1)
                        {
                            /* generate plane rotation to annihilate a(i+k-1,i) */
                            /* within the band */
                            slartg_(&ab[k - 1 + i__ * ab_dim1], &ab[k + i__ * ab_dim1], &d__[i__ + k - 1], &work[i__ + k - 1], &temp);
                            ab[k - 1 + i__ * ab_dim1] = temp;
                            /* apply rotation from the left */
                            i__2 = k - 3;
                            i__3 = *ldab - 1;
                            i__4 = *ldab - 1;
                            srot_(&i__2, &ab[k - 2 + (i__ + 1) * ab_dim1], & i__3, &ab[k - 1 + (i__ + 1) * ab_dim1], & i__4, &d__[i__ + k - 1], &work[i__ + k - 1]);
                        }
                        ++nr;
                        j1 = j1 - kdn - 1;
                    }
                    /* apply plane rotations from both sides to diagonal */
                    /* blocks */
                    if (nr > 0)
                    {
                        slar2v_(&nr, &ab[(j1 - 1) * ab_dim1 + 1], &ab[j1 * ab_dim1 + 1], &ab[(j1 - 1) * ab_dim1 + 2], & inca, &d__[j1], &work[j1], &kd1);
                    }
                    /* apply plane rotations from the right */
                    /* Dependent on the the number of diagonals either */
                    /* SLARTV or SROT is used */
                    if (nr > 0)
                    {
                        if (nr > (*kd << 1) - 1)
                        {
                            i__2 = *kd - 1;
                            for (l = 1;
                                    l <= i__2;
                                    ++l)
                            {
                                if (j2 + l > *n)
                                {
                                    nrt = nr - 1;
                                }
                                else
                                {
                                    nrt = nr;
                                }
                                if (nrt > 0)
                                {
                                    slartv_(&nrt, &ab[l + 2 + (j1 - 1) * ab_dim1], &inca, &ab[l + 1 + j1 * ab_dim1], &inca, &d__[j1], &work[ j1], &kd1);
                                }
                                /* L150: */
                            }
                        }
                        else
                        {
                            j1end = j1 + kd1 * (nr - 2);
                            if (j1end >= j1)
                            {
                                i__2 = j1end;
                                i__3 = kd1;
                                for (j1inc = j1;
                                        i__3 < 0 ? j1inc >= i__2 : j1inc <= i__2;
                                        j1inc += i__3)
                                {
                                    srot_(&kdm1, &ab[(j1inc - 1) * ab_dim1 + 3], &c__1, &ab[j1inc * ab_dim1 + 2], &c__1, &d__[j1inc], &work[ j1inc]);
                                    /* L160: */
                                }
                            }
                            /* Computing MIN */
                            i__3 = kdm1;
                            i__2 = *n - j2; // , expr subst
                            lend = min(i__3,i__2);
                            last = j1end + kd1;
                            if (lend > 0)
                            {
                                srot_(&lend, &ab[(last - 1) * ab_dim1 + 3], & c__1, &ab[last * ab_dim1 + 2], &c__1, &d__[last], &work[last]);
                            }
                        }
                    }
                    if (wantq)
                    {
                        /* accumulate product of plane rotations in Q */
                        if (initq)
                        {
                            /* take advantage of the fact that Q was */
                            /* initially the Identity matrix */
                            iqend = max(iqend,j2);
                            /* Computing MAX */
                            i__3 = 0;
                            i__2 = k - 3; // , expr subst
                            i2 = max(i__3,i__2);
                            iqaend = i__ * *kd + 1;
                            if (k == 2)
                            {
                                iqaend += *kd;
                            }
                            iqaend = min(iqaend,iqend);
                            i__3 = j2;
                            i__2 = kd1;
                            for (j = j1;
                                    i__2 < 0 ? j >= i__3 : j <= i__3;
                                    j += i__2)
                            {
                                ibl = i__ - i2 / kdm1;
                                ++i2;
                                /* Computing MAX */
                                i__4 = 1;
                                i__5 = j - ibl; // , expr subst
                                iqb = max(i__4,i__5);
                                nq = iqaend + 1 - iqb;
                                /* Computing MIN */
                                i__4 = iqaend + *kd;
                                iqaend = min(i__4,iqend);
                                srot_(&nq, &q[iqb + (j - 1) * q_dim1], &c__1, &q[iqb + j * q_dim1], &c__1, &d__[j], &work[j]);
                                /* L170: */
                            }
                        }
                        else
                        {
                            i__2 = j2;
                            i__3 = kd1;
                            for (j = j1;
                                    i__3 < 0 ? j >= i__2 : j <= i__2;
                                    j += i__3)
                            {
                                srot_(n, &q[(j - 1) * q_dim1 + 1], &c__1, &q[ j * q_dim1 + 1], &c__1, &d__[j], & work[j]);
                                /* L180: */
                            }
                        }
                    }
                    if (j2 + kdn > *n)
                    {
                        /* adjust J2 to keep within the bounds of the matrix */
                        --nr;
                        j2 = j2 - kdn - 1;
                    }
                    i__3 = j2;
                    i__2 = kd1;
                    for (j = j1;
                            i__2 < 0 ? j >= i__3 : j <= i__3;
                            j += i__2)
                    {
                        /* create nonzero element a(j+kd,j-1) outside the */
                        /* band and store it in WORK */
                        work[j + *kd] = work[j] * ab[kd1 + j * ab_dim1];
                        ab[kd1 + j * ab_dim1] = d__[j] * ab[kd1 + j * ab_dim1] ;
                        /* L190: */
                    }
                    /* L200: */
                }
                /* L210: */
            }
        }
        if (*kd > 0)
        {
            /* copy off-diagonal elements to E */
            i__1 = *n - 1;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                e[i__] = ab[i__ * ab_dim1 + 2];
                /* L220: */
            }
        }
        else
        {
            /* set E to zero if original matrix was diagonal */
            i__1 = *n - 1;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                e[i__] = 0.f;
                /* L230: */
            }
        }
        /* copy diagonal elements to D */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            d__[i__] = ab[i__ * ab_dim1 + 1];
            /* L240: */
        }
    }
    return 0;
    /* End of SSBTRD */
}
/* ssbtrd_ */
