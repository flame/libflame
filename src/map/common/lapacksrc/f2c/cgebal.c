/* ../netlib/cgebal.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CGEBAL */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGEBAL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgebal. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgebal. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgebal. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOB */
/* INTEGER IHI, ILO, INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* REAL SCALE( * ) */
/* COMPLEX A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEBAL balances a general complex matrix A. This involves, first, */
/* > permuting A by a similarity transformation to isolate eigenvalues */
/* > in the first 1 to ILO-1 and last IHI+1 to N elements on the */
/* > diagonal;
and second, applying a diagonal similarity transformation */
/* > to rows and columns ILO to IHI to make the rows and columns as */
/* > close in norm as possible. Both steps are optional. */
/* > */
/* > Balancing may reduce the 1-norm of the matrix, and improve the */
/* > accuracy of the computed eigenvalues and/or eigenvectors. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOB */
/* > \verbatim */
/* > JOB is CHARACTER*1 */
/* > Specifies the operations to be performed on A: */
/* > = 'N': none: simply set ILO = 1, IHI = N, SCALE(I) = 1.0 */
/* > for i = 1,...,N;
*/
/* > = 'P': permute only;
*/
/* > = 'S': scale only;
*/
/* > = 'B': both permute and scale. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the input matrix A. */
/* > On exit, A is overwritten by the balanced matrix. */
/* > If JOB = 'N', A is not referenced. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > \param[out] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > ILO and IHI are set to integers such that on exit */
/* > A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N. */
/* > If JOB = 'N' or 'S', ILO = 1 and IHI = N. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE */
/* > \verbatim */
/* > SCALE is REAL array, dimension (N) */
/* > Details of the permutations and scaling factors applied to */
/* > A. If P(j) is the index of the row and column interchanged */
/* > with row and column j and D(j) is the scaling factor */
/* > applied to row and column j, then */
/* > SCALE(j) = P(j) for j = 1,...,ILO-1 */
/* > = D(j) for j = ILO,...,IHI */
/* > = P(j) for j = IHI+1,...,N. */
/* > The order in which the interchanges are made is N to IHI+1, */
/* > then 1 to ILO-1. */
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
/* > \date November 2013 */
/* > \ingroup complexGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The permutations consist of row and column interchanges which put */
/* > the matrix in the form */
/* > */
/* > ( T1 X Y ) */
/* > P A P = ( 0 B Z ) */
/* > ( 0 0 T2 ) */
/* > */
/* > where T1 and T2 are upper triangular matrices whose eigenvalues lie */
/* > along the diagonal. The column indices ILO and IHI mark the starting */
/* > and ending columns of the submatrix B. Balancing consists of applying */
/* > a diagonal similarity transformation inv(D) * B * D to make the */
/* > 1-norms of each row of B and its corresponding column nearly equal. */
/* > The output matrix is */
/* > */
/* > ( T1 X*D Y ) */
/* > ( 0 inv(D)*B*D inv(D)*Z ). */
/* > ( 0 0 T2 ) */
/* > */
/* > Information about the permutations P and the diagonal matrix D is */
/* > returned in the vector SCALE. */
/* > */
/* > This subroutine is based on the EISPACK routine CBAL. */
/* > */
/* > Modified by Tzu-Yi Chen, Computer Science Division, University of */
/* > California at Berkeley, USA */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int cgebal_(char *job, integer *n, complex *a, integer *lda, integer *ilo, integer *ihi, real *scale, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2;
    /* Builtin functions */
    double r_imag(complex *), c_abs(complex *);
    /* Local variables */
    real c__, f, g;
    integer i__, j, k, l, m;
    real r__, s, ca, ra;
    integer ica, ira, iexc;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int cswap_(integer *, complex *, integer *, complex *, integer *);
    real sfmin1, sfmin2, sfmax1, sfmax2;
    extern real scnrm2_(integer *, complex *, integer *);
    extern integer icamax_(integer *, complex *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int csscal_(integer *, real *, complex *, integer *), xerbla_(char *, integer *);
    extern logical sisnan_(real *);
    logical noconv;
    /* -- LAPACK computational routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --scale;
    /* Function Body */
    *info = 0;
    if (! lsame_(job, "N") && ! lsame_(job, "P") && ! lsame_(job, "S") && ! lsame_(job, "B"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < max(1,*n))
    {
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGEBAL", &i__1);
        return 0;
    }
    k = 1;
    l = *n;
    if (*n == 0)
    {
        goto L210;
    }
    if (lsame_(job, "N"))
    {
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            scale[i__] = 1.f;
            /* L10: */
        }
        goto L210;
    }
    if (lsame_(job, "S"))
    {
        goto L120;
    }
    /* Permutation to isolate eigenvalues if possible */
    goto L50;
    /* Row and column exchange. */
L20:
    scale[m] = (real) j;
    if (j == m)
    {
        goto L30;
    }
    cswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
    i__1 = *n - k + 1;
    cswap_(&i__1, &a[j + k * a_dim1], lda, &a[m + k * a_dim1], lda);
L30:
    switch (iexc)
    {
    case 1:
        goto L40;
    case 2:
        goto L80;
    }
    /* Search for rows isolating an eigenvalue and push them down. */
L40:
    if (l == 1)
    {
        goto L210;
    }
    --l;
L50:
    for (j = l;
            j >= 1;
            --j)
    {
        i__1 = l;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            if (i__ == j)
            {
                goto L60;
            }
            i__2 = j + i__ * a_dim1;
            if (a[i__2].r != 0.f || r_imag(&a[j + i__ * a_dim1]) != 0.f)
            {
                goto L70;
            }
L60:
            ;
        }
        m = l;
        iexc = 1;
        goto L20;
L70:
        ;
    }
    goto L90;
    /* Search for columns isolating an eigenvalue and push them left. */
L80:
    ++k;
L90:
    i__1 = l;
    for (j = k;
            j <= i__1;
            ++j)
    {
        i__2 = l;
        for (i__ = k;
                i__ <= i__2;
                ++i__)
        {
            if (i__ == j)
            {
                goto L100;
            }
            i__3 = i__ + j * a_dim1;
            if (a[i__3].r != 0.f || r_imag(&a[i__ + j * a_dim1]) != 0.f)
            {
                goto L110;
            }
L100:
            ;
        }
        m = k;
        iexc = 2;
        goto L20;
L110:
        ;
    }
L120:
    i__1 = l;
    for (i__ = k;
            i__ <= i__1;
            ++i__)
    {
        scale[i__] = 1.f;
        /* L130: */
    }
    if (lsame_(job, "P"))
    {
        goto L210;
    }
    /* Balance the submatrix in rows K to L. */
    /* Iterative loop for norm reduction */
    sfmin1 = slamch_("S") / slamch_("P");
    sfmax1 = 1.f / sfmin1;
    sfmin2 = sfmin1 * 2.f;
    sfmax2 = 1.f / sfmin2;
L140:
    noconv = FALSE_;
    i__1 = l;
    for (i__ = k;
            i__ <= i__1;
            ++i__)
    {
        i__2 = l - k + 1;
        c__ = scnrm2_(&i__2, &a[k + i__ * a_dim1], &c__1);
        i__2 = l - k + 1;
        r__ = scnrm2_(&i__2, &a[i__ + k * a_dim1], lda);
        ica = icamax_(&l, &a[i__ * a_dim1 + 1], &c__1);
        ca = c_abs(&a[ica + i__ * a_dim1]);
        i__2 = *n - k + 1;
        ira = icamax_(&i__2, &a[i__ + k * a_dim1], lda);
        ra = c_abs(&a[i__ + (ira + k - 1) * a_dim1]);
        /* Guard against zero C or R due to underflow. */
        if (c__ == 0.f || r__ == 0.f)
        {
            goto L200;
        }
        g = r__ / 2.f;
        f = 1.f;
        s = c__ + r__;
L160: /* Computing MAX */
        r__1 = max(f,c__);
        /* Computing MIN */
        r__2 = min(r__,g);
        if (c__ >= g || max(r__1,ca) >= sfmax2 || min(r__2,ra) <= sfmin2)
        {
            goto L170;
        }
        r__1 = c__ + f + ca + r__ + g + ra;
        if (sisnan_(&r__1))
        {
            /* Exit if NaN to avoid infinite loop */
            *info = -3;
            i__2 = -(*info);
            xerbla_("CGEBAL", &i__2);
            return 0;
        }
        f *= 2.f;
        c__ *= 2.f;
        ca *= 2.f;
        r__ /= 2.f;
        g /= 2.f;
        ra /= 2.f;
        goto L160;
L170:
        g = c__ / 2.f;
L180: /* Computing MIN */
        r__1 = min(f,c__);
        r__1 = min(r__1,g); // , expr subst
        if (g < r__ || max(r__,ra) >= sfmax2 || min(r__1,ca) <= sfmin2)
        {
            goto L190;
        }
        f /= 2.f;
        c__ /= 2.f;
        g /= 2.f;
        ca /= 2.f;
        r__ *= 2.f;
        ra *= 2.f;
        goto L180;
        /* Now balance. */
L190:
        if (c__ + r__ >= s * .95f)
        {
            goto L200;
        }
        if (f < 1.f && scale[i__] < 1.f)
        {
            if (f * scale[i__] <= sfmin1)
            {
                goto L200;
            }
        }
        if (f > 1.f && scale[i__] > 1.f)
        {
            if (scale[i__] >= sfmax1 / f)
            {
                goto L200;
            }
        }
        g = 1.f / f;
        scale[i__] *= f;
        noconv = TRUE_;
        i__2 = *n - k + 1;
        csscal_(&i__2, &g, &a[i__ + k * a_dim1], lda);
        csscal_(&l, &f, &a[i__ * a_dim1 + 1], &c__1);
L200:
        ;
    }
    if (noconv)
    {
        goto L140;
    }
L210:
    *ilo = k;
    *ihi = l;
    return 0;
    /* End of CGEBAL */
}
/* cgebal_ */
