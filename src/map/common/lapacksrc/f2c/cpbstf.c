/* ../netlib/cpbstf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b9 = -1.f;
/* > \brief \b CPBSTF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CPBSTF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbstf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbstf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbstf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CPBSTF( UPLO, N, KD, AB, LDAB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, KD, LDAB, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX AB( LDAB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPBSTF computes a split Cholesky factorization of a complex */
/* > Hermitian positive definite band matrix A. */
/* > */
/* > This routine is designed to be used in conjunction with CHBGST. */
/* > */
/* > The factorization has the form A = S**H*S where S is a band matrix */
/* > of the same bandwidth as A and the following structure: */
/* > */
/* > S = ( U ) */
/* > ( M L ) */
/* > */
/* > where U is upper triangular of order m = (n+kd)/2, and L is lower */
/* > triangular of order n-m. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
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
/* > AB is COMPLEX array, dimension (LDAB,N) */
/* > On entry, the upper or lower triangle of the Hermitian band */
/* > matrix A, stored in the first kd+1 rows of the array. The */
/* > j-th column of A is stored in the j-th column of the array AB */
/* > as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*/
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd). */
/* > */
/* > On exit, if INFO = 0, the factor S from the split Cholesky */
/* > factorization A = S**H*S. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the factorization could not be completed, */
/* > because the updated element a(i,i) was negative;
the */
/* > matrix A is not positive definite. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The band storage scheme is illustrated by the following example, when */
/* > N = 7, KD = 2: */
/* > */
/* > S = ( s11 s12 s13 ) */
/* > ( s22 s23 s24 ) */
/* > ( s33 s34 ) */
/* > ( s44 ) */
/* > ( s53 s54 s55 ) */
/* > ( s64 s65 s66 ) */
/* > ( s75 s76 s77 ) */
/* > */
/* > If UPLO = 'U', the array AB holds: */
/* > */
/* > on entry: on exit: */
/* > */
/* > * * a13 a24 a35 a46 a57 * * s13 s24 s53**H s64**H s75**H */
/* > * a12 a23 a34 a45 a56 a67 * s12 s23 s34 s54**H s65**H s76**H */
/* > a11 a22 a33 a44 a55 a66 a77 s11 s22 s33 s44 s55 s66 s77 */
/* > */
/* > If UPLO = 'L', the array AB holds: */
/* > */
/* > on entry: on exit: */
/* > */
/* > a11 a22 a33 a44 a55 a66 a77 s11 s22 s33 s44 s55 s66 s77 */
/* > a21 a32 a43 a54 a65 a76 * s12**H s23**H s34**H s54 s65 s76 * */
/* > a31 a42 a53 a64 a64 * * s13**H s24**H s53 s64 s75 * * */
/* > */
/* > Array elements marked * are not used by the routine;
s12**H denotes */
/* > conjg(s12);
the diagonal elements of S are real. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int cpbstf_(char *uplo, integer *n, integer *kd, complex *ab, integer *ldab, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3;
    real r__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer j, m, km;
    real ajj;
    integer kld;
    extern /* Subroutine */
    int cher_(char *, integer *, real *, complex *, integer *, complex *, integer *);
    extern logical lsame_(char *, char *);
    logical upper;
    extern /* Subroutine */
    int clacgv_(integer *, complex *, integer *), csscal_(integer *, real *, complex *, integer *), xerbla_(char *, integer *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*kd < 0)
    {
        *info = -3;
    }
    else if (*ldab < *kd + 1)
    {
        *info = -5;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CPBSTF", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Computing MAX */
    i__1 = 1;
    i__2 = *ldab - 1; // , expr subst
    kld = max(i__1,i__2);
    /* Set the splitting point m. */
    m = (*n + *kd) / 2;
    if (upper)
    {
        /* Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m). */
        i__1 = m + 1;
        for (j = *n;
                j >= i__1;
                --j)
        {
            /* Compute s(j,j) and test for non-positive-definiteness. */
            i__2 = *kd + 1 + j * ab_dim1;
            ajj = ab[i__2].r;
            if (ajj <= 0.f)
            {
                i__2 = *kd + 1 + j * ab_dim1;
                ab[i__2].r = ajj;
                ab[i__2].i = 0.f; // , expr subst
                goto L50;
            }
            ajj = sqrt(ajj);
            i__2 = *kd + 1 + j * ab_dim1;
            ab[i__2].r = ajj;
            ab[i__2].i = 0.f; // , expr subst
            /* Computing MIN */
            i__2 = j - 1;
            km = min(i__2,*kd);
            /* Compute elements j-km:j-1 of the j-th column and update the */
            /* the leading submatrix within the band. */
            r__1 = 1.f / ajj;
            csscal_(&km, &r__1, &ab[*kd + 1 - km + j * ab_dim1], &c__1);
            cher_("Upper", &km, &c_b9, &ab[*kd + 1 - km + j * ab_dim1], &c__1, &ab[*kd + 1 + (j - km) * ab_dim1], &kld);
            /* L10: */
        }
        /* Factorize the updated submatrix A(1:m,1:m) as U**H*U. */
        i__1 = m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Compute s(j,j) and test for non-positive-definiteness. */
            i__2 = *kd + 1 + j * ab_dim1;
            ajj = ab[i__2].r;
            if (ajj <= 0.f)
            {
                i__2 = *kd + 1 + j * ab_dim1;
                ab[i__2].r = ajj;
                ab[i__2].i = 0.f; // , expr subst
                goto L50;
            }
            ajj = sqrt(ajj);
            i__2 = *kd + 1 + j * ab_dim1;
            ab[i__2].r = ajj;
            ab[i__2].i = 0.f; // , expr subst
            /* Computing MIN */
            i__2 = *kd;
            i__3 = m - j; // , expr subst
            km = min(i__2,i__3);
            /* Compute elements j+1:j+km of the j-th row and update the */
            /* trailing submatrix within the band. */
            if (km > 0)
            {
                r__1 = 1.f / ajj;
                csscal_(&km, &r__1, &ab[*kd + (j + 1) * ab_dim1], &kld);
                clacgv_(&km, &ab[*kd + (j + 1) * ab_dim1], &kld);
                cher_("Upper", &km, &c_b9, &ab[*kd + (j + 1) * ab_dim1], &kld, &ab[*kd + 1 + (j + 1) * ab_dim1], &kld);
                clacgv_(&km, &ab[*kd + (j + 1) * ab_dim1], &kld);
            }
            /* L20: */
        }
    }
    else
    {
        /* Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m). */
        i__1 = m + 1;
        for (j = *n;
                j >= i__1;
                --j)
        {
            /* Compute s(j,j) and test for non-positive-definiteness. */
            i__2 = j * ab_dim1 + 1;
            ajj = ab[i__2].r;
            if (ajj <= 0.f)
            {
                i__2 = j * ab_dim1 + 1;
                ab[i__2].r = ajj;
                ab[i__2].i = 0.f; // , expr subst
                goto L50;
            }
            ajj = sqrt(ajj);
            i__2 = j * ab_dim1 + 1;
            ab[i__2].r = ajj;
            ab[i__2].i = 0.f; // , expr subst
            /* Computing MIN */
            i__2 = j - 1;
            km = min(i__2,*kd);
            /* Compute elements j-km:j-1 of the j-th row and update the */
            /* trailing submatrix within the band. */
            r__1 = 1.f / ajj;
            csscal_(&km, &r__1, &ab[km + 1 + (j - km) * ab_dim1], &kld);
            clacgv_(&km, &ab[km + 1 + (j - km) * ab_dim1], &kld);
            cher_("Lower", &km, &c_b9, &ab[km + 1 + (j - km) * ab_dim1], &kld, &ab[(j - km) * ab_dim1 + 1], &kld);
            clacgv_(&km, &ab[km + 1 + (j - km) * ab_dim1], &kld);
            /* L30: */
        }
        /* Factorize the updated submatrix A(1:m,1:m) as U**H*U. */
        i__1 = m;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Compute s(j,j) and test for non-positive-definiteness. */
            i__2 = j * ab_dim1 + 1;
            ajj = ab[i__2].r;
            if (ajj <= 0.f)
            {
                i__2 = j * ab_dim1 + 1;
                ab[i__2].r = ajj;
                ab[i__2].i = 0.f; // , expr subst
                goto L50;
            }
            ajj = sqrt(ajj);
            i__2 = j * ab_dim1 + 1;
            ab[i__2].r = ajj;
            ab[i__2].i = 0.f; // , expr subst
            /* Computing MIN */
            i__2 = *kd;
            i__3 = m - j; // , expr subst
            km = min(i__2,i__3);
            /* Compute elements j+1:j+km of the j-th column and update the */
            /* trailing submatrix within the band. */
            if (km > 0)
            {
                r__1 = 1.f / ajj;
                csscal_(&km, &r__1, &ab[j * ab_dim1 + 2], &c__1);
                cher_("Lower", &km, &c_b9, &ab[j * ab_dim1 + 2], &c__1, &ab[( j + 1) * ab_dim1 + 1], &kld);
            }
            /* L40: */
        }
    }
    return 0;
L50:
    *info = j;
    return 0;
    /* End of CPBSTF */
}
/* cpbstf_ */
