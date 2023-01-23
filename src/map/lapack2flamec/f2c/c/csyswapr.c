/* csyswapr.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CSYSWAPR */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSYSWAPR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyswap r.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyswap r.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyswap r.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSYSWAPR( UPLO, N, A, LDA, I1, I2) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER I1, I2, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, N ) */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYSWAPR applies an elementary permutation on the rows and the columns of */
/* > a symmetric matrix. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**T;
*/
/* > = 'L': Lower triangular, form is A = L*D*L**T. */
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
/* > On entry, the N-by-N matrix A. On exit, the permuted matrix */
/* > where the rows I1 and I2 and columns I1 and I2 are interchanged. */
/* > If UPLO = 'U', the interchanges are applied to the upper */
/* > triangular part and the strictly lower triangular part of A is */
/* > not referenced;
if UPLO = 'L', the interchanges are applied to */
/* > the lower triangular part and the part of A above the diagonal */
/* > is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] I1 */
/* > \verbatim */
/* > I1 is INTEGER */
/* > Index of the first row to swap */
/* > \endverbatim */
/* > */
/* > \param[in] I2 */
/* > \verbatim */
/* > I2 is INTEGER */
/* > Index of the second row to swap */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup complexSYauxiliary */
/* ===================================================================== */
/* Subroutine */
int csyswapr_(char *uplo, integer *n, complex *a, integer * lda, integer *i1, integer *i2)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("csyswapr inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS ", i1 %" FLA_IS ", i2 %" FLA_IS "",*uplo, *n, *lda, *i1, *i2);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    complex tmp;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int cswap_(integer *, complex *, integer *, complex *, integer *);
    logical upper;
    /* -- LAPACK auxiliary routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* ===================================================================== */
    /* .. */
    /* .. Local Scalars .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    /* Function Body */
    upper = lsame_(uplo, "U");
    if (upper)
    {
        /* UPPER */
        /* first swap */
        /* - swap column I1 and I2 from I1 to I1-1 */
        i__1 = *i1 - 1;
        cswap_(&i__1, &a[*i1 * a_dim1 + 1], &c__1, &a[*i2 * a_dim1 + 1], & c__1);
        /* second swap : */
        /* - swap A(I1,I1) and A(I2,I2) */
        /* - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1 */
        i__1 = *i1 + *i1 * a_dim1;
        tmp.r = a[i__1].r;
        tmp.i = a[i__1].i; // , expr subst
        i__1 = *i1 + *i1 * a_dim1;
        i__2 = *i2 + *i2 * a_dim1;
        a[i__1].r = a[i__2].r;
        a[i__1].i = a[i__2].i; // , expr subst
        i__1 = *i2 + *i2 * a_dim1;
        a[i__1].r = tmp.r;
        a[i__1].i = tmp.i; // , expr subst
        i__1 = *i2 - *i1 - 1;
        cswap_(&i__1, &a[*i1 + (*i1 + 1) * a_dim1], lda, &a[*i1 + 1 + *i2 * a_dim1], &c__1);
        /* third swap */
        /* - swap row I1 and I2 from I2+1 to N */
        if (*i2 < *n)
        {
            i__1 = *n - *i2;
            cswap_(&i__1, &a[*i1 + (*i2 + 1) * a_dim1], lda, &a[*i2 + (*i2 + 1) * a_dim1], lda);
        }
    }
    else
    {
        /* LOWER */
        /* first swap */
        /* - swap row I1 and I2 from I1 to I1-1 */
        i__1 = *i1 - 1;
        cswap_(&i__1, &a[*i1 + a_dim1], lda, &a[*i2 + a_dim1], lda);
        /* second swap : */
        /* - swap A(I1,I1) and A(I2,I2) */
        /* - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1 */
        i__1 = *i1 + *i1 * a_dim1;
        tmp.r = a[i__1].r;
        tmp.i = a[i__1].i; // , expr subst
        i__1 = *i1 + *i1 * a_dim1;
        i__2 = *i2 + *i2 * a_dim1;
        a[i__1].r = a[i__2].r;
        a[i__1].i = a[i__2].i; // , expr subst
        i__1 = *i2 + *i2 * a_dim1;
        a[i__1].r = tmp.r;
        a[i__1].i = tmp.i; // , expr subst
        i__1 = *i2 - *i1 - 1;
        cswap_(&i__1, &a[*i1 + 1 + *i1 * a_dim1], &c__1, &a[*i2 + (*i1 + 1) * a_dim1], lda);
        /* third swap */
        /* - swap col I1 and I2 from I2+1 to N */
        if (*i2 < *n)
        {
            i__1 = *n - *i2;
            cswap_(&i__1, &a[*i2 + 1 + *i1 * a_dim1], &c__1, &a[*i2 + 1 + *i2 * a_dim1], &c__1);
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
}
/* csyswapr_ */
