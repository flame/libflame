/* ../netlib/dsyswapr.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DSYSWAPR applies an elementary permutation on the rows and columns of a symmetric matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSYSWAPR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyswap r.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyswap r.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyswap r.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSYSWAPR( UPLO, N, A, LDA, I1, I2) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER I1, I2, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, N ) */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYSWAPR applies an elementary permutation on the rows and the columns of */
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
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the NB diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by DSYTRF. */
/* > */
/* > On exit, if INFO = 0, the (symmetric) inverse of the original */
/* > matrix. If UPLO = 'U', the upper triangular part of the */
/* > inverse is formed and the part of A below the diagonal is not */
/* > referenced;
if UPLO = 'L' the lower triangular part of the */
/* > inverse is formed and the part of A above the diagonal is */
/* > not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
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
/* > \date September 2012 */
/* > \ingroup doubleSYauxiliary */
/* ===================================================================== */
/* Subroutine */
int dsyswapr_(char *uplo, integer *n, doublereal *a, integer *lda, integer *i1, integer *i2)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    integer i__;
    doublereal tmp;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical upper;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
        dswap_(&i__1, &a[*i1 * a_dim1 + 1], &c__1, &a[*i2 * a_dim1 + 1], & c__1);
        /* second swap : */
        /* - swap A(I1,I1) and A(I2,I2) */
        /* - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1 */
        tmp = a[*i1 + *i1 * a_dim1];
        a[*i1 + *i1 * a_dim1] = a[*i2 + *i2 * a_dim1];
        a[*i2 + *i2 * a_dim1] = tmp;
        i__1 = *i2 - *i1 - 1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            tmp = a[*i1 + (*i1 + i__) * a_dim1];
            a[*i1 + (*i1 + i__) * a_dim1] = a[*i1 + i__ + *i2 * a_dim1];
            a[*i1 + i__ + *i2 * a_dim1] = tmp;
        }
        /* third swap */
        /* - swap row I1 and I2 from I2+1 to N */
        i__1 = *n;
        for (i__ = *i2 + 1;
                i__ <= i__1;
                ++i__)
        {
            tmp = a[*i1 + i__ * a_dim1];
            a[*i1 + i__ * a_dim1] = a[*i2 + i__ * a_dim1];
            a[*i2 + i__ * a_dim1] = tmp;
        }
    }
    else
    {
        /* LOWER */
        /* first swap */
        /* - swap row I1 and I2 from I1 to I1-1 */
        i__1 = *i1 - 1;
        dswap_(&i__1, &a[*i1 + a_dim1], lda, &a[*i2 + a_dim1], lda);
        /* second swap : */
        /* - swap A(I1,I1) and A(I2,I2) */
        /* - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1 */
        tmp = a[*i1 + *i1 * a_dim1];
        a[*i1 + *i1 * a_dim1] = a[*i2 + *i2 * a_dim1];
        a[*i2 + *i2 * a_dim1] = tmp;
        i__1 = *i2 - *i1 - 1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            tmp = a[*i1 + i__ + *i1 * a_dim1];
            a[*i1 + i__ + *i1 * a_dim1] = a[*i2 + (*i1 + i__) * a_dim1];
            a[*i2 + (*i1 + i__) * a_dim1] = tmp;
        }
        /* third swap */
        /* - swap col I1 and I2 from I2+1 to N */
        i__1 = *n;
        for (i__ = *i2 + 1;
                i__ <= i__1;
                ++i__)
        {
            tmp = a[i__ + *i1 * a_dim1];
            a[i__ + *i1 * a_dim1] = a[i__ + *i2 * a_dim1];
            a[i__ + *i2 * a_dim1] = tmp;
        }
    }
    return 0;
}
/* dsyswapr_ */
