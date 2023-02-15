/* ../netlib/claset.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given val ues. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLASET + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claset. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claset. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claset. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLASET( UPLO, M, N, ALPHA, BETA, A, LDA ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER LDA, M, N */
/* COMPLEX ALPHA, BETA */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLASET initializes a 2-D array A to BETA on the diagonal and */
/* > ALPHA on the offdiagonals. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies the part of the matrix A to be set. */
/* > = 'U': Upper triangular part is set. The lower triangle */
/* > is unchanged. */
/* > = 'L': Lower triangular part is set. The upper triangle */
/* > is unchanged. */
/* > Otherwise: All of the matrix A is set. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > On entry, M specifies the number of rows of A. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > On entry, N specifies the number of columns of A. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is COMPLEX */
/* > All the offdiagonal array elements are set to ALPHA. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* > BETA is COMPLEX */
/* > All the diagonal array elements are set to BETA. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the m by n matrix A. */
/* > On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
*/
/* > A(i,i) = BETA , 1 <= i <= min(m,n) */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int claset_(char *uplo, integer *m, integer *n, complex * alpha, complex *beta, complex *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    integer i__, j;
    extern logical lsame_(char *, char *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    /* Function Body */
    if (lsame_(uplo, "U"))
    {
        /* Set the diagonal to BETA and the strictly upper triangular */
        /* part of the array to ALPHA. */
        i__1 = *n;
        for (j = 2;
                j <= i__1;
                ++j)
        {
            /* Computing MIN */
            i__3 = j - 1;
            i__2 = min(i__3,*m);
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * a_dim1;
                a[i__3].r = alpha->r;
                a[i__3].i = alpha->i; // , expr subst
                /* L10: */
            }
            /* L20: */
        }
        i__1 = min(*n,*m);
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__ + i__ * a_dim1;
            a[i__2].r = beta->r;
            a[i__2].i = beta->i; // , expr subst
            /* L30: */
        }
    }
    else if (lsame_(uplo, "L"))
    {
        /* Set the diagonal to BETA and the strictly lower triangular */
        /* part of the array to ALPHA. */
        i__1 = min(*m,*n);
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = j + 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * a_dim1;
                a[i__3].r = alpha->r;
                a[i__3].i = alpha->i; // , expr subst
                /* L40: */
            }
            /* L50: */
        }
        i__1 = min(*n,*m);
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__ + i__ * a_dim1;
            a[i__2].r = beta->r;
            a[i__2].i = beta->i; // , expr subst
            /* L60: */
        }
    }
    else
    {
        /* Set the array to BETA on the diagonal and ALPHA on the */
        /* offdiagonal. */
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * a_dim1;
                a[i__3].r = alpha->r;
                a[i__3].i = alpha->i; // , expr subst
                /* L70: */
            }
            /* L80: */
        }
        i__1 = min(*m,*n);
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = i__ + i__ * a_dim1;
            a[i__2].r = beta->r;
            a[i__2].i = beta->i; // , expr subst
            /* L90: */
        }
    }
    return 0;
    /* End of CLASET */
}
/* claset_ */
