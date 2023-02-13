/* ../netlib/zlat2c.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLAT2C converts a double complex triangular matrix to a complex triangular matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLAT2C + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlat2c. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlat2c. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlat2c. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLAT2C( UPLO, N, A, LDA, SA, LDSA, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDSA, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX SA( LDSA, * ) */
/* COMPLEX*16 A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAT2C converts a COMPLEX*16 triangular matrix, SA, to a COMPLEX */
/* > triangular matrix, A. */
/* > */
/* > RMAX is the overflow for the SINGLE PRECISION arithmetic */
/* > ZLAT2C checks that all the entries of A are between -RMAX and */
/* > RMAX. If not the convertion is aborted and a flag is raised. */
/* > */
/* > This is an auxiliary routine so there is no argument checking. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': A is upper triangular;
*/
/* > = 'L': A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of rows and columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the N-by-N triangular coefficient matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] SA */
/* > \verbatim */
/* > SA is COMPLEX array, dimension (LDSA,N) */
/* > Only the UPLO part of SA is referenced. On exit, if INFO=0, */
/* > the N-by-N coefficient matrix SA;
if INFO>0, the content of */
/* > the UPLO part of SA is unspecified. */
/* > \endverbatim */
/* > */
/* > \param[in] LDSA */
/* > \verbatim */
/* > LDSA is INTEGER */
/* > The leading dimension of the array SA. LDSA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > = 1: an entry of the matrix A is greater than the SINGLE */
/* > PRECISION overflow threshold, in this case, the content */
/* > of the UPLO part of SA in exit is unspecified. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int zlat2c_(char *uplo, integer *n, doublecomplex *a, integer *lda, complex *sa, integer *ldsa, integer *info)
{
    /* System generated locals */
    integer sa_dim1, sa_offset, a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    integer i__, j;
    doublereal rmax;
    extern logical lsame_(char *, char *);
    logical upper;
    extern real slamch_(char *);
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    sa_dim1 = *ldsa;
    sa_offset = 1 + sa_dim1;
    sa -= sa_offset;
    /* Function Body */
    rmax = slamch_("O");
    upper = lsame_(uplo, "U");
    if (upper)
    {
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * a_dim1;
                i__4 = i__ + j * a_dim1;
                if (a[i__3].r < -rmax || a[i__4].r > rmax || d_imag(&a[i__ + j * a_dim1]) < -rmax || d_imag(&a[i__ + j * a_dim1]) > rmax)
                {
                    *info = 1;
                    goto L50;
                }
                i__3 = i__ + j * sa_dim1;
                i__4 = i__ + j * a_dim1;
                sa[i__3].r = a[i__4].r;
                sa[i__3].i = a[i__4].i; // , expr subst
                /* L10: */
            }
            /* L20: */
        }
    }
    else
    {
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *n;
            for (i__ = j;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * a_dim1;
                i__4 = i__ + j * a_dim1;
                if (a[i__3].r < -rmax || a[i__4].r > rmax || d_imag(&a[i__ + j * a_dim1]) < -rmax || d_imag(&a[i__ + j * a_dim1]) > rmax)
                {
                    *info = 1;
                    goto L50;
                }
                i__3 = i__ + j * sa_dim1;
                i__4 = i__ + j * a_dim1;
                sa[i__3].r = a[i__4].r;
                sa[i__3].i = a[i__4].i; // , expr subst
                /* L30: */
            }
            /* L40: */
        }
    }
L50:
    return 0;
    /* End of ZLAT2C */
}
/* zlat2c_ */
