/* ../netlib/zlacp2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLACP2 copies all or part of a real two-dimensional array to a complex array. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLACP2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacp2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacp2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacp2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLACP2( UPLO, M, N, A, LDA, B, LDB ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER LDA, LDB, M, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ) */
/* COMPLEX*16 B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLACP2 copies all or part of a real two-dimensional matrix A to a */
/* > complex matrix B. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies the part of the matrix A to be copied to B. */
/* > = 'U': Upper triangular part */
/* > = 'L': Lower triangular part */
/* > Otherwise: All of the matrix A */
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
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > The m by n matrix A. If UPLO = 'U', only the upper trapezium */
/* > is accessed;
if UPLO = 'L', only the lower trapezium is */
/* > accessed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,N) */
/* > On exit, B = A in the locations specified by UPLO. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,M). */
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
int zlacp2_(char *uplo, integer *m, integer *n, doublereal * a, integer *lda, doublecomplex *b, integer *ldb)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlacp2 inputs: uplo %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "",*uplo, *m, *n, *lda, *ldb);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    if (lsame_(uplo, "U"))
    {
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = fla_min(j,*m);
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * a_dim1;
                b[i__3].r = a[i__4];
                b[i__3].i = 0.; // , expr subst
                /* L10: */
            }
            /* L20: */
        }
    }
    else if (lsame_(uplo, "L"))
    {
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = j;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * a_dim1;
                b[i__3].r = a[i__4];
                b[i__3].i = 0.; // , expr subst
                /* L30: */
            }
            /* L40: */
        }
    }
    else
    {
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
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * a_dim1;
                b[i__3].r = a[i__4];
                b[i__3].i = 0.; // , expr subst
                /* L50: */
            }
            /* L60: */
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of ZLACP2 */
}
/* zlacp2_ */
