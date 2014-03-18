/* ../netlib/dlaqge.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLAQGE scales a general rectangular matrix, using row and column scaling factors computed by sg eequ. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAQGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqge. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqge. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqge. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/* EQUED ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED */
/* INTEGER LDA, M, N */
/* DOUBLE PRECISION AMAX, COLCND, ROWCND */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), C( * ), R( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQGE equilibrates a general M by N matrix A using the row and */
/* > column scaling factors in the vectors R and C. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the M by N matrix A. */
/* > On exit, the equilibrated matrix. See EQUED for the form of */
/* > the equilibrated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(M,1). */
/* > \endverbatim */
/* > */
/* > \param[in] R */
/* > \verbatim */
/* > R is DOUBLE PRECISION array, dimension (M) */
/* > The row scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (N) */
/* > The column scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] ROWCND */
/* > \verbatim */
/* > ROWCND is DOUBLE PRECISION */
/* > Ratio of the smallest R(i) to the largest R(i). */
/* > \endverbatim */
/* > */
/* > \param[in] COLCND */
/* > \verbatim */
/* > COLCND is DOUBLE PRECISION */
/* > Ratio of the smallest C(i) to the largest C(i). */
/* > \endverbatim */
/* > */
/* > \param[in] AMAX */
/* > \verbatim */
/* > AMAX is DOUBLE PRECISION */
/* > Absolute value of largest matrix entry. */
/* > \endverbatim */
/* > */
/* > \param[out] EQUED */
/* > \verbatim */
/* > EQUED is CHARACTER*1 */
/* > Specifies the form of equilibration that was done. */
/* > = 'N': No equilibration */
/* > = 'R': Row equilibration, i.e., A has been premultiplied by */
/* > diag(R). */
/* > = 'C': Column equilibration, i.e., A has been postmultiplied */
/* > by diag(C). */
/* > = 'B': Both row and column equilibration, i.e., A has been */
/* > replaced by diag(R) * A * diag(C). */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > THRESH is a threshold value used to decide if row or column scaling */
/* > should be done based on the ratio of the row or column scaling */
/* > factors. If ROWCND < THRESH, row scaling is done, and if */
/* > COLCND < THRESH, column scaling is done. */
/* > */
/* > LARGE and SMALL are threshold values used to decide if row scaling */
/* > should be done based on the absolute size of the largest matrix */
/* > element. If AMAX > LARGE or AMAX < SMALL, row scaling is done. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleGEauxiliary */
/* ===================================================================== */
/* Subroutine */
int dlaqge_(integer *m, integer *n, doublereal *a, integer * lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    integer i__, j;
    doublereal cj, large, small;
    extern doublereal dlamch_(char *);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --r__;
    --c__;
    /* Function Body */
    if (*m <= 0 || *n <= 0)
    {
        *(unsigned char *)equed = 'N';
        return 0;
    }
    /* Initialize LARGE and SMALL. */
    small = dlamch_("Safe minimum") / dlamch_("Precision");
    large = 1. / small;
    if (*rowcnd >= .1 && *amax >= small && *amax <= large)
    {
        /* No row scaling */
        if (*colcnd >= .1)
        {
            /* No column scaling */
            *(unsigned char *)equed = 'N';
        }
        else
        {
            /* Column scaling */
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                cj = c__[j];
                i__2 = *m;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    a[i__ + j * a_dim1] = cj * a[i__ + j * a_dim1];
                    /* L10: */
                }
                /* L20: */
            }
            *(unsigned char *)equed = 'C';
        }
    }
    else if (*colcnd >= .1)
    {
        /* Row scaling, no column scaling */
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
                a[i__ + j * a_dim1] = r__[i__] * a[i__ + j * a_dim1];
                /* L30: */
            }
            /* L40: */
        }
        *(unsigned char *)equed = 'R';
    }
    else
    {
        /* Row and column scaling */
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            cj = c__[j];
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + j * a_dim1] = cj * r__[i__] * a[i__ + j * a_dim1];
                /* L50: */
            }
            /* L60: */
        }
        *(unsigned char *)equed = 'B';
    }
    return 0;
    /* End of DLAQGE */
}
/* dlaqge_ */
