/* ../netlib/claqge.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLAQGE scales a general rectangular matrix, using row and column scaling factors computed by sg eequ. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAQGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqge. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqge. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqge. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, */
/* EQUED ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED */
/* INTEGER LDA, M, N */
/* REAL AMAX, COLCND, ROWCND */
/* .. */
/* .. Array Arguments .. */
/* REAL C( * ), R( * ) */
/* COMPLEX A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQGE equilibrates a general M by N matrix A using the row and */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the M by N matrix A. */
/* > On exit, the equilibrated matrix. See EQUED for the form of */
/* > the equilibrated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(M,1). */
/* > \endverbatim */
/* > */
/* > \param[in] R */
/* > \verbatim */
/* > R is REAL array, dimension (M) */
/* > The row scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is REAL array, dimension (N) */
/* > The column scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] ROWCND */
/* > \verbatim */
/* > ROWCND is REAL */
/* > Ratio of the smallest R(i) to the largest R(i). */
/* > \endverbatim */
/* > */
/* > \param[in] COLCND */
/* > \verbatim */
/* > COLCND is REAL */
/* > Ratio of the smallest C(i) to the largest C(i). */
/* > \endverbatim */
/* > */
/* > \param[in] AMAX */
/* > \verbatim */
/* > AMAX is REAL */
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
/* > \ingroup complexGEauxiliary */
/* ===================================================================== */
/* Subroutine */
int claqge_(integer *m, integer *n, complex *a, integer *lda, real *r__, real *c__, real *rowcnd, real *colcnd, real *amax, char * equed)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"claqge inputs: m %lld, n %lld, lda %lld",*m, *n, *lda);
#else
    snprintf(buffer, 256,"claqge inputs: m %d, n %d, lda %d",*m, *n, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1;
    complex q__1;
    /* Local variables */
    integer i__, j;
    real cj, large, small_val;
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
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Initialize LARGE and SMALL. */
    small_val = slamch_("Safe minimum") / slamch_("Precision");
    large = 1.f / small_val;
    if (*rowcnd >= .1f && *amax >= small_val && *amax <= large)
    {
        /* No row scaling */
        if (*colcnd >= .1f)
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
                    i__3 = i__ + j * a_dim1;
                    i__4 = i__ + j * a_dim1;
                    q__1.r = cj * a[i__4].r;
                    q__1.i = cj * a[i__4].i; // , expr subst
                    a[i__3].r = q__1.r;
                    a[i__3].i = q__1.i; // , expr subst
                    /* L10: */
                }
                /* L20: */
            }
            *(unsigned char *)equed = 'C';
        }
    }
    else if (*colcnd >= .1f)
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
                i__3 = i__ + j * a_dim1;
                i__4 = i__;
                i__5 = i__ + j * a_dim1;
                q__1.r = r__[i__4] * a[i__5].r;
                q__1.i = r__[i__4] * a[i__5] .i; // , expr subst
                a[i__3].r = q__1.r;
                a[i__3].i = q__1.i; // , expr subst
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
                i__3 = i__ + j * a_dim1;
                r__1 = cj * r__[i__];
                i__4 = i__ + j * a_dim1;
                q__1.r = r__1 * a[i__4].r;
                q__1.i = r__1 * a[i__4].i; // , expr subst
                a[i__3].r = q__1.r;
                a[i__3].i = q__1.i; // , expr subst
                /* L50: */
            }
            /* L60: */
        }
        *(unsigned char *)equed = 'B';
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CLAQGE */
}
/* claqge_ */
