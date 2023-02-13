/* ../netlib/claqgb.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLAQGB scales a general band matrix, using row and column scaling factors computed by sgbequ. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAQGB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqgb. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqgb. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqgb. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAQGB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/* AMAX, EQUED ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED */
/* INTEGER KL, KU, LDAB, M, N */
/* REAL AMAX, COLCND, ROWCND */
/* .. */
/* .. Array Arguments .. */
/* REAL C( * ), R( * ) */
/* COMPLEX AB( LDAB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQGB equilibrates a general M by N band matrix A with KL */
/* > subdiagonals and KU superdiagonals using the row and scaling factors */
/* > in the vectors R and C. */
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
/* > \param[in] KL */
/* > \verbatim */
/* > KL is INTEGER */
/* > The number of subdiagonals within the band of A. KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* > KU is INTEGER */
/* > The number of superdiagonals within the band of A. KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB,N) */
/* > On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* > The j-th column of A is stored in the j-th column of the */
/* > array AB as follows: */
/* > AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl) */
/* > */
/* > On exit, the equilibrated matrix, in the same storage format */
/* > as A. See EQUED for the form of the equilibrated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDA >= KL+KU+1. */
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
/* > \ingroup complexGBauxiliary */
/* ===================================================================== */
/* Subroutine */
int claqgb_(integer *m, integer *n, integer *kl, integer *ku, complex *ab, integer *ldab, real *r__, real *c__, real *rowcnd, real *colcnd, real *amax, char *equed)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1;
    complex q__1;
    /* Local variables */
    integer i__, j;
    real cj, large, small;
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --r__;
    --c__;
    /* Function Body */
    if (*m <= 0 || *n <= 0)
    {
        *(unsigned char *)equed = 'N';
        return 0;
    }
    /* Initialize LARGE and SMALL. */
    small = slamch_("Safe minimum") / slamch_("Precision");
    large = 1.f / small;
    if (*rowcnd >= .1f && *amax >= small && *amax <= large)
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
                /* Computing MAX */
                i__2 = 1;
                i__3 = j - *ku; // , expr subst
                /* Computing MIN */
                i__5 = *m;
                i__6 = j + *kl; // , expr subst
                i__4 = min(i__5,i__6);
                for (i__ = max(i__2,i__3);
                        i__ <= i__4;
                        ++i__)
                {
                    i__2 = *ku + 1 + i__ - j + j * ab_dim1;
                    i__3 = *ku + 1 + i__ - j + j * ab_dim1;
                    q__1.r = cj * ab[i__3].r;
                    q__1.i = cj * ab[i__3].i; // , expr subst
                    ab[i__2].r = q__1.r;
                    ab[i__2].i = q__1.i; // , expr subst
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
            /* Computing MAX */
            i__4 = 1;
            i__2 = j - *ku; // , expr subst
            /* Computing MIN */
            i__5 = *m;
            i__6 = j + *kl; // , expr subst
            i__3 = min(i__5,i__6);
            for (i__ = max(i__4,i__2);
                    i__ <= i__3;
                    ++i__)
            {
                i__4 = *ku + 1 + i__ - j + j * ab_dim1;
                i__2 = i__;
                i__5 = *ku + 1 + i__ - j + j * ab_dim1;
                q__1.r = r__[i__2] * ab[i__5].r;
                q__1.i = r__[i__2] * ab[i__5] .i; // , expr subst
                ab[i__4].r = q__1.r;
                ab[i__4].i = q__1.i; // , expr subst
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
            /* Computing MAX */
            i__3 = 1;
            i__4 = j - *ku; // , expr subst
            /* Computing MIN */
            i__5 = *m;
            i__6 = j + *kl; // , expr subst
            i__2 = min(i__5,i__6);
            for (i__ = max(i__3,i__4);
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = *ku + 1 + i__ - j + j * ab_dim1;
                r__1 = cj * r__[i__];
                i__4 = *ku + 1 + i__ - j + j * ab_dim1;
                q__1.r = r__1 * ab[i__4].r;
                q__1.i = r__1 * ab[i__4].i; // , expr subst
                ab[i__3].r = q__1.r;
                ab[i__3].i = q__1.i; // , expr subst
                /* L50: */
            }
            /* L60: */
        }
        *(unsigned char *)equed = 'B';
    }
    return 0;
    /* End of CLAQGB */
}
/* claqgb_ */
