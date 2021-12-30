/* ../netlib/dlaqsp.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLAQSP scales a symmetric/Hermitian matrix in packed storage, using scaling factors computed by sppequ. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAQSP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqsp. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqsp. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqsp. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAQSP( UPLO, N, AP, S, SCOND, AMAX, EQUED ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED, UPLO */
/* INTEGER N */
/* DOUBLE PRECISION AMAX, SCOND */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION AP( * ), S( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQSP equilibrates a symmetric matrix A using the scaling factors */
/* > in the vector S. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > symmetric matrix A is stored. */
/* > = 'U': Upper triangular */
/* > = 'L': Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the symmetric matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* > On exit, the equilibrated matrix: diag(S) * A * diag(S), in */
/* > the same storage format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is DOUBLE PRECISION array, dimension (N) */
/* > The scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[in] SCOND */
/* > \verbatim */
/* > SCOND is DOUBLE PRECISION */
/* > Ratio of the smallest S(i) to the largest S(i). */
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
/* > Specifies whether or not equilibration was done. */
/* > = 'N': No equilibration. */
/* > = 'Y': Equilibration was done, i.e., A has been replaced by */
/* > diag(S) * A * diag(S). */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > THRESH is a threshold value used to decide if scaling should be done */
/* > based on the ratio of the scaling factors. If SCOND < THRESH, */
/* > scaling is done. */
/* > */
/* > LARGE and SMALL are threshold values used to decide if scaling should */
/* > be done based on the absolute size of the largest matrix element. */
/* > If AMAX > LARGE or AMAX < SMALL, scaling is done. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int dlaqsp_(char *uplo, integer *n, doublereal *ap, doublereal *s, doublereal *scond, doublereal *amax, char *equed)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE 
    char buffer[256]; 
    snprintf(buffer, 256,"dlaqsp inputs: uplo %c, n %" FLA_IS "",*uplo, *n);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    integer i__, j, jc;
    doublereal cj, large;
    extern logical lsame_(char *, char *);
    doublereal small;
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
    --s;
    --ap;
    /* Function Body */
    if (*n <= 0)
    {
        *(unsigned char *)equed = 'N';
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Initialize LARGE and SMALL. */
    small = dlamch_("Safe minimum") / dlamch_("Precision");
    large = 1. / small;
    if (*scond >= .1 && *amax >= small && *amax <= large)
    {
        /* No equilibration */
        *(unsigned char *)equed = 'N';
    }
    else
    {
        /* Replace A by diag(S) * A * diag(S). */
        if (lsame_(uplo, "U"))
        {
            /* Upper triangle of A is stored. */
            jc = 1;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                cj = s[j];
                i__2 = j;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    ap[jc + i__ - 1] = cj * s[i__] * ap[jc + i__ - 1];
                    /* L10: */
                }
                jc += j;
                /* L20: */
            }
        }
        else
        {
            /* Lower triangle of A is stored. */
            jc = 1;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                cj = s[j];
                i__2 = *n;
                for (i__ = j;
                        i__ <= i__2;
                        ++i__)
                {
                    ap[jc + i__ - j] = cj * s[i__] * ap[jc + i__ - j];
                    /* L30: */
                }
                jc = jc + *n - j + 1;
                /* L40: */
            }
        }
        *(unsigned char *)equed = 'Y';
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of DLAQSP */
}
/* dlaqsp_ */
