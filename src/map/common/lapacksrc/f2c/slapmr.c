/* ../netlib/slapmr.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLAPMR rearranges rows of a matrix as specified by a permutation vector. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAPMR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapmr. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapmr. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapmr. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAPMR( FORWRD, M, N, X, LDX, K ) */
/* .. Scalar Arguments .. */
/* LOGICAL FORWRD */
/* INTEGER LDX, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER K( * ) */
/* REAL X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAPMR rearranges the rows of the M by N matrix X as specified */
/* > by the permutation K(1),K(2),...,K(M) of the integers 1,...,M. */
/* > If FORWRD = .TRUE., forward permutation: */
/* > */
/* > X(K(I),*) is moved X(I,*) for I = 1,2,...,M. */
/* > */
/* > If FORWRD = .FALSE., backward permutation: */
/* > */
/* > X(I,*) is moved to X(K(I),*) for I = 1,2,...,M. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] FORWRD */
/* > \verbatim */
/* > FORWRD is LOGICAL */
/* > = .TRUE., forward permutation */
/* > = .FALSE., backward permutation */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix X. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix X. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is REAL array, dimension (LDX,N) */
/* > On entry, the M by N matrix X. */
/* > On exit, X contains the permuted matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X, LDX >= MAX(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] K */
/* > \verbatim */
/* > K is INTEGER array, dimension (M) */
/* > On entry, K contains the permutation vector. K is used as */
/* > internal workspace, but reset to its original value on */
/* > output. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
int slapmr_(logical *forwrd, integer *m, integer *n, real *x, integer *ldx, integer *k)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;
    /* Local variables */
    integer i__, j, jj, in;
    real temp;
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
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --k;
    /* Function Body */
    if (*m <= 1)
    {
        return 0;
    }
    i__1 = *m;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        k[i__] = -k[i__];
        /* L10: */
    }
    if (*forwrd)
    {
        /* Forward permutation */
        i__1 = *m;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            if (k[i__] > 0)
            {
                goto L40;
            }
            j = i__;
            k[j] = -k[j];
            in = k[j];
L20:
            if (k[in] > 0)
            {
                goto L40;
            }
            i__2 = *n;
            for (jj = 1;
                    jj <= i__2;
                    ++jj)
            {
                temp = x[j + jj * x_dim1];
                x[j + jj * x_dim1] = x[in + jj * x_dim1];
                x[in + jj * x_dim1] = temp;
                /* L30: */
            }
            k[in] = -k[in];
            j = in;
            in = k[in];
            goto L20;
L40: /* L50: */
            ;
        }
    }
    else
    {
        /* Backward permutation */
        i__1 = *m;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            if (k[i__] > 0)
            {
                goto L80;
            }
            k[i__] = -k[i__];
            j = k[i__];
L60:
            if (j == i__)
            {
                goto L80;
            }
            i__2 = *n;
            for (jj = 1;
                    jj <= i__2;
                    ++jj)
            {
                temp = x[i__ + jj * x_dim1];
                x[i__ + jj * x_dim1] = x[j + jj * x_dim1];
                x[j + jj * x_dim1] = temp;
                /* L70: */
            }
            k[j] = -k[j];
            j = k[j];
            goto L60;
L80: /* L90: */
            ;
        }
    }
    return 0;
    /* End of ZLAPMT */
}
/* slapmr_ */
