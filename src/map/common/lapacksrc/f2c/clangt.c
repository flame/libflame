/* ../netlib/clangt.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CLANGT returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general tridiagonal matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLANGT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clangt. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clangt. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clangt. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION CLANGT( NORM, N, DL, D, DU ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM */
/* INTEGER N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX D( * ), DL( * ), DU( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANGT returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > complex tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANGT */
/* > \verbatim */
/* > */
/* > CLANGT = ( max(f2c_abs(A(i,j))), NORM = 'M' or 'm' */
/* > ( */
/* > ( norm1(A), NORM = '1', 'O' or 'o' */
/* > ( */
/* > ( normI(A), NORM = 'I' or 'i' */
/* > ( */
/* > ( normF(A), NORM = 'F', 'f', 'E' or 'e' */
/* > */
/* > where norm1 denotes the one norm of a matrix (maximum column sum), */
/* > normI denotes the infinity norm of a matrix (maximum row sum) and */
/* > normF denotes the Frobenius norm of a matrix (square root of sum of */
/* > squares). Note that max(f2c_abs(A(i,j))) is not a consistent matrix norm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] NORM */
/* > \verbatim */
/* > NORM is CHARACTER*1 */
/* > Specifies the value to be returned in CLANGT as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. When N = 0, CLANGT is */
/* > set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* > DL is COMPLEX array, dimension (N-1) */
/* > The (n-1) sub-diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is COMPLEX array, dimension (N) */
/* > The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* > DU is COMPLEX array, dimension (N-1) */
/* > The (n-1) super-diagonal elements of A. */
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
real clangt_(char *norm, integer *n, complex *dl, complex *d__, complex *du)
{
    /* System generated locals */
    integer i__1;
    real ret_val, r__1;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    integer i__;
    real sum, temp, scale;
    extern logical lsame_(char *, char *);
    real anorm;
    extern /* Subroutine */
    int classq_(integer *, complex *, integer *, real *, real *);
    extern logical sisnan_(real *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --du;
    --d__;
    --dl;
    /* Function Body */
    if (*n <= 0)
    {
        anorm = 0.f;
    }
    else if (lsame_(norm, "M"))
    {
        /* Find max(f2c_abs(A(i,j))). */
        anorm = c_abs(&d__[*n]);
        i__1 = *n - 1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            r__1 = c_abs(&dl[i__]);
            if (anorm < c_abs(&dl[i__]) || sisnan_(&r__1))
            {
                anorm = c_abs(&dl[i__]);
            }
            r__1 = c_abs(&d__[i__]);
            if (anorm < c_abs(&d__[i__]) || sisnan_(&r__1))
            {
                anorm = c_abs(&d__[i__]);
            }
            r__1 = c_abs(&du[i__]);
            if (anorm < c_abs(&du[i__]) || sisnan_(&r__1))
            {
                anorm = c_abs(&du[i__]);
            }
            /* L10: */
        }
    }
    else if (lsame_(norm, "O") || *(unsigned char *) norm == '1')
    {
        /* Find norm1(A). */
        if (*n == 1)
        {
            anorm = c_abs(&d__[1]);
        }
        else
        {
            anorm = c_abs(&d__[1]) + c_abs(&dl[1]);
            temp = c_abs(&d__[*n]) + c_abs(&du[*n - 1]);
            if (anorm < temp || sisnan_(&temp))
            {
                anorm = temp;
            }
            i__1 = *n - 1;
            for (i__ = 2;
                    i__ <= i__1;
                    ++i__)
            {
                temp = c_abs(&d__[i__]) + c_abs(&dl[i__]) + c_abs(&du[i__ - 1] );
                if (anorm < temp || sisnan_(&temp))
                {
                    anorm = temp;
                }
                /* L20: */
            }
        }
    }
    else if (lsame_(norm, "I"))
    {
        /* Find normI(A). */
        if (*n == 1)
        {
            anorm = c_abs(&d__[1]);
        }
        else
        {
            anorm = c_abs(&d__[1]) + c_abs(&du[1]);
            temp = c_abs(&d__[*n]) + c_abs(&dl[*n - 1]);
            if (anorm < temp || sisnan_(&temp))
            {
                anorm = temp;
            }
            i__1 = *n - 1;
            for (i__ = 2;
                    i__ <= i__1;
                    ++i__)
            {
                temp = c_abs(&d__[i__]) + c_abs(&du[i__]) + c_abs(&dl[i__ - 1] );
                if (anorm < temp || sisnan_(&temp))
                {
                    anorm = temp;
                }
                /* L30: */
            }
        }
    }
    else if (lsame_(norm, "F") || lsame_(norm, "E"))
    {
        /* Find normF(A). */
        scale = 0.f;
        sum = 1.f;
        classq_(n, &d__[1], &c__1, &scale, &sum);
        if (*n > 1)
        {
            i__1 = *n - 1;
            classq_(&i__1, &dl[1], &c__1, &scale, &sum);
            i__1 = *n - 1;
            classq_(&i__1, &du[1], &c__1, &scale, &sum);
        }
        anorm = scale * sqrt(sum);
    }
    ret_val = anorm;
    return ret_val;
    /* End of CLANGT */
}
/* clangt_ */
