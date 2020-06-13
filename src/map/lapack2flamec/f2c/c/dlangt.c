/* ../netlib/dlangt.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DLANGT returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general tridiagonal matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLANGT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlangt. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlangt. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlangt. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* DOUBLE PRECISION FUNCTION DLANGT( NORM, N, DL, D, DU ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM */
/* INTEGER N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION D( * ), DL( * ), DU( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANGT returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > real tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \return DLANGT */
/* > \verbatim */
/* > */
/* > DLANGT = ( max(f2c_dabs(A(i,j))), NORM = 'M' or 'm' */
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
/* > squares). Note that max(f2c_dabs(A(i,j))) is not a consistent matrix norm. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] NORM */
/* > \verbatim */
/* > NORM is CHARACTER*1 */
/* > Specifies the value to be returned in DLANGT as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. When N = 0, DLANGT is */
/* > set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* > DL is DOUBLE PRECISION array, dimension (N-1) */
/* > The (n-1) sub-diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* > DU is DOUBLE PRECISION array, dimension (N-1) */
/* > The (n-1) super-diagonal elements of A. */
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
doublereal dlangt_(char *norm, integer *n, doublereal *dl, doublereal *d__, doublereal *du)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__;
    doublereal sum, temp, scale;
    extern logical lsame_(char *, char *);
    doublereal anorm;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */
    int dlassq_(integer *, doublereal *, integer *, doublereal *, doublereal *);
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
        anorm = 0.;
    }
    else if (lsame_(norm, "M"))
    {
        /* Find max(f2c_dabs(A(i,j))). */
        anorm = (d__1 = d__[*n], f2c_dabs(d__1));
        i__1 = *n - 1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            d__3 = (d__2 = dl[i__], f2c_dabs(d__2));
            if (anorm < (d__1 = dl[i__], f2c_dabs(d__1)) || disnan_(&d__3))
            {
                anorm = (d__4 = dl[i__], f2c_dabs(d__4));
            }
            d__3 = (d__2 = d__[i__], f2c_dabs(d__2));
            if (anorm < (d__1 = d__[i__], f2c_dabs(d__1)) || disnan_(&d__3))
            {
                anorm = (d__4 = d__[i__], f2c_dabs(d__4));
            }
            d__3 = (d__2 = du[i__], f2c_dabs(d__2));
            if (anorm < (d__1 = du[i__], f2c_dabs(d__1)) || disnan_(&d__3))
            {
                anorm = (d__4 = du[i__], f2c_dabs(d__4));
            }
            /* L10: */
        }
    }
    else if (lsame_(norm, "O") || *(unsigned char *) norm == '1')
    {
        /* Find norm1(A). */
        if (*n == 1)
        {
            anorm = f2c_dabs(d__[1]);
        }
        else
        {
            anorm = f2c_dabs(d__[1]) + f2c_dabs(dl[1]);
            temp = (d__1 = d__[*n], f2c_dabs(d__1)) + (d__2 = du[*n - 1], f2c_dabs(d__2) );
            if (anorm < temp || disnan_(&temp))
            {
                anorm = temp;
            }
            i__1 = *n - 1;
            for (i__ = 2;
                    i__ <= i__1;
                    ++i__)
            {
                temp = (d__1 = d__[i__], f2c_dabs(d__1)) + (d__2 = dl[i__], f2c_dabs( d__2)) + (d__3 = du[i__ - 1], f2c_dabs(d__3));
                if (anorm < temp || disnan_(&temp))
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
            anorm = f2c_dabs(d__[1]);
        }
        else
        {
            anorm = f2c_dabs(d__[1]) + f2c_dabs(du[1]);
            temp = (d__1 = d__[*n], f2c_dabs(d__1)) + (d__2 = dl[*n - 1], f2c_dabs(d__2) );
            if (anorm < temp || disnan_(&temp))
            {
                anorm = temp;
            }
            i__1 = *n - 1;
            for (i__ = 2;
                    i__ <= i__1;
                    ++i__)
            {
                temp = (d__1 = d__[i__], f2c_dabs(d__1)) + (d__2 = du[i__], f2c_dabs( d__2)) + (d__3 = dl[i__ - 1], f2c_dabs(d__3));
                if (anorm < temp || disnan_(&temp))
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
        scale = 0.;
        sum = 1.;
        dlassq_(n, &d__[1], &c__1, &scale, &sum);
        if (*n > 1)
        {
            i__1 = *n - 1;
            dlassq_(&i__1, &dl[1], &c__1, &scale, &sum);
            i__1 = *n - 1;
            dlassq_(&i__1, &du[1], &c__1, &scale, &sum);
        }
        anorm = scale * sqrt(sum);
    }
    ret_val = anorm;
    return ret_val;
    /* End of DLANGT */
}
/* dlangt_ */
