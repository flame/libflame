/* ../netlib/clanht.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b CLANHT returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele ment of largest absolute value of a complex Hermitian tridiagonal matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLANHT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clanht. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clanht. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clanht. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* REAL FUNCTION CLANHT( NORM, N, D, E ) */
/* .. Scalar Arguments .. */
/* CHARACTER NORM */
/* INTEGER N */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ) */
/* COMPLEX E( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLANHT returns the value of the one norm, or the Frobenius norm, or */
/* > the infinity norm, or the element of largest absolute value of a */
/* > complex Hermitian tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \return CLANHT */
/* > \verbatim */
/* > */
/* > CLANHT = ( max(f2c_abs(A(i,j))), NORM = 'M' or 'm' */
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
/* > Specifies the value to be returned in CLANHT as described */
/* > above. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. When N = 0, CLANHT is */
/* > set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is COMPLEX array, dimension (N-1) */
/* > The (n-1) sub-diagonal or super-diagonal elements of A. */
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
real clanht_(char *norm, integer *n, real *d__, complex *e)
{
    /* System generated locals */
    integer i__1;
    real ret_val, r__1;
    /* Builtin functions */
    double c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    integer i__;
    real sum, scale;
    extern logical lsame_(char *, char *);
    real anorm;
    extern /* Subroutine */
    int classq_(integer *, complex *, integer *, real *, real *);
    extern logical sisnan_(real *);
    extern /* Subroutine */
    int slassq_(integer *, real *, integer *, real *, real *);
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
    --e;
    --d__;
    /* Function Body */
    if (*n <= 0)
    {
        anorm = 0.f;
    }
    else if (lsame_(norm, "M"))
    {
        /* Find max(f2c_abs(A(i,j))). */
        anorm = (r__1 = d__[*n], f2c_abs(r__1));
        i__1 = *n - 1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            sum = (r__1 = d__[i__], f2c_abs(r__1));
            if (anorm < sum || sisnan_(&sum))
            {
                anorm = sum;
            }
            sum = c_abs(&e[i__]);
            if (anorm < sum || sisnan_(&sum))
            {
                anorm = sum;
            }
            /* L10: */
        }
    }
    else if (lsame_(norm, "O") || *(unsigned char *) norm == '1' || lsame_(norm, "I"))
    {
        /* Find norm1(A). */
        if (*n == 1)
        {
            anorm = f2c_abs(d__[1]);
        }
        else
        {
            anorm = f2c_abs(d__[1]) + c_abs(&e[1]);
            sum = c_abs(&e[*n - 1]) + (r__1 = d__[*n], f2c_abs(r__1));
            if (anorm < sum || sisnan_(&sum))
            {
                anorm = sum;
            }
            i__1 = *n - 1;
            for (i__ = 2;
                    i__ <= i__1;
                    ++i__)
            {
                sum = (r__1 = d__[i__], f2c_abs(r__1)) + c_abs(&e[i__]) + c_abs(& e[i__ - 1]);
                if (anorm < sum || sisnan_(&sum))
                {
                    anorm = sum;
                }
                /* L20: */
            }
        }
    }
    else if (lsame_(norm, "F") || lsame_(norm, "E"))
    {
        /* Find normF(A). */
        scale = 0.f;
        sum = 1.f;
        if (*n > 1)
        {
            i__1 = *n - 1;
            classq_(&i__1, &e[1], &c__1, &scale, &sum);
            sum *= 2;
        }
        slassq_(n, &d__[1], &c__1, &scale, &sum);
        anorm = scale * sqrt(sum);
    }
    ret_val = anorm;
    return ret_val;
    /* End of CLANHT */
}
/* clanht_ */
