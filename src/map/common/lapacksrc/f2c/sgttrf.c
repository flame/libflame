/* ../netlib/sgttrf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SGTTRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgttrf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgttrf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgttrf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGTTRF( N, DL, D, DU, DU2, IPIV, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL D( * ), DL( * ), DU( * ), DU2( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGTTRF computes an LU factorization of a real tridiagonal matrix A */
/* > using elimination with partial pivoting and row interchanges. */
/* > */
/* > The factorization has the form */
/* > A = L * U */
/* > where L is a product of permutation and unit lower bidiagonal */
/* > matrices and U is upper triangular with nonzeros in only the main */
/* > diagonal and first two superdiagonals. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DL */
/* > \verbatim */
/* > DL is REAL array, dimension (N-1) */
/* > On entry, DL must contain the (n-1) sub-diagonal elements of */
/* > A. */
/* > */
/* > On exit, DL is overwritten by the (n-1) multipliers that */
/* > define the matrix L from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, D must contain the diagonal elements of A. */
/* > */
/* > On exit, D is overwritten by the n diagonal elements of the */
/* > upper triangular matrix U from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* > DU is REAL array, dimension (N-1) */
/* > On entry, DU must contain the (n-1) super-diagonal elements */
/* > of A. */
/* > */
/* > On exit, DU is overwritten by the (n-1) elements of the first */
/* > super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] DU2 */
/* > \verbatim */
/* > DU2 is REAL array, dimension (N-2) */
/* > On exit, DU2 is overwritten by the (n-2) elements of the */
/* > second super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices;
for 1 <= i <= n, row i of the matrix was */
/* > interchanged with row IPIV(i). IPIV(i) will always be either */
/* > i or i+1;
IPIV(i) = i indicates a row interchange was not */
/* > required. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -k, the k-th argument had an illegal value */
/* > > 0: if INFO = k, U(k,k) is exactly zero. The factorization */
/* > has been completed, but the factor U is exactly */
/* > singular, and division by zero will occur if it is used */
/* > to solve a system of equations. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realGTcomputational */
/* ===================================================================== */
/* Subroutine */
int sgttrf_(integer *n, real *dl, real *d__, real *du, real * du2, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Local variables */
    integer i__;
    real fact, temp;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    /* -- LAPACK computational routine (version 3.4.2) -- */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --ipiv;
    --du2;
    --du;
    --d__;
    --dl;
    /* Function Body */
    *info = 0;
    if (*n < 0)
    {
        *info = -1;
        i__1 = -(*info);
        xerbla_("SGTTRF", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Initialize IPIV(i) = i and DU2(I) = 0 */
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        ipiv[i__] = i__;
        /* L10: */
    }
    i__1 = *n - 2;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        du2[i__] = 0.f;
        /* L20: */
    }
    i__1 = *n - 2;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        if ((r__1 = d__[i__], f2c_abs(r__1)) >= (r__2 = dl[i__], f2c_abs(r__2)))
        {
            /* No row interchange required, eliminate DL(I) */
            if (d__[i__] != 0.f)
            {
                fact = dl[i__] / d__[i__];
                dl[i__] = fact;
                d__[i__ + 1] -= fact * du[i__];
            }
        }
        else
        {
            /* Interchange rows I and I+1, eliminate DL(I) */
            fact = d__[i__] / dl[i__];
            d__[i__] = dl[i__];
            dl[i__] = fact;
            temp = du[i__];
            du[i__] = d__[i__ + 1];
            d__[i__ + 1] = temp - fact * d__[i__ + 1];
            du2[i__] = du[i__ + 1];
            du[i__ + 1] = -fact * du[i__ + 1];
            ipiv[i__] = i__ + 1;
        }
        /* L30: */
    }
    if (*n > 1)
    {
        i__ = *n - 1;
        if ((r__1 = d__[i__], f2c_abs(r__1)) >= (r__2 = dl[i__], f2c_abs(r__2)))
        {
            if (d__[i__] != 0.f)
            {
                fact = dl[i__] / d__[i__];
                dl[i__] = fact;
                d__[i__ + 1] -= fact * du[i__];
            }
        }
        else
        {
            fact = d__[i__] / dl[i__];
            d__[i__] = dl[i__];
            dl[i__] = fact;
            temp = du[i__];
            du[i__] = d__[i__ + 1];
            d__[i__ + 1] = temp - fact * d__[i__ + 1];
            ipiv[i__] = i__ + 1;
        }
    }
    /* Check for a zero on the diagonal of U. */
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        if (d__[i__] == 0.f)
        {
            *info = i__;
            goto L50;
        }
        /* L40: */
    }
L50:
    return 0;
    /* End of SGTTRF */
}
/* sgttrf_ */
