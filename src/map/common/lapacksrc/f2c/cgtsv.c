/* ../netlib/cgtsv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief <b> CGTSV computes the solution to system of linear equations A * X = B for GT matrices <b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGTSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgtsv.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgtsv.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgtsv.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGTSV( N, NRHS, DL, D, DU, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX B( LDB, * ), D( * ), DL( * ), DU( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGTSV solves the equation */
/* > */
/* > A*X = B, */
/* > */
/* > where A is an N-by-N tridiagonal matrix, by Gaussian elimination with */
/* > partial pivoting. */
/* > */
/* > Note that the equation A**T *X = B may be solved by interchanging the */
/* > order of the arguments DU and DL. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DL */
/* > \verbatim */
/* > DL is COMPLEX array, dimension (N-1) */
/* > On entry, DL must contain the (n-1) subdiagonal elements of */
/* > A. */
/* > On exit, DL is overwritten by the (n-2) elements of the */
/* > second superdiagonal of the upper triangular matrix U from */
/* > the LU factorization of A, in DL(1), ..., DL(n-2). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is COMPLEX array, dimension (N) */
/* > On entry, D must contain the diagonal elements of A. */
/* > On exit, D is overwritten by the n diagonal elements of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* > DU is COMPLEX array, dimension (N-1) */
/* > On entry, DU must contain the (n-1) superdiagonal elements */
/* > of A. */
/* > On exit, DU is overwritten by the (n-1) elements of the first */
/* > superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > On entry, the N-by-NRHS right hand side matrix B. */
/* > On exit, if INFO = 0, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, U(i,i) is exactly zero, and the solution */
/* > has not been computed. The factorization has not been */
/* > completed unless i = N. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexGTsolve */
/* ===================================================================== */
/* Subroutine */
int cgtsv_(integer *n, integer *nrhs, complex *dl, complex * d__, complex *du, complex *b, integer *ldb, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3, q__4, q__5;
    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    integer j, k;
    complex temp, mult;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    /* -- LAPACK driver routine (version 3.4.2) -- */
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    *info = 0;
    if (*n < 0)
    {
        *info = -1;
    }
    else if (*nrhs < 0)
    {
        *info = -2;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -7;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGTSV ", &i__1);
        return 0;
    }
    if (*n == 0)
    {
        return 0;
    }
    i__1 = *n - 1;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        i__2 = k;
        if (dl[i__2].r == 0.f && dl[i__2].i == 0.f)
        {
            /* Subdiagonal is zero, no elimination is required. */
            i__2 = k;
            if (d__[i__2].r == 0.f && d__[i__2].i == 0.f)
            {
                /* Diagonal is zero: set INFO = K and return;
                a unique */
                /* solution can not be found. */
                *info = k;
                return 0;
            }
        }
        else /* if(complicated condition) */
        {
            i__2 = k;
            i__3 = k;
            if ((r__1 = d__[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&d__[k]), f2c_abs(r__2)) >= (r__3 = dl[i__3].r, f2c_abs(r__3)) + (r__4 = r_imag(&dl[k]), f2c_abs(r__4)))
            {
                /* No row interchange required */
                c_div(&q__1, &dl[k], &d__[k]);
                mult.r = q__1.r;
                mult.i = q__1.i; // , expr subst
                i__2 = k + 1;
                i__3 = k + 1;
                i__4 = k;
                q__2.r = mult.r * du[i__4].r - mult.i * du[i__4].i;
                q__2.i = mult.r * du[i__4].i + mult.i * du[i__4].r; // , expr subst
                q__1.r = d__[i__3].r - q__2.r;
                q__1.i = d__[i__3].i - q__2.i; // , expr subst
                d__[i__2].r = q__1.r;
                d__[i__2].i = q__1.i; // , expr subst
                i__2 = *nrhs;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = k + 1 + j * b_dim1;
                    i__4 = k + 1 + j * b_dim1;
                    i__5 = k + j * b_dim1;
                    q__2.r = mult.r * b[i__5].r - mult.i * b[i__5].i;
                    q__2.i = mult.r * b[i__5].i + mult.i * b[i__5].r; // , expr subst
                    q__1.r = b[i__4].r - q__2.r;
                    q__1.i = b[i__4].i - q__2.i; // , expr subst
                    b[i__3].r = q__1.r;
                    b[i__3].i = q__1.i; // , expr subst
                    /* L10: */
                }
                if (k < *n - 1)
                {
                    i__2 = k;
                    dl[i__2].r = 0.f;
                    dl[i__2].i = 0.f; // , expr subst
                }
            }
            else
            {
                /* Interchange rows K and K+1 */
                c_div(&q__1, &d__[k], &dl[k]);
                mult.r = q__1.r;
                mult.i = q__1.i; // , expr subst
                i__2 = k;
                i__3 = k;
                d__[i__2].r = dl[i__3].r;
                d__[i__2].i = dl[i__3].i; // , expr subst
                i__2 = k + 1;
                temp.r = d__[i__2].r;
                temp.i = d__[i__2].i; // , expr subst
                i__2 = k + 1;
                i__3 = k;
                q__2.r = mult.r * temp.r - mult.i * temp.i;
                q__2.i = mult.r * temp.i + mult.i * temp.r; // , expr subst
                q__1.r = du[i__3].r - q__2.r;
                q__1.i = du[i__3].i - q__2.i; // , expr subst
                d__[i__2].r = q__1.r;
                d__[i__2].i = q__1.i; // , expr subst
                if (k < *n - 1)
                {
                    i__2 = k;
                    i__3 = k + 1;
                    dl[i__2].r = du[i__3].r;
                    dl[i__2].i = du[i__3].i; // , expr subst
                    i__2 = k + 1;
                    q__2.r = -mult.r;
                    q__2.i = -mult.i; // , expr subst
                    i__3 = k;
                    q__1.r = q__2.r * dl[i__3].r - q__2.i * dl[i__3].i;
                    q__1.i = q__2.r * dl[i__3].i + q__2.i * dl[i__3] .r; // , expr subst
                    du[i__2].r = q__1.r;
                    du[i__2].i = q__1.i; // , expr subst
                }
                i__2 = k;
                du[i__2].r = temp.r;
                du[i__2].i = temp.i; // , expr subst
                i__2 = *nrhs;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = k + j * b_dim1;
                    temp.r = b[i__3].r;
                    temp.i = b[i__3].i; // , expr subst
                    i__3 = k + j * b_dim1;
                    i__4 = k + 1 + j * b_dim1;
                    b[i__3].r = b[i__4].r;
                    b[i__3].i = b[i__4].i; // , expr subst
                    i__3 = k + 1 + j * b_dim1;
                    i__4 = k + 1 + j * b_dim1;
                    q__2.r = mult.r * b[i__4].r - mult.i * b[i__4].i;
                    q__2.i = mult.r * b[i__4].i + mult.i * b[i__4].r; // , expr subst
                    q__1.r = temp.r - q__2.r;
                    q__1.i = temp.i - q__2.i; // , expr subst
                    b[i__3].r = q__1.r;
                    b[i__3].i = q__1.i; // , expr subst
                    /* L20: */
                }
            }
        }
        /* L30: */
    }
    i__1 = *n;
    if (d__[i__1].r == 0.f && d__[i__1].i == 0.f)
    {
        *info = *n;
        return 0;
    }
    /* Back solve with the matrix U from the factorization. */
    i__1 = *nrhs;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        i__2 = *n + j * b_dim1;
        c_div(&q__1, &b[*n + j * b_dim1], &d__[*n]);
        b[i__2].r = q__1.r;
        b[i__2].i = q__1.i; // , expr subst
        if (*n > 1)
        {
            i__2 = *n - 1 + j * b_dim1;
            i__3 = *n - 1 + j * b_dim1;
            i__4 = *n - 1;
            i__5 = *n + j * b_dim1;
            q__3.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i;
            q__3.i = du[i__4].r * b[i__5].i + du[i__4].i * b[i__5].r; // , expr subst
            q__2.r = b[i__3].r - q__3.r;
            q__2.i = b[i__3].i - q__3.i; // , expr subst
            c_div(&q__1, &q__2, &d__[*n - 1]);
            b[i__2].r = q__1.r;
            b[i__2].i = q__1.i; // , expr subst
        }
        for (k = *n - 2;
                k >= 1;
                --k)
        {
            i__2 = k + j * b_dim1;
            i__3 = k + j * b_dim1;
            i__4 = k;
            i__5 = k + 1 + j * b_dim1;
            q__4.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i;
            q__4.i = du[i__4].r * b[i__5].i + du[i__4].i * b[i__5].r; // , expr subst
            q__3.r = b[i__3].r - q__4.r;
            q__3.i = b[i__3].i - q__4.i; // , expr subst
            i__6 = k;
            i__7 = k + 2 + j * b_dim1;
            q__5.r = dl[i__6].r * b[i__7].r - dl[i__6].i * b[i__7].i;
            q__5.i = dl[i__6].r * b[i__7].i + dl[i__6].i * b[i__7].r; // , expr subst
            q__2.r = q__3.r - q__5.r;
            q__2.i = q__3.i - q__5.i; // , expr subst
            c_div(&q__1, &q__2, &d__[k]);
            b[i__2].r = q__1.r;
            b[i__2].i = q__1.i; // , expr subst
            /* L40: */
        }
        /* L50: */
    }
    return 0;
    /* End of CGTSV */
}
/* cgtsv_ */
