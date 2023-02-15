/* ../netlib/zgtts2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZGTTS2 solves a system of linear equations with a tridiagonal matrix using the LU factorization computed by sgttrf. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGTTS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgtts2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgtts2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgtts2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB ) */
/* .. Scalar Arguments .. */
/* INTEGER ITRANS, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGTTS2 solves one of the systems of equations */
/* > A * X = B, A**T * X = B, or A**H * X = B, */
/* > with a tridiagonal matrix A using the LU factorization computed */
/* > by ZGTTRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ITRANS */
/* > \verbatim */
/* > ITRANS is INTEGER */
/* > Specifies the form of the system of equations. */
/* > = 0: A * X = B (No transpose) */
/* > = 1: A**T * X = B (Transpose) */
/* > = 2: A**H * X = B (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* > DL is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) multipliers that define the matrix L from the */
/* > LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* > D is COMPLEX*16 array, dimension (N) */
/* > The n diagonal elements of the upper triangular matrix U from */
/* > the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* > DU is COMPLEX*16 array, dimension (N-1) */
/* > The (n-1) elements of the first super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* > DU2 is COMPLEX*16 array, dimension (N-2) */
/* > The (n-2) elements of the second super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
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
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > On entry, the matrix of right hand side vectors B. */
/* > On exit, B is overwritten by the solution vectors X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16GTcomputational */
/* ===================================================================== */
/* Subroutine */
int zgtts2_(integer *itrans, integer *n, integer *nrhs, doublecomplex *dl, doublecomplex *d__, doublecomplex *du, doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8;
    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg( doublecomplex *, doublecomplex *);
    /* Local variables */
    integer i__, j;
    doublecomplex temp;
    /* -- LAPACK computational routine (version 3.4.2) -- */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    --du2;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    if (*n == 0 || *nrhs == 0)
    {
        return 0;
    }
    if (*itrans == 0)
    {
        /* Solve A*X = B using the LU factorization of A, */
        /* overwriting each right hand side vector with its solution. */
        if (*nrhs <= 1)
        {
            j = 1;
L10: /* Solve L*x = b. */
            i__1 = *n - 1;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                if (ipiv[i__] == i__)
                {
                    i__2 = i__ + 1 + j * b_dim1;
                    i__3 = i__ + 1 + j * b_dim1;
                    i__4 = i__;
                    i__5 = i__ + j * b_dim1;
                    z__2.r = dl[i__4].r * b[i__5].r - dl[i__4].i * b[i__5].i;
                    z__2.i = dl[i__4].r * b[i__5].i + dl[i__4].i * b[ i__5].r; // , expr subst
                    z__1.r = b[i__3].r - z__2.r;
                    z__1.i = b[i__3].i - z__2.i; // , expr subst
                    b[i__2].r = z__1.r;
                    b[i__2].i = z__1.i; // , expr subst
                }
                else
                {
                    i__2 = i__ + j * b_dim1;
                    temp.r = b[i__2].r;
                    temp.i = b[i__2].i; // , expr subst
                    i__2 = i__ + j * b_dim1;
                    i__3 = i__ + 1 + j * b_dim1;
                    b[i__2].r = b[i__3].r;
                    b[i__2].i = b[i__3].i; // , expr subst
                    i__2 = i__ + 1 + j * b_dim1;
                    i__3 = i__;
                    i__4 = i__ + j * b_dim1;
                    z__2.r = dl[i__3].r * b[i__4].r - dl[i__3].i * b[i__4].i;
                    z__2.i = dl[i__3].r * b[i__4].i + dl[i__3].i * b[ i__4].r; // , expr subst
                    z__1.r = temp.r - z__2.r;
                    z__1.i = temp.i - z__2.i; // , expr subst
                    b[i__2].r = z__1.r;
                    b[i__2].i = z__1.i; // , expr subst
                }
                /* L20: */
            }
            /* Solve U*x = b. */
            i__1 = *n + j * b_dim1;
            z_div(&z__1, &b[*n + j * b_dim1], &d__[*n]);
            b[i__1].r = z__1.r;
            b[i__1].i = z__1.i; // , expr subst
            if (*n > 1)
            {
                i__1 = *n - 1 + j * b_dim1;
                i__2 = *n - 1 + j * b_dim1;
                i__3 = *n - 1;
                i__4 = *n + j * b_dim1;
                z__3.r = du[i__3].r * b[i__4].r - du[i__3].i * b[i__4].i;
                z__3.i = du[i__3].r * b[i__4].i + du[i__3].i * b[i__4] .r; // , expr subst
                z__2.r = b[i__2].r - z__3.r;
                z__2.i = b[i__2].i - z__3.i; // , expr subst
                z_div(&z__1, &z__2, &d__[*n - 1]);
                b[i__1].r = z__1.r;
                b[i__1].i = z__1.i; // , expr subst
            }
            for (i__ = *n - 2;
                    i__ >= 1;
                    --i__)
            {
                i__1 = i__ + j * b_dim1;
                i__2 = i__ + j * b_dim1;
                i__3 = i__;
                i__4 = i__ + 1 + j * b_dim1;
                z__4.r = du[i__3].r * b[i__4].r - du[i__3].i * b[i__4].i;
                z__4.i = du[i__3].r * b[i__4].i + du[i__3].i * b[i__4] .r; // , expr subst
                z__3.r = b[i__2].r - z__4.r;
                z__3.i = b[i__2].i - z__4.i; // , expr subst
                i__5 = i__;
                i__6 = i__ + 2 + j * b_dim1;
                z__5.r = du2[i__5].r * b[i__6].r - du2[i__5].i * b[i__6].i;
                z__5.i = du2[i__5].r * b[i__6].i + du2[i__5].i * b[ i__6].r; // , expr subst
                z__2.r = z__3.r - z__5.r;
                z__2.i = z__3.i - z__5.i; // , expr subst
                z_div(&z__1, &z__2, &d__[i__]);
                b[i__1].r = z__1.r;
                b[i__1].i = z__1.i; // , expr subst
                /* L30: */
            }
            if (j < *nrhs)
            {
                ++j;
                goto L10;
            }
        }
        else
        {
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Solve L*x = b. */
                i__2 = *n - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    if (ipiv[i__] == i__)
                    {
                        i__3 = i__ + 1 + j * b_dim1;
                        i__4 = i__ + 1 + j * b_dim1;
                        i__5 = i__;
                        i__6 = i__ + j * b_dim1;
                        z__2.r = dl[i__5].r * b[i__6].r - dl[i__5].i * b[i__6] .i;
                        z__2.i = dl[i__5].r * b[i__6].i + dl[i__5] .i * b[i__6].r; // , expr subst
                        z__1.r = b[i__4].r - z__2.r;
                        z__1.i = b[i__4].i - z__2.i; // , expr subst
                        b[i__3].r = z__1.r;
                        b[i__3].i = z__1.i; // , expr subst
                    }
                    else
                    {
                        i__3 = i__ + j * b_dim1;
                        temp.r = b[i__3].r;
                        temp.i = b[i__3].i; // , expr subst
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + 1 + j * b_dim1;
                        b[i__3].r = b[i__4].r;
                        b[i__3].i = b[i__4].i; // , expr subst
                        i__3 = i__ + 1 + j * b_dim1;
                        i__4 = i__;
                        i__5 = i__ + j * b_dim1;
                        z__2.r = dl[i__4].r * b[i__5].r - dl[i__4].i * b[i__5] .i;
                        z__2.i = dl[i__4].r * b[i__5].i + dl[i__4] .i * b[i__5].r; // , expr subst
                        z__1.r = temp.r - z__2.r;
                        z__1.i = temp.i - z__2.i; // , expr subst
                        b[i__3].r = z__1.r;
                        b[i__3].i = z__1.i; // , expr subst
                    }
                    /* L40: */
                }
                /* Solve U*x = b. */
                i__2 = *n + j * b_dim1;
                z_div(&z__1, &b[*n + j * b_dim1], &d__[*n]);
                b[i__2].r = z__1.r;
                b[i__2].i = z__1.i; // , expr subst
                if (*n > 1)
                {
                    i__2 = *n - 1 + j * b_dim1;
                    i__3 = *n - 1 + j * b_dim1;
                    i__4 = *n - 1;
                    i__5 = *n + j * b_dim1;
                    z__3.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i;
                    z__3.i = du[i__4].r * b[i__5].i + du[i__4].i * b[ i__5].r; // , expr subst
                    z__2.r = b[i__3].r - z__3.r;
                    z__2.i = b[i__3].i - z__3.i; // , expr subst
                    z_div(&z__1, &z__2, &d__[*n - 1]);
                    b[i__2].r = z__1.r;
                    b[i__2].i = z__1.i; // , expr subst
                }
                for (i__ = *n - 2;
                        i__ >= 1;
                        --i__)
                {
                    i__2 = i__ + j * b_dim1;
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__;
                    i__5 = i__ + 1 + j * b_dim1;
                    z__4.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i;
                    z__4.i = du[i__4].r * b[i__5].i + du[i__4].i * b[ i__5].r; // , expr subst
                    z__3.r = b[i__3].r - z__4.r;
                    z__3.i = b[i__3].i - z__4.i; // , expr subst
                    i__6 = i__;
                    i__7 = i__ + 2 + j * b_dim1;
                    z__5.r = du2[i__6].r * b[i__7].r - du2[i__6].i * b[i__7] .i;
                    z__5.i = du2[i__6].r * b[i__7].i + du2[i__6] .i * b[i__7].r; // , expr subst
                    z__2.r = z__3.r - z__5.r;
                    z__2.i = z__3.i - z__5.i; // , expr subst
                    z_div(&z__1, &z__2, &d__[i__]);
                    b[i__2].r = z__1.r;
                    b[i__2].i = z__1.i; // , expr subst
                    /* L50: */
                }
                /* L60: */
            }
        }
    }
    else if (*itrans == 1)
    {
        /* Solve A**T * X = B. */
        if (*nrhs <= 1)
        {
            j = 1;
L70: /* Solve U**T * x = b. */
            i__1 = j * b_dim1 + 1;
            z_div(&z__1, &b[j * b_dim1 + 1], &d__[1]);
            b[i__1].r = z__1.r;
            b[i__1].i = z__1.i; // , expr subst
            if (*n > 1)
            {
                i__1 = j * b_dim1 + 2;
                i__2 = j * b_dim1 + 2;
                i__3 = j * b_dim1 + 1;
                z__3.r = du[1].r * b[i__3].r - du[1].i * b[i__3].i;
                z__3.i = du[1].r * b[i__3].i + du[1].i * b[i__3].r; // , expr subst
                z__2.r = b[i__2].r - z__3.r;
                z__2.i = b[i__2].i - z__3.i; // , expr subst
                z_div(&z__1, &z__2, &d__[2]);
                b[i__1].r = z__1.r;
                b[i__1].i = z__1.i; // , expr subst
            }
            i__1 = *n;
            for (i__ = 3;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = i__ + j * b_dim1;
                i__3 = i__ + j * b_dim1;
                i__4 = i__ - 1;
                i__5 = i__ - 1 + j * b_dim1;
                z__4.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i;
                z__4.i = du[i__4].r * b[i__5].i + du[i__4].i * b[i__5] .r; // , expr subst
                z__3.r = b[i__3].r - z__4.r;
                z__3.i = b[i__3].i - z__4.i; // , expr subst
                i__6 = i__ - 2;
                i__7 = i__ - 2 + j * b_dim1;
                z__5.r = du2[i__6].r * b[i__7].r - du2[i__6].i * b[i__7].i;
                z__5.i = du2[i__6].r * b[i__7].i + du2[i__6].i * b[ i__7].r; // , expr subst
                z__2.r = z__3.r - z__5.r;
                z__2.i = z__3.i - z__5.i; // , expr subst
                z_div(&z__1, &z__2, &d__[i__]);
                b[i__2].r = z__1.r;
                b[i__2].i = z__1.i; // , expr subst
                /* L80: */
            }
            /* Solve L**T * x = b. */
            for (i__ = *n - 1;
                    i__ >= 1;
                    --i__)
            {
                if (ipiv[i__] == i__)
                {
                    i__1 = i__ + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    i__3 = i__;
                    i__4 = i__ + 1 + j * b_dim1;
                    z__2.r = dl[i__3].r * b[i__4].r - dl[i__3].i * b[i__4].i;
                    z__2.i = dl[i__3].r * b[i__4].i + dl[i__3].i * b[ i__4].r; // , expr subst
                    z__1.r = b[i__2].r - z__2.r;
                    z__1.i = b[i__2].i - z__2.i; // , expr subst
                    b[i__1].r = z__1.r;
                    b[i__1].i = z__1.i; // , expr subst
                }
                else
                {
                    i__1 = i__ + 1 + j * b_dim1;
                    temp.r = b[i__1].r;
                    temp.i = b[i__1].i; // , expr subst
                    i__1 = i__ + 1 + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    i__3 = i__;
                    z__2.r = dl[i__3].r * temp.r - dl[i__3].i * temp.i;
                    z__2.i = dl[i__3].r * temp.i + dl[i__3].i * temp.r; // , expr subst
                    z__1.r = b[i__2].r - z__2.r;
                    z__1.i = b[i__2].i - z__2.i; // , expr subst
                    b[i__1].r = z__1.r;
                    b[i__1].i = z__1.i; // , expr subst
                    i__1 = i__ + j * b_dim1;
                    b[i__1].r = temp.r;
                    b[i__1].i = temp.i; // , expr subst
                }
                /* L90: */
            }
            if (j < *nrhs)
            {
                ++j;
                goto L70;
            }
        }
        else
        {
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Solve U**T * x = b. */
                i__2 = j * b_dim1 + 1;
                z_div(&z__1, &b[j * b_dim1 + 1], &d__[1]);
                b[i__2].r = z__1.r;
                b[i__2].i = z__1.i; // , expr subst
                if (*n > 1)
                {
                    i__2 = j * b_dim1 + 2;
                    i__3 = j * b_dim1 + 2;
                    i__4 = j * b_dim1 + 1;
                    z__3.r = du[1].r * b[i__4].r - du[1].i * b[i__4].i;
                    z__3.i = du[1].r * b[i__4].i + du[1].i * b[i__4] .r; // , expr subst
                    z__2.r = b[i__3].r - z__3.r;
                    z__2.i = b[i__3].i - z__3.i; // , expr subst
                    z_div(&z__1, &z__2, &d__[2]);
                    b[i__2].r = z__1.r;
                    b[i__2].i = z__1.i; // , expr subst
                }
                i__2 = *n;
                for (i__ = 3;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__ + j * b_dim1;
                    i__5 = i__ - 1;
                    i__6 = i__ - 1 + j * b_dim1;
                    z__4.r = du[i__5].r * b[i__6].r - du[i__5].i * b[i__6].i;
                    z__4.i = du[i__5].r * b[i__6].i + du[i__5].i * b[ i__6].r; // , expr subst
                    z__3.r = b[i__4].r - z__4.r;
                    z__3.i = b[i__4].i - z__4.i; // , expr subst
                    i__7 = i__ - 2;
                    i__8 = i__ - 2 + j * b_dim1;
                    z__5.r = du2[i__7].r * b[i__8].r - du2[i__7].i * b[i__8] .i;
                    z__5.i = du2[i__7].r * b[i__8].i + du2[i__7] .i * b[i__8].r; // , expr subst
                    z__2.r = z__3.r - z__5.r;
                    z__2.i = z__3.i - z__5.i; // , expr subst
                    z_div(&z__1, &z__2, &d__[i__]);
                    b[i__3].r = z__1.r;
                    b[i__3].i = z__1.i; // , expr subst
                    /* L100: */
                }
                /* Solve L**T * x = b. */
                for (i__ = *n - 1;
                        i__ >= 1;
                        --i__)
                {
                    if (ipiv[i__] == i__)
                    {
                        i__2 = i__ + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__;
                        i__5 = i__ + 1 + j * b_dim1;
                        z__2.r = dl[i__4].r * b[i__5].r - dl[i__4].i * b[i__5] .i;
                        z__2.i = dl[i__4].r * b[i__5].i + dl[i__4] .i * b[i__5].r; // , expr subst
                        z__1.r = b[i__3].r - z__2.r;
                        z__1.i = b[i__3].i - z__2.i; // , expr subst
                        b[i__2].r = z__1.r;
                        b[i__2].i = z__1.i; // , expr subst
                    }
                    else
                    {
                        i__2 = i__ + 1 + j * b_dim1;
                        temp.r = b[i__2].r;
                        temp.i = b[i__2].i; // , expr subst
                        i__2 = i__ + 1 + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__;
                        z__2.r = dl[i__4].r * temp.r - dl[i__4].i * temp.i;
                        z__2.i = dl[i__4].r * temp.i + dl[i__4].i * temp.r; // , expr subst
                        z__1.r = b[i__3].r - z__2.r;
                        z__1.i = b[i__3].i - z__2.i; // , expr subst
                        b[i__2].r = z__1.r;
                        b[i__2].i = z__1.i; // , expr subst
                        i__2 = i__ + j * b_dim1;
                        b[i__2].r = temp.r;
                        b[i__2].i = temp.i; // , expr subst
                    }
                    /* L110: */
                }
                /* L120: */
            }
        }
    }
    else
    {
        /* Solve A**H * X = B. */
        if (*nrhs <= 1)
        {
            j = 1;
L130: /* Solve U**H * x = b. */
            i__1 = j * b_dim1 + 1;
            d_cnjg(&z__2, &d__[1]);
            z_div(&z__1, &b[j * b_dim1 + 1], &z__2);
            b[i__1].r = z__1.r;
            b[i__1].i = z__1.i; // , expr subst
            if (*n > 1)
            {
                i__1 = j * b_dim1 + 2;
                i__2 = j * b_dim1 + 2;
                d_cnjg(&z__4, &du[1]);
                i__3 = j * b_dim1 + 1;
                z__3.r = z__4.r * b[i__3].r - z__4.i * b[i__3].i;
                z__3.i = z__4.r * b[i__3].i + z__4.i * b[i__3].r; // , expr subst
                z__2.r = b[i__2].r - z__3.r;
                z__2.i = b[i__2].i - z__3.i; // , expr subst
                d_cnjg(&z__5, &d__[2]);
                z_div(&z__1, &z__2, &z__5);
                b[i__1].r = z__1.r;
                b[i__1].i = z__1.i; // , expr subst
            }
            i__1 = *n;
            for (i__ = 3;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = i__ + j * b_dim1;
                i__3 = i__ + j * b_dim1;
                d_cnjg(&z__5, &du[i__ - 1]);
                i__4 = i__ - 1 + j * b_dim1;
                z__4.r = z__5.r * b[i__4].r - z__5.i * b[i__4].i;
                z__4.i = z__5.r * b[i__4].i + z__5.i * b[i__4].r; // , expr subst
                z__3.r = b[i__3].r - z__4.r;
                z__3.i = b[i__3].i - z__4.i; // , expr subst
                d_cnjg(&z__7, &du2[i__ - 2]);
                i__5 = i__ - 2 + j * b_dim1;
                z__6.r = z__7.r * b[i__5].r - z__7.i * b[i__5].i;
                z__6.i = z__7.r * b[i__5].i + z__7.i * b[i__5].r; // , expr subst
                z__2.r = z__3.r - z__6.r;
                z__2.i = z__3.i - z__6.i; // , expr subst
                d_cnjg(&z__8, &d__[i__]);
                z_div(&z__1, &z__2, &z__8);
                b[i__2].r = z__1.r;
                b[i__2].i = z__1.i; // , expr subst
                /* L140: */
            }
            /* Solve L**H * x = b. */
            for (i__ = *n - 1;
                    i__ >= 1;
                    --i__)
            {
                if (ipiv[i__] == i__)
                {
                    i__1 = i__ + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    d_cnjg(&z__3, &dl[i__]);
                    i__3 = i__ + 1 + j * b_dim1;
                    z__2.r = z__3.r * b[i__3].r - z__3.i * b[i__3].i;
                    z__2.i = z__3.r * b[i__3].i + z__3.i * b[i__3].r; // , expr subst
                    z__1.r = b[i__2].r - z__2.r;
                    z__1.i = b[i__2].i - z__2.i; // , expr subst
                    b[i__1].r = z__1.r;
                    b[i__1].i = z__1.i; // , expr subst
                }
                else
                {
                    i__1 = i__ + 1 + j * b_dim1;
                    temp.r = b[i__1].r;
                    temp.i = b[i__1].i; // , expr subst
                    i__1 = i__ + 1 + j * b_dim1;
                    i__2 = i__ + j * b_dim1;
                    d_cnjg(&z__3, &dl[i__]);
                    z__2.r = z__3.r * temp.r - z__3.i * temp.i;
                    z__2.i = z__3.r * temp.i + z__3.i * temp.r; // , expr subst
                    z__1.r = b[i__2].r - z__2.r;
                    z__1.i = b[i__2].i - z__2.i; // , expr subst
                    b[i__1].r = z__1.r;
                    b[i__1].i = z__1.i; // , expr subst
                    i__1 = i__ + j * b_dim1;
                    b[i__1].r = temp.r;
                    b[i__1].i = temp.i; // , expr subst
                }
                /* L150: */
            }
            if (j < *nrhs)
            {
                ++j;
                goto L130;
            }
        }
        else
        {
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                /* Solve U**H * x = b. */
                i__2 = j * b_dim1 + 1;
                d_cnjg(&z__2, &d__[1]);
                z_div(&z__1, &b[j * b_dim1 + 1], &z__2);
                b[i__2].r = z__1.r;
                b[i__2].i = z__1.i; // , expr subst
                if (*n > 1)
                {
                    i__2 = j * b_dim1 + 2;
                    i__3 = j * b_dim1 + 2;
                    d_cnjg(&z__4, &du[1]);
                    i__4 = j * b_dim1 + 1;
                    z__3.r = z__4.r * b[i__4].r - z__4.i * b[i__4].i;
                    z__3.i = z__4.r * b[i__4].i + z__4.i * b[i__4].r; // , expr subst
                    z__2.r = b[i__3].r - z__3.r;
                    z__2.i = b[i__3].i - z__3.i; // , expr subst
                    d_cnjg(&z__5, &d__[2]);
                    z_div(&z__1, &z__2, &z__5);
                    b[i__2].r = z__1.r;
                    b[i__2].i = z__1.i; // , expr subst
                }
                i__2 = *n;
                for (i__ = 3;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__ + j * b_dim1;
                    d_cnjg(&z__5, &du[i__ - 1]);
                    i__5 = i__ - 1 + j * b_dim1;
                    z__4.r = z__5.r * b[i__5].r - z__5.i * b[i__5].i;
                    z__4.i = z__5.r * b[i__5].i + z__5.i * b[i__5].r; // , expr subst
                    z__3.r = b[i__4].r - z__4.r;
                    z__3.i = b[i__4].i - z__4.i; // , expr subst
                    d_cnjg(&z__7, &du2[i__ - 2]);
                    i__6 = i__ - 2 + j * b_dim1;
                    z__6.r = z__7.r * b[i__6].r - z__7.i * b[i__6].i;
                    z__6.i = z__7.r * b[i__6].i + z__7.i * b[i__6].r; // , expr subst
                    z__2.r = z__3.r - z__6.r;
                    z__2.i = z__3.i - z__6.i; // , expr subst
                    d_cnjg(&z__8, &d__[i__]);
                    z_div(&z__1, &z__2, &z__8);
                    b[i__3].r = z__1.r;
                    b[i__3].i = z__1.i; // , expr subst
                    /* L160: */
                }
                /* Solve L**H * x = b. */
                for (i__ = *n - 1;
                        i__ >= 1;
                        --i__)
                {
                    if (ipiv[i__] == i__)
                    {
                        i__2 = i__ + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        d_cnjg(&z__3, &dl[i__]);
                        i__4 = i__ + 1 + j * b_dim1;
                        z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i;
                        z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4] .r; // , expr subst
                        z__1.r = b[i__3].r - z__2.r;
                        z__1.i = b[i__3].i - z__2.i; // , expr subst
                        b[i__2].r = z__1.r;
                        b[i__2].i = z__1.i; // , expr subst
                    }
                    else
                    {
                        i__2 = i__ + 1 + j * b_dim1;
                        temp.r = b[i__2].r;
                        temp.i = b[i__2].i; // , expr subst
                        i__2 = i__ + 1 + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        d_cnjg(&z__3, &dl[i__]);
                        z__2.r = z__3.r * temp.r - z__3.i * temp.i;
                        z__2.i = z__3.r * temp.i + z__3.i * temp.r; // , expr subst
                        z__1.r = b[i__3].r - z__2.r;
                        z__1.i = b[i__3].i - z__2.i; // , expr subst
                        b[i__2].r = z__1.r;
                        b[i__2].i = z__1.i; // , expr subst
                        i__2 = i__ + j * b_dim1;
                        b[i__2].r = temp.r;
                        b[i__2].i = temp.i; // , expr subst
                    }
                    /* L170: */
                }
                /* L180: */
            }
        }
    }
    /* End of ZGTTS2 */
    return 0;
}
/* zgtts2_ */
