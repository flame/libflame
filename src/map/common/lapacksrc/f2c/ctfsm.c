/* ../netlib/ctfsm.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    1.f,0.f
}
;
/* > \brief \b CTFSM solves a matrix equation (one operand is a triangular matrix in RFP format). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CTFSM + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctfsm.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctfsm.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctfsm.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, */
/* B, LDB ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANSR, DIAG, SIDE, TRANS, UPLO */
/* INTEGER LDB, M, N */
/* COMPLEX ALPHA */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( 0: * ), B( 0: LDB-1, 0: * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for A in RFP Format. */
/* > */
/* > CTFSM solves the matrix equation */
/* > */
/* > op( A )*X = alpha*B or X*op( A ) = alpha*B */
/* > */
/* > where alpha is a scalar, X and B are m by n matrices, A is a unit, or */
/* > non-unit, upper or lower triangular matrix and op( A ) is one of */
/* > */
/* > op( A ) = A or op( A ) = A**H. */
/* > */
/* > A is in Rectangular Full Packed (RFP) Format. */
/* > */
/* > The matrix X is overwritten on B. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANSR */
/* > \verbatim */
/* > TRANSR is CHARACTER*1 */
/* > = 'N': The Normal Form of RFP A is stored;
*/
/* > = 'C': The Conjugate-transpose Form of RFP A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > On entry, SIDE specifies whether op( A ) appears on the left */
/* > or right of X as follows: */
/* > */
/* > SIDE = 'L' or 'l' op( A )*X = alpha*B. */
/* > */
/* > SIDE = 'R' or 'r' X*op( A ) = alpha*B. */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > On entry, UPLO specifies whether the RFP matrix A came from */
/* > an upper or lower triangular matrix as follows: */
/* > UPLO = 'U' or 'u' RFP A came from an upper triangular matrix */
/* > UPLO = 'L' or 'l' RFP A came from a lower triangular matrix */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > On entry, TRANS specifies the form of op( A ) to be used */
/* > in the matrix multiplication as follows: */
/* > */
/* > TRANS = 'N' or 'n' op( A ) = A. */
/* > */
/* > TRANS = 'C' or 'c' op( A ) = conjg( A' ). */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* > DIAG is CHARACTER*1 */
/* > On entry, DIAG specifies whether or not RFP A is unit */
/* > triangular as follows: */
/* > */
/* > DIAG = 'U' or 'u' A is assumed to be unit triangular. */
/* > */
/* > DIAG = 'N' or 'n' A is not assumed to be unit */
/* > triangular. */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > On entry, M specifies the number of rows of B. M must be at */
/* > least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > On entry, N specifies the number of columns of B. N must be */
/* > at least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is COMPLEX */
/* > On entry, ALPHA specifies the scalar alpha. When alpha is */
/* > zero then A is not referenced and B need not be set before */
/* > entry. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (N*(N+1)/2) */
/* > NT = N*(N+1)/2. On entry, the matrix A in RFP Format. */
/* > RFP Format is described by TRANSR, UPLO and N as follows: */
/* > If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;
*/
/* > K=N/2. RFP A is (0:N-1,0:K) when N is odd;
K=N/2. If */
/* > TRANSR = 'C' then RFP is the Conjugate-transpose of RFP A as */
/* > defined when TRANSR = 'N'. The contents of RFP A are defined */
/* > by UPLO as follows: If UPLO = 'U' the RFP A contains the NT */
/* > elements of upper packed A either in normal or */
/* > conjugate-transpose Format. If UPLO = 'L' the RFP A contains */
/* > the NT elements of lower packed A either in normal or */
/* > conjugate-transpose Format. The LDA of RFP A is (N+1)/2 when */
/* > TRANSR = 'C'. When TRANSR is 'N' the LDA is N+1 when N is */
/* > even and is N when is odd. */
/* > See the Note below for more details. Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,N) */
/* > Before entry, the leading m by n part of the array B must */
/* > contain the right-hand side matrix B, and on exit is */
/* > overwritten by the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > On entry, LDB specifies the first dimension of B as declared */
/* > in the calling (sub) program. LDB must be at least */
/* > max( 1, m ). */
/* > Unchanged on exit. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > We first consider Standard Packed Format when N is even. */
/* > We give an example where N = 6. */
/* > */
/* > AP is Upper AP is Lower */
/* > */
/* > 00 01 02 03 04 05 00 */
/* > 11 12 13 14 15 10 11 */
/* > 22 23 24 25 20 21 22 */
/* > 33 34 35 30 31 32 33 */
/* > 44 45 40 41 42 43 44 */
/* > 55 50 51 52 53 54 55 */
/* > */
/* > */
/* > Let TRANSR = 'N'. RFP holds AP as follows: */
/* > For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last */
/* > three columns of AP upper. The lower triangle A(4:6,0:2) consists of */
/* > conjugate-transpose of the first three columns of AP upper. */
/* > For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first */
/* > three columns of AP lower. The upper triangle A(0:2,0:2) consists of */
/* > conjugate-transpose of the last three columns of AP lower. */
/* > To denote conjugate we place -- above the element. This covers the */
/* > case N even and TRANSR = 'N'. */
/* > */
/* > RFP A RFP A */
/* > */
/* > -- -- -- */
/* > 03 04 05 33 43 53 */
/* > -- -- */
/* > 13 14 15 00 44 54 */
/* > -- */
/* > 23 24 25 10 11 55 */
/* > */
/* > 33 34 35 20 21 22 */
/* > -- */
/* > 00 44 45 30 31 32 */
/* > -- -- */
/* > 01 11 55 40 41 42 */
/* > -- -- -- */
/* > 02 12 22 50 51 52 */
/* > */
/* > Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate- */
/* > transpose of RFP A above. One therefore gets: */
/* > */
/* > */
/* > RFP A RFP A */
/* > */
/* > -- -- -- -- -- -- -- -- -- -- */
/* > 03 13 23 33 00 01 02 33 00 10 20 30 40 50 */
/* > -- -- -- -- -- -- -- -- -- -- */
/* > 04 14 24 34 44 11 12 43 44 11 21 31 41 51 */
/* > -- -- -- -- -- -- -- -- -- -- */
/* > 05 15 25 35 45 55 22 53 54 55 22 32 42 52 */
/* > */
/* > */
/* > We next consider Standard Packed Format when N is odd. */
/* > We give an example where N = 5. */
/* > */
/* > AP is Upper AP is Lower */
/* > */
/* > 00 01 02 03 04 00 */
/* > 11 12 13 14 10 11 */
/* > 22 23 24 20 21 22 */
/* > 33 34 30 31 32 33 */
/* > 44 40 41 42 43 44 */
/* > */
/* > */
/* > Let TRANSR = 'N'. RFP holds AP as follows: */
/* > For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last */
/* > three columns of AP upper. The lower triangle A(3:4,0:1) consists of */
/* > conjugate-transpose of the first two columns of AP upper. */
/* > For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first */
/* > three columns of AP lower. The upper triangle A(0:1,1:2) consists of */
/* > conjugate-transpose of the last two columns of AP lower. */
/* > To denote conjugate we place -- above the element. This covers the */
/* > case N odd and TRANSR = 'N'. */
/* > */
/* > RFP A RFP A */
/* > */
/* > -- -- */
/* > 02 03 04 00 33 43 */
/* > -- */
/* > 12 13 14 10 11 44 */
/* > */
/* > 22 23 24 20 21 22 */
/* > -- */
/* > 00 33 34 30 31 32 */
/* > -- -- */
/* > 01 11 44 40 41 42 */
/* > */
/* > Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate- */
/* > transpose of RFP A above. One therefore gets: */
/* > */
/* > */
/* > RFP A RFP A */
/* > */
/* > -- -- -- -- -- -- -- -- -- */
/* > 02 12 22 00 01 00 10 20 30 40 50 */
/* > -- -- -- -- -- -- -- -- -- */
/* > 03 13 23 33 11 33 11 21 31 41 51 */
/* > -- -- -- -- -- -- -- -- -- */
/* > 04 14 24 34 44 43 44 22 32 42 52 */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int ctfsm_(char *transr, char *side, char *uplo, char *trans, char *diag, integer *m, integer *n, complex *alpha, complex *a, complex *b, integer *ldb)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3;
    complex q__1;
    /* Local variables */
    integer i__, j, k, m1, m2, n1, n2, info;
    logical normaltransr;
    extern /* Subroutine */
    int cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
    logical lside;
    extern logical lsame_(char *, char *);
    logical lower;
    extern /* Subroutine */
    int ctrsm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *), xerbla_(char *, integer *);
    logical misodd, nisodd, notrans;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
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
    /* Test the input parameters. */
    /* Parameter adjustments */
    b_dim1 = *ldb - 1 - 0 + 1;
    b_offset = 0 + b_dim1 * 0;
    b -= b_offset;
    /* Function Body */
    info = 0;
    normaltransr = lsame_(transr, "N");
    lside = lsame_(side, "L");
    lower = lsame_(uplo, "L");
    notrans = lsame_(trans, "N");
    if (! normaltransr && ! lsame_(transr, "C"))
    {
        info = -1;
    }
    else if (! lside && ! lsame_(side, "R"))
    {
        info = -2;
    }
    else if (! lower && ! lsame_(uplo, "U"))
    {
        info = -3;
    }
    else if (! notrans && ! lsame_(trans, "C"))
    {
        info = -4;
    }
    else if (! lsame_(diag, "N") && ! lsame_(diag, "U"))
    {
        info = -5;
    }
    else if (*m < 0)
    {
        info = -6;
    }
    else if (*n < 0)
    {
        info = -7;
    }
    else if (*ldb < max(1,*m))
    {
        info = -11;
    }
    if (info != 0)
    {
        i__1 = -info;
        xerbla_("CTFSM ", &i__1);
        return 0;
    }
    /* Quick return when ( (N.EQ.0).OR.(M.EQ.0) ) */
    if (*m == 0 || *n == 0)
    {
        return 0;
    }
    /* Quick return when ALPHA.EQ.(0E+0,0E+0) */
    if (alpha->r == 0.f && alpha->i == 0.f)
    {
        i__1 = *n - 1;
        for (j = 0;
                j <= i__1;
                ++j)
        {
            i__2 = *m - 1;
            for (i__ = 0;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * b_dim1;
                b[i__3].r = 0.f;
                b[i__3].i = 0.f; // , expr subst
                /* L10: */
            }
            /* L20: */
        }
        return 0;
    }
    if (lside)
    {
        /* SIDE = 'L' */
        /* A is M-by-M. */
        /* If M is odd, set NISODD = .TRUE., and M1 and M2. */
        /* If M is even, NISODD = .FALSE., and M. */
        if (*m % 2 == 0)
        {
            misodd = FALSE_;
            k = *m / 2;
        }
        else
        {
            misodd = TRUE_;
            if (lower)
            {
                m2 = *m / 2;
                m1 = *m - m2;
            }
            else
            {
                m1 = *m / 2;
                m2 = *m - m1;
            }
        }
        if (misodd)
        {
            /* SIDE = 'L' and N is odd */
            if (normaltransr)
            {
                /* SIDE = 'L', N is odd, and TRANSR = 'N' */
                if (lower)
                {
                    /* SIDE ='L', N is odd, TRANSR = 'N', and UPLO = 'L' */
                    if (notrans)
                    {
                        /* SIDE ='L', N is odd, TRANSR = 'N', UPLO = 'L', and */
                        /* TRANS = 'N' */
                        if (*m == 1)
                        {
                            ctrsm_("L", "L", "N", diag, &m1, n, alpha, a, m, & b[b_offset], ldb);
                        }
                        else
                        {
                            ctrsm_("L", "L", "N", diag, &m1, n, alpha, a, m, & b[b_offset], ldb);
                            q__1.r = -1.f;
                            q__1.i = -0.f; // , expr subst
                            cgemm_("N", "N", &m2, n, &m1, &q__1, &a[m1], m, & b[b_offset], ldb, alpha, &b[m1], ldb);
                            ctrsm_("L", "U", "C", diag, &m2, n, &c_b1, &a[*m], m, &b[m1], ldb);
                        }
                    }
                    else
                    {
                        /* SIDE ='L', N is odd, TRANSR = 'N', UPLO = 'L', and */
                        /* TRANS = 'C' */
                        if (*m == 1)
                        {
                            ctrsm_("L", "L", "C", diag, &m1, n, alpha, a, m, & b[b_offset], ldb);
                        }
                        else
                        {
                            ctrsm_("L", "U", "N", diag, &m2, n, alpha, &a[*m], m, &b[m1], ldb);
                            q__1.r = -1.f;
                            q__1.i = -0.f; // , expr subst
                            cgemm_("C", "N", &m1, n, &m2, &q__1, &a[m1], m, & b[m1], ldb, alpha, &b[b_offset], ldb);
                            ctrsm_("L", "L", "C", diag, &m1, n, &c_b1, a, m, & b[b_offset], ldb);
                        }
                    }
                }
                else
                {
                    /* SIDE ='L', N is odd, TRANSR = 'N', and UPLO = 'U' */
                    if (! notrans)
                    {
                        /* SIDE ='L', N is odd, TRANSR = 'N', UPLO = 'U', and */
                        /* TRANS = 'N' */
                        ctrsm_("L", "L", "N", diag, &m1, n, alpha, &a[m2], m, &b[b_offset], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("C", "N", &m2, n, &m1, &q__1, a, m, &b[ b_offset], ldb, alpha, &b[m1], ldb);
                        ctrsm_("L", "U", "C", diag, &m2, n, &c_b1, &a[m1], m, &b[m1], ldb);
                    }
                    else
                    {
                        /* SIDE ='L', N is odd, TRANSR = 'N', UPLO = 'U', and */
                        /* TRANS = 'C' */
                        ctrsm_("L", "U", "N", diag, &m2, n, alpha, &a[m1], m, &b[m1], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "N", &m1, n, &m2, &q__1, a, m, &b[m1], ldb, alpha, &b[b_offset], ldb);
                        ctrsm_("L", "L", "C", diag, &m1, n, &c_b1, &a[m2], m, &b[b_offset], ldb);
                    }
                }
            }
            else
            {
                /* SIDE = 'L', N is odd, and TRANSR = 'C' */
                if (lower)
                {
                    /* SIDE ='L', N is odd, TRANSR = 'C', and UPLO = 'L' */
                    if (notrans)
                    {
                        /* SIDE ='L', N is odd, TRANSR = 'C', UPLO = 'L', and */
                        /* TRANS = 'N' */
                        if (*m == 1)
                        {
                            ctrsm_("L", "U", "C", diag, &m1, n, alpha, a, &m1, &b[b_offset], ldb);
                        }
                        else
                        {
                            ctrsm_("L", "U", "C", diag, &m1, n, alpha, a, &m1, &b[b_offset], ldb);
                            q__1.r = -1.f;
                            q__1.i = -0.f; // , expr subst
                            cgemm_("C", "N", &m2, n, &m1, &q__1, &a[m1 * m1], &m1, &b[b_offset], ldb, alpha, &b[m1], ldb);
                            ctrsm_("L", "L", "N", diag, &m2, n, &c_b1, &a[1], &m1, &b[m1], ldb);
                        }
                    }
                    else
                    {
                        /* SIDE ='L', N is odd, TRANSR = 'C', UPLO = 'L', and */
                        /* TRANS = 'C' */
                        if (*m == 1)
                        {
                            ctrsm_("L", "U", "N", diag, &m1, n, alpha, a, &m1, &b[b_offset], ldb);
                        }
                        else
                        {
                            ctrsm_("L", "L", "C", diag, &m2, n, alpha, &a[1], &m1, &b[m1], ldb);
                            q__1.r = -1.f;
                            q__1.i = -0.f; // , expr subst
                            cgemm_("N", "N", &m1, n, &m2, &q__1, &a[m1 * m1], &m1, &b[m1], ldb, alpha, &b[b_offset], ldb);
                            ctrsm_("L", "U", "N", diag, &m1, n, &c_b1, a, &m1, &b[b_offset], ldb);
                        }
                    }
                }
                else
                {
                    /* SIDE ='L', N is odd, TRANSR = 'C', and UPLO = 'U' */
                    if (! notrans)
                    {
                        /* SIDE ='L', N is odd, TRANSR = 'C', UPLO = 'U', and */
                        /* TRANS = 'N' */
                        ctrsm_("L", "U", "C", diag, &m1, n, alpha, &a[m2 * m2] , &m2, &b[b_offset], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "N", &m2, n, &m1, &q__1, a, &m2, &b[ b_offset], ldb, alpha, &b[m1], ldb);
                        ctrsm_("L", "L", "N", diag, &m2, n, &c_b1, &a[m1 * m2] , &m2, &b[m1], ldb);
                    }
                    else
                    {
                        /* SIDE ='L', N is odd, TRANSR = 'C', UPLO = 'U', and */
                        /* TRANS = 'C' */
                        ctrsm_("L", "L", "C", diag, &m2, n, alpha, &a[m1 * m2] , &m2, &b[m1], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("C", "N", &m1, n, &m2, &q__1, a, &m2, &b[m1], ldb, alpha, &b[b_offset], ldb);
                        ctrsm_("L", "U", "N", diag, &m1, n, &c_b1, &a[m2 * m2] , &m2, &b[b_offset], ldb);
                    }
                }
            }
        }
        else
        {
            /* SIDE = 'L' and N is even */
            if (normaltransr)
            {
                /* SIDE = 'L', N is even, and TRANSR = 'N' */
                if (lower)
                {
                    /* SIDE ='L', N is even, TRANSR = 'N', and UPLO = 'L' */
                    if (notrans)
                    {
                        /* SIDE ='L', N is even, TRANSR = 'N', UPLO = 'L', */
                        /* and TRANS = 'N' */
                        i__1 = *m + 1;
                        ctrsm_("L", "L", "N", diag, &k, n, alpha, &a[1], & i__1, &b[b_offset], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        i__1 = *m + 1;
                        cgemm_("N", "N", &k, n, &k, &q__1, &a[k + 1], &i__1, & b[b_offset], ldb, alpha, &b[k], ldb);
                        i__1 = *m + 1;
                        ctrsm_("L", "U", "C", diag, &k, n, &c_b1, a, &i__1, & b[k], ldb);
                    }
                    else
                    {
                        /* SIDE ='L', N is even, TRANSR = 'N', UPLO = 'L', */
                        /* and TRANS = 'C' */
                        i__1 = *m + 1;
                        ctrsm_("L", "U", "N", diag, &k, n, alpha, a, &i__1, & b[k], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        i__1 = *m + 1;
                        cgemm_("C", "N", &k, n, &k, &q__1, &a[k + 1], &i__1, & b[k], ldb, alpha, &b[b_offset], ldb);
                        i__1 = *m + 1;
                        ctrsm_("L", "L", "C", diag, &k, n, &c_b1, &a[1], & i__1, &b[b_offset], ldb);
                    }
                }
                else
                {
                    /* SIDE ='L', N is even, TRANSR = 'N', and UPLO = 'U' */
                    if (! notrans)
                    {
                        /* SIDE ='L', N is even, TRANSR = 'N', UPLO = 'U', */
                        /* and TRANS = 'N' */
                        i__1 = *m + 1;
                        ctrsm_("L", "L", "N", diag, &k, n, alpha, &a[k + 1], & i__1, &b[b_offset], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        i__1 = *m + 1;
                        cgemm_("C", "N", &k, n, &k, &q__1, a, &i__1, &b[ b_offset], ldb, alpha, &b[k], ldb);
                        i__1 = *m + 1;
                        ctrsm_("L", "U", "C", diag, &k, n, &c_b1, &a[k], & i__1, &b[k], ldb);
                    }
                    else
                    {
                        /* SIDE ='L', N is even, TRANSR = 'N', UPLO = 'U', */
                        /* and TRANS = 'C' */
                        i__1 = *m + 1;
                        ctrsm_("L", "U", "N", diag, &k, n, alpha, &a[k], & i__1, &b[k], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        i__1 = *m + 1;
                        cgemm_("N", "N", &k, n, &k, &q__1, a, &i__1, &b[k], ldb, alpha, &b[b_offset], ldb);
                        i__1 = *m + 1;
                        ctrsm_("L", "L", "C", diag, &k, n, &c_b1, &a[k + 1], & i__1, &b[b_offset], ldb);
                    }
                }
            }
            else
            {
                /* SIDE = 'L', N is even, and TRANSR = 'C' */
                if (lower)
                {
                    /* SIDE ='L', N is even, TRANSR = 'C', and UPLO = 'L' */
                    if (notrans)
                    {
                        /* SIDE ='L', N is even, TRANSR = 'C', UPLO = 'L', */
                        /* and TRANS = 'N' */
                        ctrsm_("L", "U", "C", diag, &k, n, alpha, &a[k], &k, & b[b_offset], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("C", "N", &k, n, &k, &q__1, &a[k * (k + 1)], & k, &b[b_offset], ldb, alpha, &b[k], ldb);
                        ctrsm_("L", "L", "N", diag, &k, n, &c_b1, a, &k, &b[k] , ldb);
                    }
                    else
                    {
                        /* SIDE ='L', N is even, TRANSR = 'C', UPLO = 'L', */
                        /* and TRANS = 'C' */
                        ctrsm_("L", "L", "C", diag, &k, n, alpha, a, &k, &b[k] , ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "N", &k, n, &k, &q__1, &a[k * (k + 1)], & k, &b[k], ldb, alpha, &b[b_offset], ldb);
                        ctrsm_("L", "U", "N", diag, &k, n, &c_b1, &a[k], &k, & b[b_offset], ldb);
                    }
                }
                else
                {
                    /* SIDE ='L', N is even, TRANSR = 'C', and UPLO = 'U' */
                    if (! notrans)
                    {
                        /* SIDE ='L', N is even, TRANSR = 'C', UPLO = 'U', */
                        /* and TRANS = 'N' */
                        ctrsm_("L", "U", "C", diag, &k, n, alpha, &a[k * (k + 1)], &k, &b[b_offset], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "N", &k, n, &k, &q__1, a, &k, &b[b_offset] , ldb, alpha, &b[k], ldb);
                        ctrsm_("L", "L", "N", diag, &k, n, &c_b1, &a[k * k], & k, &b[k], ldb);
                    }
                    else
                    {
                        /* SIDE ='L', N is even, TRANSR = 'C', UPLO = 'U', */
                        /* and TRANS = 'C' */
                        ctrsm_("L", "L", "C", diag, &k, n, alpha, &a[k * k], & k, &b[k], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("C", "N", &k, n, &k, &q__1, a, &k, &b[k], ldb, alpha, &b[b_offset], ldb);
                        ctrsm_("L", "U", "N", diag, &k, n, &c_b1, &a[k * (k + 1)], &k, &b[b_offset], ldb);
                    }
                }
            }
        }
    }
    else
    {
        /* SIDE = 'R' */
        /* A is N-by-N. */
        /* If N is odd, set NISODD = .TRUE., and N1 and N2. */
        /* If N is even, NISODD = .FALSE., and K. */
        if (*n % 2 == 0)
        {
            nisodd = FALSE_;
            k = *n / 2;
        }
        else
        {
            nisodd = TRUE_;
            if (lower)
            {
                n2 = *n / 2;
                n1 = *n - n2;
            }
            else
            {
                n1 = *n / 2;
                n2 = *n - n1;
            }
        }
        if (nisodd)
        {
            /* SIDE = 'R' and N is odd */
            if (normaltransr)
            {
                /* SIDE = 'R', N is odd, and TRANSR = 'N' */
                if (lower)
                {
                    /* SIDE ='R', N is odd, TRANSR = 'N', and UPLO = 'L' */
                    if (notrans)
                    {
                        /* SIDE ='R', N is odd, TRANSR = 'N', UPLO = 'L', and */
                        /* TRANS = 'N' */
                        ctrsm_("R", "U", "C", diag, m, &n2, alpha, &a[*n], n, &b[n1 * b_dim1], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "N", m, &n1, &n2, &q__1, &b[n1 * b_dim1], ldb, &a[n1], n, alpha, b, ldb);
                        ctrsm_("R", "L", "N", diag, m, &n1, &c_b1, a, n, b, ldb);
                    }
                    else
                    {
                        /* SIDE ='R', N is odd, TRANSR = 'N', UPLO = 'L', and */
                        /* TRANS = 'C' */
                        ctrsm_("R", "L", "C", diag, m, &n1, alpha, a, n, b, ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "C", m, &n2, &n1, &q__1, b, ldb, &a[n1], n, alpha, &b[n1 * b_dim1], ldb);
                        ctrsm_("R", "U", "N", diag, m, &n2, &c_b1, &a[*n], n, &b[n1 * b_dim1], ldb);
                    }
                }
                else
                {
                    /* SIDE ='R', N is odd, TRANSR = 'N', and UPLO = 'U' */
                    if (notrans)
                    {
                        /* SIDE ='R', N is odd, TRANSR = 'N', UPLO = 'U', and */
                        /* TRANS = 'N' */
                        ctrsm_("R", "L", "C", diag, m, &n1, alpha, &a[n2], n, b, ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "N", m, &n2, &n1, &q__1, b, ldb, a, n, alpha, &b[n1 * b_dim1], ldb);
                        ctrsm_("R", "U", "N", diag, m, &n2, &c_b1, &a[n1], n, &b[n1 * b_dim1], ldb);
                    }
                    else
                    {
                        /* SIDE ='R', N is odd, TRANSR = 'N', UPLO = 'U', and */
                        /* TRANS = 'C' */
                        ctrsm_("R", "U", "C", diag, m, &n2, alpha, &a[n1], n, &b[n1 * b_dim1], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "C", m, &n1, &n2, &q__1, &b[n1 * b_dim1], ldb, a, n, alpha, b, ldb);
                        ctrsm_("R", "L", "N", diag, m, &n1, &c_b1, &a[n2], n, b, ldb);
                    }
                }
            }
            else
            {
                /* SIDE = 'R', N is odd, and TRANSR = 'C' */
                if (lower)
                {
                    /* SIDE ='R', N is odd, TRANSR = 'C', and UPLO = 'L' */
                    if (notrans)
                    {
                        /* SIDE ='R', N is odd, TRANSR = 'C', UPLO = 'L', and */
                        /* TRANS = 'N' */
                        ctrsm_("R", "L", "N", diag, m, &n2, alpha, &a[1], &n1, &b[n1 * b_dim1], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "C", m, &n1, &n2, &q__1, &b[n1 * b_dim1], ldb, &a[n1 * n1], &n1, alpha, b, ldb);
                        ctrsm_("R", "U", "C", diag, m, &n1, &c_b1, a, &n1, b, ldb);
                    }
                    else
                    {
                        /* SIDE ='R', N is odd, TRANSR = 'C', UPLO = 'L', and */
                        /* TRANS = 'C' */
                        ctrsm_("R", "U", "N", diag, m, &n1, alpha, a, &n1, b, ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "N", m, &n2, &n1, &q__1, b, ldb, &a[n1 * n1], &n1, alpha, &b[n1 * b_dim1], ldb);
                        ctrsm_("R", "L", "C", diag, m, &n2, &c_b1, &a[1], &n1, &b[n1 * b_dim1], ldb);
                    }
                }
                else
                {
                    /* SIDE ='R', N is odd, TRANSR = 'C', and UPLO = 'U' */
                    if (notrans)
                    {
                        /* SIDE ='R', N is odd, TRANSR = 'C', UPLO = 'U', and */
                        /* TRANS = 'N' */
                        ctrsm_("R", "U", "N", diag, m, &n1, alpha, &a[n2 * n2] , &n2, b, ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "C", m, &n2, &n1, &q__1, b, ldb, a, &n2, alpha, &b[n1 * b_dim1], ldb);
                        ctrsm_("R", "L", "C", diag, m, &n2, &c_b1, &a[n1 * n2] , &n2, &b[n1 * b_dim1], ldb);
                    }
                    else
                    {
                        /* SIDE ='R', N is odd, TRANSR = 'C', UPLO = 'U', and */
                        /* TRANS = 'C' */
                        ctrsm_("R", "L", "N", diag, m, &n2, alpha, &a[n1 * n2] , &n2, &b[n1 * b_dim1], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "N", m, &n1, &n2, &q__1, &b[n1 * b_dim1], ldb, a, &n2, alpha, b, ldb);
                        ctrsm_("R", "U", "C", diag, m, &n1, &c_b1, &a[n2 * n2] , &n2, b, ldb);
                    }
                }
            }
        }
        else
        {
            /* SIDE = 'R' and N is even */
            if (normaltransr)
            {
                /* SIDE = 'R', N is even, and TRANSR = 'N' */
                if (lower)
                {
                    /* SIDE ='R', N is even, TRANSR = 'N', and UPLO = 'L' */
                    if (notrans)
                    {
                        /* SIDE ='R', N is even, TRANSR = 'N', UPLO = 'L', */
                        /* and TRANS = 'N' */
                        i__1 = *n + 1;
                        ctrsm_("R", "U", "C", diag, m, &k, alpha, a, &i__1, & b[k * b_dim1], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        i__1 = *n + 1;
                        cgemm_("N", "N", m, &k, &k, &q__1, &b[k * b_dim1], ldb, &a[k + 1], &i__1, alpha, b, ldb);
                        i__1 = *n + 1;
                        ctrsm_("R", "L", "N", diag, m, &k, &c_b1, &a[1], & i__1, b, ldb);
                    }
                    else
                    {
                        /* SIDE ='R', N is even, TRANSR = 'N', UPLO = 'L', */
                        /* and TRANS = 'C' */
                        i__1 = *n + 1;
                        ctrsm_("R", "L", "C", diag, m, &k, alpha, &a[1], & i__1, b, ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        i__1 = *n + 1;
                        cgemm_("N", "C", m, &k, &k, &q__1, b, ldb, &a[k + 1], &i__1, alpha, &b[k * b_dim1], ldb);
                        i__1 = *n + 1;
                        ctrsm_("R", "U", "N", diag, m, &k, &c_b1, a, &i__1, & b[k * b_dim1], ldb);
                    }
                }
                else
                {
                    /* SIDE ='R', N is even, TRANSR = 'N', and UPLO = 'U' */
                    if (notrans)
                    {
                        /* SIDE ='R', N is even, TRANSR = 'N', UPLO = 'U', */
                        /* and TRANS = 'N' */
                        i__1 = *n + 1;
                        ctrsm_("R", "L", "C", diag, m, &k, alpha, &a[k + 1], & i__1, b, ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        i__1 = *n + 1;
                        cgemm_("N", "N", m, &k, &k, &q__1, b, ldb, a, &i__1, alpha, &b[k * b_dim1], ldb);
                        i__1 = *n + 1;
                        ctrsm_("R", "U", "N", diag, m, &k, &c_b1, &a[k], & i__1, &b[k * b_dim1], ldb);
                    }
                    else
                    {
                        /* SIDE ='R', N is even, TRANSR = 'N', UPLO = 'U', */
                        /* and TRANS = 'C' */
                        i__1 = *n + 1;
                        ctrsm_("R", "U", "C", diag, m, &k, alpha, &a[k], & i__1, &b[k * b_dim1], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        i__1 = *n + 1;
                        cgemm_("N", "C", m, &k, &k, &q__1, &b[k * b_dim1], ldb, a, &i__1, alpha, b, ldb);
                        i__1 = *n + 1;
                        ctrsm_("R", "L", "N", diag, m, &k, &c_b1, &a[k + 1], & i__1, b, ldb);
                    }
                }
            }
            else
            {
                /* SIDE = 'R', N is even, and TRANSR = 'C' */
                if (lower)
                {
                    /* SIDE ='R', N is even, TRANSR = 'C', and UPLO = 'L' */
                    if (notrans)
                    {
                        /* SIDE ='R', N is even, TRANSR = 'C', UPLO = 'L', */
                        /* and TRANS = 'N' */
                        ctrsm_("R", "L", "N", diag, m, &k, alpha, a, &k, &b[k * b_dim1], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "C", m, &k, &k, &q__1, &b[k * b_dim1], ldb, &a[(k + 1) * k], &k, alpha, b, ldb);
                        ctrsm_("R", "U", "C", diag, m, &k, &c_b1, &a[k], &k, b, ldb);
                    }
                    else
                    {
                        /* SIDE ='R', N is even, TRANSR = 'C', UPLO = 'L', */
                        /* and TRANS = 'C' */
                        ctrsm_("R", "U", "N", diag, m, &k, alpha, &a[k], &k, b, ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "N", m, &k, &k, &q__1, b, ldb, &a[(k + 1) * k], &k, alpha, &b[k * b_dim1], ldb);
                        ctrsm_("R", "L", "C", diag, m, &k, &c_b1, a, &k, &b[k * b_dim1], ldb);
                    }
                }
                else
                {
                    /* SIDE ='R', N is even, TRANSR = 'C', and UPLO = 'U' */
                    if (notrans)
                    {
                        /* SIDE ='R', N is even, TRANSR = 'C', UPLO = 'U', */
                        /* and TRANS = 'N' */
                        ctrsm_("R", "U", "N", diag, m, &k, alpha, &a[(k + 1) * k], &k, b, ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "C", m, &k, &k, &q__1, b, ldb, a, &k, alpha, &b[k * b_dim1], ldb);
                        ctrsm_("R", "L", "C", diag, m, &k, &c_b1, &a[k * k], & k, &b[k * b_dim1], ldb);
                    }
                    else
                    {
                        /* SIDE ='R', N is even, TRANSR = 'C', UPLO = 'U', */
                        /* and TRANS = 'C' */
                        ctrsm_("R", "L", "N", diag, m, &k, alpha, &a[k * k], & k, &b[k * b_dim1], ldb);
                        q__1.r = -1.f;
                        q__1.i = -0.f; // , expr subst
                        cgemm_("N", "N", m, &k, &k, &q__1, &b[k * b_dim1], ldb, a, &k, alpha, b, ldb);
                        ctrsm_("R", "U", "C", diag, m, &k, &c_b1, &a[(k + 1) * k], &k, b, ldb);
                    }
                }
            }
        }
    }
    return 0;
    /* End of CTFSM */
}
/* ctfsm_ */
