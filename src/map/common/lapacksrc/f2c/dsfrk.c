/* ../netlib/dsfrk.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DSFRK performs a symmetric rank-k operation for matrix in RFP format. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSFRK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsfrk.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsfrk.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsfrk.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, */
/* C ) */
/* .. Scalar Arguments .. */
/* DOUBLE PRECISION ALPHA, BETA */
/* INTEGER K, LDA, N */
/* CHARACTER TRANS, TRANSR, UPLO */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), C( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > Level 3 BLAS like routine for C in RFP Format. */
/* > */
/* > DSFRK performs one of the symmetric rank--k operations */
/* > */
/* > C := alpha*A*A**T + beta*C, */
/* > */
/* > or */
/* > */
/* > C := alpha*A**T*A + beta*C, */
/* > */
/* > where alpha and beta are real scalars, C is an n--by--n symmetric */
/* > matrix and A is an n--by--k matrix in the first case and a k--by--n */
/* > matrix in the second case. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANSR */
/* > \verbatim */
/* > TRANSR is CHARACTER*1 */
/* > = 'N': The Normal Form of RFP A is stored;
*/
/* > = 'T': The Transpose Form of RFP A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > On entry, UPLO specifies whether the upper or lower */
/* > triangular part of the array C is to be referenced as */
/* > follows: */
/* > */
/* > UPLO = 'U' or 'u' Only the upper triangular part of C */
/* > is to be referenced. */
/* > */
/* > UPLO = 'L' or 'l' Only the lower triangular part of C */
/* > is to be referenced. */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > On entry, TRANS specifies the operation to be performed as */
/* > follows: */
/* > */
/* > TRANS = 'N' or 'n' C := alpha*A*A**T + beta*C. */
/* > */
/* > TRANS = 'T' or 't' C := alpha*A**T*A + beta*C. */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > On entry, N specifies the order of the matrix C. N must be */
/* > at least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > On entry with TRANS = 'N' or 'n', K specifies the number */
/* > of columns of the matrix A, and on entry with TRANS = 'T' */
/* > or 't', K specifies the number of rows of the matrix A. K */
/* > must be at least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is DOUBLE PRECISION */
/* > On entry, ALPHA specifies the scalar alpha. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,ka) */
/* > where KA */
/* > is K when TRANS = 'N' or 'n', and is N otherwise. Before */
/* > entry with TRANS = 'N' or 'n', the leading N--by--K part of */
/* > the array A must contain the matrix A, otherwise the leading */
/* > K--by--N part of the array A must contain the matrix A. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > On entry, LDA specifies the first dimension of A as declared */
/* > in the calling (sub) program. When TRANS = 'N' or 'n' */
/* > then LDA must be at least max( 1, n ), otherwise LDA must */
/* > be at least max( 1, k ). */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* > BETA is DOUBLE PRECISION */
/* > On entry, BETA specifies the scalar beta. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (NT) */
/* > NT = N*(N+1)/2. On entry, the symmetric matrix C in RFP */
/* > Format. RFP Format is described by TRANSR, UPLO and N. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int dsfrk_(char *transr, char *uplo, char *trans, integer *n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *beta, doublereal *c__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    integer j, n1, n2, nk, info;
    logical normaltransr;
    extern /* Subroutine */
    int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *);
    integer nrowa;
    logical lower;
    extern /* Subroutine */
    int dsyrk_(char *, char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
    logical nisodd, notrans;
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --c__;
    /* Function Body */
    info = 0;
    normaltransr = lsame_(transr, "N");
    lower = lsame_(uplo, "L");
    notrans = lsame_(trans, "N");
    if (notrans)
    {
        nrowa = *n;
    }
    else
    {
        nrowa = *k;
    }
    if (! normaltransr && ! lsame_(transr, "T"))
    {
        info = -1;
    }
    else if (! lower && ! lsame_(uplo, "U"))
    {
        info = -2;
    }
    else if (! notrans && ! lsame_(trans, "T"))
    {
        info = -3;
    }
    else if (*n < 0)
    {
        info = -4;
    }
    else if (*k < 0)
    {
        info = -5;
    }
    else if (*lda < max(1,nrowa))
    {
        info = -8;
    }
    if (info != 0)
    {
        i__1 = -info;
        xerbla_("DSFRK ", &i__1);
        return 0;
    }
    /* Quick return if possible. */
    /* The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not */
    /* done (it is in DSYRK for example) and left in the general case. */
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.)
    {
        return 0;
    }
    if (*alpha == 0. && *beta == 0.)
    {
        i__1 = *n * (*n + 1) / 2;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            c__[j] = 0.;
        }
        return 0;
    }
    /* C is N-by-N. */
    /* If N is odd, set NISODD = .TRUE., and N1 and N2. */
    /* If N is even, NISODD = .FALSE., and NK. */
    if (*n % 2 == 0)
    {
        nisodd = FALSE_;
        nk = *n / 2;
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
        /* N is odd */
        if (normaltransr)
        {
            /* N is odd and TRANSR = 'N' */
            if (lower)
            {
                /* N is odd, TRANSR = 'N', and UPLO = 'L' */
                if (notrans)
                {
                    /* N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */
                    dsyrk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[1], n);
                    dsyrk_("U", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, beta, &c__[*n + 1], n);
                    dgemm_("N", "T", &n2, &n1, k, alpha, &a[n1 + 1 + a_dim1], lda, &a[a_dim1 + 1], lda, beta, &c__[n1 + 1], n);
                }
                else
                {
                    /* N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'T' */
                    dsyrk_("L", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[1], n);
                    dsyrk_("U", "T", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1], lda, beta, &c__[*n + 1], n) ;
                    dgemm_("T", "N", &n2, &n1, k, alpha, &a[(n1 + 1) * a_dim1 + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[n1 + 1] , n);
                }
            }
            else
            {
                /* N is odd, TRANSR = 'N', and UPLO = 'U' */
                if (notrans)
                {
                    /* N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */
                    dsyrk_("L", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[n2 + 1], n);
                    dsyrk_("U", "N", &n2, k, alpha, &a[n2 + a_dim1], lda, beta, &c__[n1 + 1], n);
                    dgemm_("N", "T", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, &a[n2 + a_dim1], lda, beta, &c__[1], n);
                }
                else
                {
                    /* N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'T' */
                    dsyrk_("L", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[n2 + 1], n);
                    dsyrk_("U", "T", &n2, k, alpha, &a[n2 * a_dim1 + 1], lda, beta, &c__[n1 + 1], n);
                    dgemm_("T", "N", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, &a[n2 * a_dim1 + 1], lda, beta, &c__[1], n);
                }
            }
        }
        else
        {
            /* N is odd, and TRANSR = 'T' */
            if (lower)
            {
                /* N is odd, TRANSR = 'T', and UPLO = 'L' */
                if (notrans)
                {
                    /* N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'N' */
                    dsyrk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[1], &n1);
                    dsyrk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, beta, &c__[2], &n1);
                    dgemm_("N", "T", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, &a[n1 + 1 + a_dim1], lda, beta, &c__[n1 * n1 + 1], &n1);
                }
                else
                {
                    /* N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'T' */
                    dsyrk_("U", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[1], &n1);
                    dsyrk_("L", "T", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1], lda, beta, &c__[2], &n1);
                    dgemm_("T", "N", &n1, &n2, k, alpha, &a[a_dim1 + 1], lda, &a[(n1 + 1) * a_dim1 + 1], lda, beta, &c__[n1 * n1 + 1], &n1);
                }
            }
            else
            {
                /* N is odd, TRANSR = 'T', and UPLO = 'U' */
                if (notrans)
                {
                    /* N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'N' */
                    dsyrk_("U", "N", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[n2 * n2 + 1], &n2);
                    dsyrk_("L", "N", &n2, k, alpha, &a[n1 + 1 + a_dim1], lda, beta, &c__[n1 * n2 + 1], &n2);
                    dgemm_("N", "T", &n2, &n1, k, alpha, &a[n1 + 1 + a_dim1], lda, &a[a_dim1 + 1], lda, beta, &c__[1], &n2);
                }
                else
                {
                    /* N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'T' */
                    dsyrk_("U", "T", &n1, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[n2 * n2 + 1], &n2);
                    dsyrk_("L", "T", &n2, k, alpha, &a[(n1 + 1) * a_dim1 + 1], lda, beta, &c__[n1 * n2 + 1], &n2);
                    dgemm_("T", "N", &n2, &n1, k, alpha, &a[(n1 + 1) * a_dim1 + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[1], & n2);
                }
            }
        }
    }
    else
    {
        /* N is even */
        if (normaltransr)
        {
            /* N is even and TRANSR = 'N' */
            if (lower)
            {
                /* N is even, TRANSR = 'N', and UPLO = 'L' */
                if (notrans)
                {
                    /* N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N' */
                    i__1 = *n + 1;
                    dsyrk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[2], &i__1);
                    i__1 = *n + 1;
                    dsyrk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, beta, &c__[1], &i__1);
                    i__1 = *n + 1;
                    dgemm_("N", "T", &nk, &nk, k, alpha, &a[nk + 1 + a_dim1], lda, &a[a_dim1 + 1], lda, beta, &c__[nk + 2], & i__1);
                }
                else
                {
                    /* N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'T' */
                    i__1 = *n + 1;
                    dsyrk_("L", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[2], &i__1);
                    i__1 = *n + 1;
                    dsyrk_("U", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[1], &i__1);
                    i__1 = *n + 1;
                    dgemm_("T", "N", &nk, &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[nk + 2] , &i__1);
                }
            }
            else
            {
                /* N is even, TRANSR = 'N', and UPLO = 'U' */
                if (notrans)
                {
                    /* N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N' */
                    i__1 = *n + 1;
                    dsyrk_("L", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[nk + 2], &i__1);
                    i__1 = *n + 1;
                    dsyrk_("U", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, beta, &c__[nk + 1], &i__1);
                    i__1 = *n + 1;
                    dgemm_("N", "T", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, &a[nk + 1 + a_dim1], lda, beta, &c__[1], &i__1);
                }
                else
                {
                    /* N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'T' */
                    i__1 = *n + 1;
                    dsyrk_("L", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[nk + 2], &i__1);
                    i__1 = *n + 1;
                    dsyrk_("U", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[nk + 1], &i__1);
                    i__1 = *n + 1;
                    dgemm_("T", "N", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[1], & i__1);
                }
            }
        }
        else
        {
            /* N is even, and TRANSR = 'T' */
            if (lower)
            {
                /* N is even, TRANSR = 'T', and UPLO = 'L' */
                if (notrans)
                {
                    /* N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'N' */
                    dsyrk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[nk + 1], &nk);
                    dsyrk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, beta, &c__[1], &nk);
                    dgemm_("N", "T", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, &a[nk + 1 + a_dim1], lda, beta, &c__[(nk + 1) * nk + 1], &nk);
                }
                else
                {
                    /* N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'T' */
                    dsyrk_("U", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[nk + 1], &nk);
                    dsyrk_("L", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[1], &nk);
                    dgemm_("T", "N", &nk, &nk, k, alpha, &a[a_dim1 + 1], lda, &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[(nk + 1) * nk + 1], &nk);
                }
            }
            else
            {
                /* N is even, TRANSR = 'T', and UPLO = 'U' */
                if (notrans)
                {
                    /* N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'N' */
                    dsyrk_("U", "N", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[nk * (nk + 1) + 1], &nk);
                    dsyrk_("L", "N", &nk, k, alpha, &a[nk + 1 + a_dim1], lda, beta, &c__[nk * nk + 1], &nk);
                    dgemm_("N", "T", &nk, &nk, k, alpha, &a[nk + 1 + a_dim1], lda, &a[a_dim1 + 1], lda, beta, &c__[1], &nk);
                }
                else
                {
                    /* N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'T' */
                    dsyrk_("U", "T", &nk, k, alpha, &a[a_dim1 + 1], lda, beta, &c__[nk * (nk + 1) + 1], &nk);
                    dsyrk_("L", "T", &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1], lda, beta, &c__[nk * nk + 1], &nk);
                    dgemm_("T", "N", &nk, &nk, k, alpha, &a[(nk + 1) * a_dim1 + 1], lda, &a[a_dim1 + 1], lda, beta, &c__[1], & nk);
                }
            }
        }
    }
    return 0;
    /* End of DSFRK */
}
/* dsfrk_ */
