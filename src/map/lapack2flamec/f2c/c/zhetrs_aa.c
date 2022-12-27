/* ../netlib/v3.9.0/zhetrs_aa.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b9 =
{
    1.,0.
}
;
static integer c__1 = 1;
/* > \brief \b ZHETRS_AA */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHETRS_AA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs_ aa.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs_ aa.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs_ aa.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
/* WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER N, NRHS, LDA, LDB, LWORK, INFO */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRS_AA solves a system of linear equations A*X = B with a complex */
/* > hermitian matrix A using the factorization A = U**H*T*U or */
/* > A = L*T*L**H computed by ZHETRF_AA. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U**H*T*U;
*/
/* > = 'L': Lower triangular, form is A = L*T*L**H. */
/* > \endverbatim */
/* > */
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
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > Details of factors computed by ZHETRF_AA. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges as computed by ZHETRF_AA. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > On entry, the right hand side matrix B. */
/* > On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,3*N-2). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup complex16HEcomputational */
/* ===================================================================== */
/* Subroutine */
int zhetrs_aa_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhetrs_aa inputs: uplo %c, n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "",*uplo, *n, *nrhs, *lda, *ldb);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    /* Local variables */
    integer k, kp;
    extern logical lsame_(char *, char *);
    logical upper;
    extern /* Subroutine */
    int zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), zgtsv_(integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *, integer *, integer *), ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(char *, integer *), zlacgv_(integer *, doublecomplex *, integer *), zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    integer lwkopt;
    logical lquery;
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    lquery = *lwork == -1;
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*nrhs < 0)
    {
        *info = -3;
    }
    else if (*lda < fla_max(1,*n))
    {
        *info = -5;
    }
    else if (*ldb < fla_max(1,*n))
    {
        *info = -8;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n * 3 - 2; // , expr subst
        if (*lwork < fla_max(i__1,i__2) && ! lquery)
        {
            *info = -10;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHETRS_AA", &i__1);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (lquery)
    {
        lwkopt = *n * 3 - 2;
        work[1].r = (doublereal) lwkopt;
        work[1].i = 0.; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *nrhs == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (upper)
    {
        /* Solve A*X = B, where A = U**H*T*U. */
        /* 1) Forward substitution with U**H */
        if (*n > 1)
        {
            /* Pivot, P**T * B -> B */
            i__1 = *n;
            for (k = 1;
                    k <= i__1;
                    ++k)
            {
                kp = ipiv[k];
                if (kp != k)
                {
                    zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
            }
            /* Compute U**H \ B -> B [ (U**H \P**T * B) ] */
            i__1 = *n - 1;
            ztrsm_("L", "U", "C", "U", &i__1, nrhs, &c_b9, &a[(a_dim1 << 1) + 1], lda, &b[b_dim1 + 2], ldb);
        }
        /* 2) Solve with triangular matrix T */
        /* Compute T \ B -> B [ T \ (U**H \P**T * B) ] */
        i__1 = *lda + 1;
        zlacpy_("F", &c__1, n, &a[a_dim1 + 1], &i__1, &work[*n], &c__1);
        if (*n > 1)
        {
            i__1 = *n - 1;
            i__2 = *lda + 1;
            zlacpy_("F", &c__1, &i__1, &a[(a_dim1 << 1) + 1], &i__2, &work[*n * 2], &c__1);
            i__1 = *n - 1;
            i__2 = *lda + 1;
            zlacpy_("F", &c__1, &i__1, &a[(a_dim1 << 1) + 1], &i__2, &work[1], &c__1);
            i__1 = *n - 1;
            zlacgv_(&i__1, &work[1], &c__1);
        }
        zgtsv_(n, nrhs, &work[1], &work[*n], &work[*n * 2], &b[b_offset], ldb, info);
        /* 3) Backward substitution with U */
        if (*n > 1)
        {
            /* Compute U \ B -> B [ U \ (T \ (U**H \P**T * B) ) ] */
            i__1 = *n - 1;
            ztrsm_("L", "U", "N", "U", &i__1, nrhs, &c_b9, &a[(a_dim1 << 1) + 1], lda, &b[b_dim1 + 2], ldb);
            /* Pivot, P * B [ P * (U**H \ (T \ (U \P**T * B) )) ] */
            for (k = *n;
                    k >= 1;
                    --k)
            {
                kp = ipiv[k];
                if (kp != k)
                {
                    zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
            }
        }
    }
    else
    {
        /* Solve A*X = B, where A = L*T*L**H. */
        /* 1) Forward substitution with L */
        if (*n > 1)
        {
            /* Pivot, P**T * B -> B */
            i__1 = *n;
            for (k = 1;
                    k <= i__1;
                    ++k)
            {
                kp = ipiv[k];
                if (kp != k)
                {
                    zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
            }
            /* Compute L \ B -> B [ (L \P**T * B) ] */
            i__1 = *n - 1;
            ztrsm_("L", "L", "N", "U", &i__1, nrhs, &c_b9, &a[a_dim1 + 2], lda, &b[b_dim1 + 2], ldb);
        }
        /* 2) Solve with triangular matrix T */
        /* Compute T \ B -> B [ T \ (L \P**T * B) ] */
        i__1 = *lda + 1;
        zlacpy_("F", &c__1, n, &a[a_dim1 + 1], &i__1, &work[*n], &c__1);
        if (*n > 1)
        {
            i__1 = *n - 1;
            i__2 = *lda + 1;
            zlacpy_("F", &c__1, &i__1, &a[a_dim1 + 2], &i__2, &work[1], &c__1);
            i__1 = *n - 1;
            i__2 = *lda + 1;
            zlacpy_("F", &c__1, &i__1, &a[a_dim1 + 2], &i__2, &work[*n * 2], & c__1);
            i__1 = *n - 1;
            zlacgv_(&i__1, &work[*n * 2], &c__1);
        }
        zgtsv_(n, nrhs, &work[1], &work[*n], &work[*n * 2], &b[b_offset], ldb, info);
        /* 3) Backward substitution with L**H */
        if (*n > 1)
        {
            /* Compute L**H \ B -> B [ L**H \ (T \ (L \P**T * B) ) ] */
            i__1 = *n - 1;
            ztrsm_("L", "L", "C", "U", &i__1, nrhs, &c_b9, &a[a_dim1 + 2], lda, &b[b_dim1 + 2], ldb);
            /* Pivot, P * B [ P * (L**H \ (T \ (L \P**T * B) )) ] */
            for (k = *n;
                    k >= 1;
                    --k)
            {
                kp = ipiv[k];
                if (kp != k)
                {
                    zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of ZHETRS_AA */
}
/* zhetrs_aa__ */

