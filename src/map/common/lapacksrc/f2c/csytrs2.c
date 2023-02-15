/* ../netlib/csytrs2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    1.f,0.f
}
;
/* > \brief \b CSYTRS2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSYTRS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrs2 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrs2 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrs2 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
/* WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYTRS2 solves a system of linear equations A*X = B with a COMPLEX */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by CSYTRF and converted by CSYCONV. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**T;
*/
/* > = 'L': Lower triangular, form is A = L*D*L**T. */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > The block diagonal matrix D and the multipliers used to */
/* > obtain the factor U or L as computed by CSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by CSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > On entry, the right hand side matrix B. */
/* > On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N) */
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
/* > \date November 2011 */
/* > \ingroup complexSYcomputational */
/* ===================================================================== */
/* Subroutine */
int csytrs2_(char *uplo, integer *n, integer *nrhs, complex * a, integer *lda, integer *ipiv, complex *b, integer *ldb, complex * work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    integer i__, j, k;
    complex ak, bk;
    integer kp;
    complex akm1, bkm1, akm1k;
    extern /* Subroutine */
    int cscal_(integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *);
    complex denom;
    integer iinfo;
    extern /* Subroutine */
    int cswap_(integer *, complex *, integer *, complex *, integer *), ctrsm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *);
    logical upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *), csyconv_( char *, char *, integer *, complex *, integer *, integer *, complex *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -8;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CSYTRS2", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *nrhs == 0)
    {
        return 0;
    }
    /* Convert A */
    csyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo);
    if (upper)
    {
        /* Solve A*X = B, where A = U*D*U**T. */
        /* P**T * B */
        k = *n;
        while(k >= 1)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal block */
                /* Interchange rows K and IPIV(K). */
                kp = ipiv[k];
                if (kp != k)
                {
                    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                --k;
            }
            else
            {
                /* 2 x 2 diagonal block */
                /* Interchange rows K-1 and -IPIV(K). */
                kp = -ipiv[k];
                if (kp == -ipiv[k - 1])
                {
                    cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += -2;
            }
        }
        /* Compute (U \P**T * B) -> B [ (U \P**T * B) ] */
        ctrsm_("L", "U", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[ b_offset], ldb);
        /* Compute D \ B -> B [ D \ (U \P**T * B) ] */
        i__ = *n;
        while(i__ >= 1)
        {
            if (ipiv[i__] > 0)
            {
                c_div(&q__1, &c_b1, &a[i__ + i__ * a_dim1]);
                cscal_(nrhs, &q__1, &b[i__ + b_dim1], ldb);
            }
            else if (i__ > 1)
            {
                if (ipiv[i__ - 1] == ipiv[i__])
                {
                    i__1 = i__;
                    akm1k.r = work[i__1].r;
                    akm1k.i = work[i__1].i; // , expr subst
                    c_div(&q__1, &a[i__ - 1 + (i__ - 1) * a_dim1], &akm1k);
                    akm1.r = q__1.r;
                    akm1.i = q__1.i; // , expr subst
                    c_div(&q__1, &a[i__ + i__ * a_dim1], &akm1k);
                    ak.r = q__1.r;
                    ak.i = q__1.i; // , expr subst
                    q__2.r = akm1.r * ak.r - akm1.i * ak.i;
                    q__2.i = akm1.r * ak.i + akm1.i * ak.r; // , expr subst
                    q__1.r = q__2.r - 1.f;
                    q__1.i = q__2.i - 0.f; // , expr subst
                    denom.r = q__1.r;
                    denom.i = q__1.i; // , expr subst
                    i__1 = *nrhs;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        c_div(&q__1, &b[i__ - 1 + j * b_dim1], &akm1k);
                        bkm1.r = q__1.r;
                        bkm1.i = q__1.i; // , expr subst
                        c_div(&q__1, &b[i__ + j * b_dim1], &akm1k);
                        bk.r = q__1.r;
                        bk.i = q__1.i; // , expr subst
                        i__2 = i__ - 1 + j * b_dim1;
                        q__3.r = ak.r * bkm1.r - ak.i * bkm1.i;
                        q__3.i = ak.r * bkm1.i + ak.i * bkm1.r; // , expr subst
                        q__2.r = q__3.r - bk.r;
                        q__2.i = q__3.i - bk.i; // , expr subst
                        c_div(&q__1, &q__2, &denom);
                        b[i__2].r = q__1.r;
                        b[i__2].i = q__1.i; // , expr subst
                        i__2 = i__ + j * b_dim1;
                        q__3.r = akm1.r * bk.r - akm1.i * bk.i;
                        q__3.i = akm1.r * bk.i + akm1.i * bk.r; // , expr subst
                        q__2.r = q__3.r - bkm1.r;
                        q__2.i = q__3.i - bkm1.i; // , expr subst
                        c_div(&q__1, &q__2, &denom);
                        b[i__2].r = q__1.r;
                        b[i__2].i = q__1.i; // , expr subst
                        /* L15: */
                    }
                    --i__;
                }
            }
            --i__;
        }
        /* Compute (U**T \ B) -> B [ U**T \ (D \ (U \P**T * B) ) ] */
        ctrsm_("L", "U", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[ b_offset], ldb);
        /* P * B [ P * (U**T \ (D \ (U \P**T * B) )) ] */
        k = 1;
        while(k <= *n)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal block */
                /* Interchange rows K and IPIV(K). */
                kp = ipiv[k];
                if (kp != k)
                {
                    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                ++k;
            }
            else
            {
                /* 2 x 2 diagonal block */
                /* Interchange rows K-1 and -IPIV(K). */
                kp = -ipiv[k];
                if (k < *n && kp == -ipiv[k + 1])
                {
                    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += 2;
            }
        }
    }
    else
    {
        /* Solve A*X = B, where A = L*D*L**T. */
        /* P**T * B */
        k = 1;
        while(k <= *n)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal block */
                /* Interchange rows K and IPIV(K). */
                kp = ipiv[k];
                if (kp != k)
                {
                    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                ++k;
            }
            else
            {
                /* 2 x 2 diagonal block */
                /* Interchange rows K and -IPIV(K+1). */
                kp = -ipiv[k + 1];
                if (kp == -ipiv[k])
                {
                    cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += 2;
            }
        }
        /* Compute (L \P**T * B) -> B [ (L \P**T * B) ] */
        ctrsm_("L", "L", "N", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[ b_offset], ldb);
        /* Compute D \ B -> B [ D \ (L \P**T * B) ] */
        i__ = 1;
        while(i__ <= *n)
        {
            if (ipiv[i__] > 0)
            {
                c_div(&q__1, &c_b1, &a[i__ + i__ * a_dim1]);
                cscal_(nrhs, &q__1, &b[i__ + b_dim1], ldb);
            }
            else
            {
                i__1 = i__;
                akm1k.r = work[i__1].r;
                akm1k.i = work[i__1].i; // , expr subst
                c_div(&q__1, &a[i__ + i__ * a_dim1], &akm1k);
                akm1.r = q__1.r;
                akm1.i = q__1.i; // , expr subst
                c_div(&q__1, &a[i__ + 1 + (i__ + 1) * a_dim1], &akm1k);
                ak.r = q__1.r;
                ak.i = q__1.i; // , expr subst
                q__2.r = akm1.r * ak.r - akm1.i * ak.i;
                q__2.i = akm1.r * ak.i + akm1.i * ak.r; // , expr subst
                q__1.r = q__2.r - 1.f;
                q__1.i = q__2.i - 0.f; // , expr subst
                denom.r = q__1.r;
                denom.i = q__1.i; // , expr subst
                i__1 = *nrhs;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    c_div(&q__1, &b[i__ + j * b_dim1], &akm1k);
                    bkm1.r = q__1.r;
                    bkm1.i = q__1.i; // , expr subst
                    c_div(&q__1, &b[i__ + 1 + j * b_dim1], &akm1k);
                    bk.r = q__1.r;
                    bk.i = q__1.i; // , expr subst
                    i__2 = i__ + j * b_dim1;
                    q__3.r = ak.r * bkm1.r - ak.i * bkm1.i;
                    q__3.i = ak.r * bkm1.i + ak.i * bkm1.r; // , expr subst
                    q__2.r = q__3.r - bk.r;
                    q__2.i = q__3.i - bk.i; // , expr subst
                    c_div(&q__1, &q__2, &denom);
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    i__2 = i__ + 1 + j * b_dim1;
                    q__3.r = akm1.r * bk.r - akm1.i * bk.i;
                    q__3.i = akm1.r * bk.i + akm1.i * bk.r; // , expr subst
                    q__2.r = q__3.r - bkm1.r;
                    q__2.i = q__3.i - bkm1.i; // , expr subst
                    c_div(&q__1, &q__2, &denom);
                    b[i__2].r = q__1.r;
                    b[i__2].i = q__1.i; // , expr subst
                    /* L25: */
                }
                ++i__;
            }
            ++i__;
        }
        /* Compute (L**T \ B) -> B [ L**T \ (D \ (L \P**T * B) ) ] */
        ctrsm_("L", "L", "T", "U", n, nrhs, &c_b1, &a[a_offset], lda, &b[ b_offset], ldb);
        /* P * B [ P * (L**T \ (D \ (L \P**T * B) )) ] */
        k = *n;
        while(k >= 1)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal block */
                /* Interchange rows K and IPIV(K). */
                kp = ipiv[k];
                if (kp != k)
                {
                    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                --k;
            }
            else
            {
                /* 2 x 2 diagonal block */
                /* Interchange rows K-1 and -IPIV(K). */
                kp = -ipiv[k];
                if (k > 1 && kp == -ipiv[k - 1])
                {
                    cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += -2;
            }
        }
    }
    /* Revert A */
    csyconv_(uplo, "R", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo);
    return 0;
    /* End of CSYTRS2 */
}
/* csytrs2_ */
