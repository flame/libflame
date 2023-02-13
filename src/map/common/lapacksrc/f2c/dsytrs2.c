/* ../netlib/dsytrs2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b10 = 1.;
/* > \brief \b DSYTRS2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSYTRS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrs2 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrs2 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrs2 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, */
/* WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* DOUBLE PRECISION A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTRS2 solves a system of linear equations A*X = B with a real */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by DSYTRF and converted by DSYCONV. */
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
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > The block diagonal matrix D and the multipliers used to */
/* > obtain the factor U or L as computed by DSYTRF. */
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
/* > as determined by DSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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
/* > WORK is REAL array, dimension (N) */
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
/* > \ingroup doubleSYcomputational */
/* ===================================================================== */
/* Subroutine */
int dsytrs2_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer * ldb, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1;
    /* Local variables */
    integer i__, j, k;
    doublereal ak, bk;
    integer kp;
    doublereal akm1, bkm1, akm1k;
    extern /* Subroutine */
    int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *);
    doublereal denom;
    integer iinfo;
    extern /* Subroutine */
    int dswap_(integer *, doublereal *, integer *, doublereal *, integer *), dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    logical upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *), dsyconv_( char *, char *, integer *, doublereal *, integer *, integer *, doublereal *, integer *);
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
        xerbla_("DSYTRS2", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *nrhs == 0)
    {
        return 0;
    }
    /* Convert A */
    dsyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo);
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
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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
                    dswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += -2;
            }
        }
        /* Compute (U \P**T * B) -> B [ (U \P**T * B) ] */
        dtrsm_("L", "U", "N", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[ b_offset], ldb);
        /* Compute D \ B -> B [ D \ (U \P**T * B) ] */
        i__ = *n;
        while(i__ >= 1)
        {
            if (ipiv[i__] > 0)
            {
                d__1 = 1. / a[i__ + i__ * a_dim1];
                dscal_(nrhs, &d__1, &b[i__ + b_dim1], ldb);
            }
            else if (i__ > 1)
            {
                if (ipiv[i__ - 1] == ipiv[i__])
                {
                    akm1k = work[i__];
                    akm1 = a[i__ - 1 + (i__ - 1) * a_dim1] / akm1k;
                    ak = a[i__ + i__ * a_dim1] / akm1k;
                    denom = akm1 * ak - 1.;
                    i__1 = *nrhs;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        bkm1 = b[i__ - 1 + j * b_dim1] / akm1k;
                        bk = b[i__ + j * b_dim1] / akm1k;
                        b[i__ - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
                        b[i__ + j * b_dim1] = (akm1 * bk - bkm1) / denom;
                        /* L15: */
                    }
                    --i__;
                }
            }
            --i__;
        }
        /* Compute (U**T \ B) -> B [ U**T \ (D \ (U \P**T * B) ) ] */
        dtrsm_("L", "U", "T", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[ b_offset], ldb);
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
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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
                    dswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += 2;
            }
        }
        /* Compute (L \P**T * B) -> B [ (L \P**T * B) ] */
        dtrsm_("L", "L", "N", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[ b_offset], ldb);
        /* Compute D \ B -> B [ D \ (L \P**T * B) ] */
        i__ = 1;
        while(i__ <= *n)
        {
            if (ipiv[i__] > 0)
            {
                d__1 = 1. / a[i__ + i__ * a_dim1];
                dscal_(nrhs, &d__1, &b[i__ + b_dim1], ldb);
            }
            else
            {
                akm1k = work[i__];
                akm1 = a[i__ + i__ * a_dim1] / akm1k;
                ak = a[i__ + 1 + (i__ + 1) * a_dim1] / akm1k;
                denom = akm1 * ak - 1.;
                i__1 = *nrhs;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    bkm1 = b[i__ + j * b_dim1] / akm1k;
                    bk = b[i__ + 1 + j * b_dim1] / akm1k;
                    b[i__ + j * b_dim1] = (ak * bkm1 - bk) / denom;
                    b[i__ + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
                    /* L25: */
                }
                ++i__;
            }
            ++i__;
        }
        /* Compute (L**T \ B) -> B [ L**T \ (D \ (L \P**T * B) ) ] */
        dtrsm_("L", "L", "T", "U", n, nrhs, &c_b10, &a[a_offset], lda, &b[ b_offset], ldb);
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
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += -2;
            }
        }
    }
    /* Revert A */
    dsyconv_(uplo, "R", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo);
    return 0;
    /* End of DSYTRS2 */
}
/* dsytrs2_ */
