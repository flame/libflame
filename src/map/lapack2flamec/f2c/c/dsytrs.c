/* ../netlib/dsytrs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b7 = -1.;
static integer c__1 = 1;
static doublereal c_b19 = 1.;
/* > \brief \b DSYTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSYTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytrs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytrs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytrs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* DOUBLE PRECISION A( LDA, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYTRS solves a system of linear equations A*X = B with a real */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by DSYTRF. */
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
int dsytrs_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer * ldb, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1;
    /* Local variables */
    integer j, k;
    doublereal ak, bk;
    integer kp;
    doublereal akm1, bkm1;
    extern /* Subroutine */
    int dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal akm1k;
    extern /* Subroutine */
    int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *);
    doublereal denom;
    extern /* Subroutine */
    int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
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
        xerbla_("DSYTRS", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0 || *nrhs == 0)
    {
        return 0;
    }
    if (upper)
    {
        /* Solve A*X = B, where A = U*D*U**T. */
        /* First solve U*D*X = B, overwriting B with X. */
        /* K is the main loop index, decreasing from N to 1 in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        k = *n;
L10: /* If K < 1, exit from loop. */
        if (k < 1)
        {
            goto L30;
        }
        if (ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Interchange rows K and IPIV(K). */
            kp = ipiv[k];
            if (kp != k)
            {
                dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(U(K)), where U(K) is the transformation */
            /* stored in column K of A. */
            i__1 = k - 1;
            dger_(&i__1, nrhs, &c_b7, &a[k * a_dim1 + 1], &c__1, &b[k + b_dim1], ldb, &b[b_dim1 + 1], ldb);
            /* Multiply by the inverse of the diagonal block. */
            d__1 = 1. / a[k + k * a_dim1];
            dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
            --k;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Interchange rows K-1 and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k - 1)
            {
                dswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(U(K)), where U(K) is the transformation */
            /* stored in columns K-1 and K of A. */
            i__1 = k - 2;
            dger_(&i__1, nrhs, &c_b7, &a[k * a_dim1 + 1], &c__1, &b[k + b_dim1], ldb, &b[b_dim1 + 1], ldb);
            i__1 = k - 2;
            dger_(&i__1, nrhs, &c_b7, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);
            /* Multiply by the inverse of the diagonal block. */
            akm1k = a[k - 1 + k * a_dim1];
            akm1 = a[k - 1 + (k - 1) * a_dim1] / akm1k;
            ak = a[k + k * a_dim1] / akm1k;
            denom = akm1 * ak - 1.;
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                bkm1 = b[k - 1 + j * b_dim1] / akm1k;
                bk = b[k + j * b_dim1] / akm1k;
                b[k - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
                b[k + j * b_dim1] = (akm1 * bk - bkm1) / denom;
                /* L20: */
            }
            k += -2;
        }
        goto L10;
L30: /* Next solve U**T *X = B, overwriting B with X. */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        k = 1;
L40: /* If K > N, exit from loop. */
        if (k > *n)
        {
            goto L50;
        }
        if (ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Multiply by inv(U**T(K)), where U(K) is the transformation */
            /* stored in column K of A. */
            i__1 = k - 1;
            dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[k * a_dim1 + 1], &c__1, &c_b19, &b[k + b_dim1], ldb);
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
            /* Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
            /* stored in columns K and K+1 of A. */
            i__1 = k - 1;
            dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[k * a_dim1 + 1], &c__1, &c_b19, &b[k + b_dim1], ldb);
            i__1 = k - 1;
            dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[(k + 1) * a_dim1 + 1], &c__1, &c_b19, &b[k + 1 + b_dim1], ldb);
            /* Interchange rows K and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k)
            {
                dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            k += 2;
        }
        goto L40;
L50:
        ;
    }
    else
    {
        /* Solve A*X = B, where A = L*D*L**T. */
        /* First solve L*D*X = B, overwriting B with X. */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        k = 1;
L60: /* If K > N, exit from loop. */
        if (k > *n)
        {
            goto L80;
        }
        if (ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Interchange rows K and IPIV(K). */
            kp = ipiv[k];
            if (kp != k)
            {
                dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(L(K)), where L(K) is the transformation */
            /* stored in column K of A. */
            if (k < *n)
            {
                i__1 = *n - k;
                dger_(&i__1, nrhs, &c_b7, &a[k + 1 + k * a_dim1], &c__1, &b[k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
            }
            /* Multiply by the inverse of the diagonal block. */
            d__1 = 1. / a[k + k * a_dim1];
            dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
            ++k;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Interchange rows K+1 and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k + 1)
            {
                dswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(L(K)), where L(K) is the transformation */
            /* stored in columns K and K+1 of A. */
            if (k < *n - 1)
            {
                i__1 = *n - k - 1;
                dger_(&i__1, nrhs, &c_b7, &a[k + 2 + k * a_dim1], &c__1, &b[k + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
                i__1 = *n - k - 1;
                dger_(&i__1, nrhs, &c_b7, &a[k + 2 + (k + 1) * a_dim1], &c__1, &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
            }
            /* Multiply by the inverse of the diagonal block. */
            akm1k = a[k + 1 + k * a_dim1];
            akm1 = a[k + k * a_dim1] / akm1k;
            ak = a[k + 1 + (k + 1) * a_dim1] / akm1k;
            denom = akm1 * ak - 1.;
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                bkm1 = b[k + j * b_dim1] / akm1k;
                bk = b[k + 1 + j * b_dim1] / akm1k;
                b[k + j * b_dim1] = (ak * bkm1 - bk) / denom;
                b[k + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
                /* L70: */
            }
            k += 2;
        }
        goto L60;
L80: /* Next solve L**T *X = B, overwriting B with X. */
        /* K is the main loop index, decreasing from N to 1 in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        k = *n;
L90: /* If K < 1, exit from loop. */
        if (k < 1)
        {
            goto L100;
        }
        if (ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Multiply by inv(L**T(K)), where L(K) is the transformation */
            /* stored in column K of A. */
            if (k < *n)
            {
                i__1 = *n - k;
                dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b19, &b[k + b_dim1], ldb);
            }
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
            /* Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
            /* stored in columns K-1 and K of A. */
            if (k < *n)
            {
                i__1 = *n - k;
                dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b19, &b[k + b_dim1], ldb);
                i__1 = *n - k;
                dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b19, &b[ k - 1 + b_dim1], ldb);
            }
            /* Interchange rows K and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k)
            {
                dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            k += -2;
        }
        goto L90;
L100:
        ;
    }
    return 0;
    /* End of DSYTRS */
}
/* dsytrs_ */
