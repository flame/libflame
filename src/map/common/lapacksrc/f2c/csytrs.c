/* ../netlib/csytrs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    1.f,0.f
}
;
static integer c__1 = 1;
/* > \brief \b CSYTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSYTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYTRS solves a system of linear equations A*X = B with a complex */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by CSYTRF. */
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
int csytrs_(char *uplo, integer *n, integer *nrhs, complex * a, integer *lda, integer *ipiv, complex *b, integer *ldb, integer * info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    integer j, k;
    complex ak, bk;
    integer kp;
    complex akm1, bkm1, akm1k;
    extern /* Subroutine */
    int cscal_(integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *);
    complex denom;
    extern /* Subroutine */
    int cgemv_(char *, integer *, integer *, complex * , complex *, integer *, complex *, integer *, complex *, complex * , integer *), cgeru_(integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, integer *), cswap_(integer *, complex *, integer *, complex *, integer *);
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
        xerbla_("CSYTRS", &i__1);
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
                cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(U(K)), where U(K) is the transformation */
            /* stored in column K of A. */
            i__1 = k - 1;
            q__1.r = -1.f;
            q__1.i = -0.f; // , expr subst
            cgeru_(&i__1, nrhs, &q__1, &a[k * a_dim1 + 1], &c__1, &b[k + b_dim1], ldb, &b[b_dim1 + 1], ldb);
            /* Multiply by the inverse of the diagonal block. */
            c_div(&q__1, &c_b1, &a[k + k * a_dim1]);
            cscal_(nrhs, &q__1, &b[k + b_dim1], ldb);
            --k;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Interchange rows K-1 and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k - 1)
            {
                cswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(U(K)), where U(K) is the transformation */
            /* stored in columns K-1 and K of A. */
            i__1 = k - 2;
            q__1.r = -1.f;
            q__1.i = -0.f; // , expr subst
            cgeru_(&i__1, nrhs, &q__1, &a[k * a_dim1 + 1], &c__1, &b[k + b_dim1], ldb, &b[b_dim1 + 1], ldb);
            i__1 = k - 2;
            q__1.r = -1.f;
            q__1.i = -0.f; // , expr subst
            cgeru_(&i__1, nrhs, &q__1, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);
            /* Multiply by the inverse of the diagonal block. */
            i__1 = k - 1 + k * a_dim1;
            akm1k.r = a[i__1].r;
            akm1k.i = a[i__1].i; // , expr subst
            c_div(&q__1, &a[k - 1 + (k - 1) * a_dim1], &akm1k);
            akm1.r = q__1.r;
            akm1.i = q__1.i; // , expr subst
            c_div(&q__1, &a[k + k * a_dim1], &akm1k);
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
                c_div(&q__1, &b[k - 1 + j * b_dim1], &akm1k);
                bkm1.r = q__1.r;
                bkm1.i = q__1.i; // , expr subst
                c_div(&q__1, &b[k + j * b_dim1], &akm1k);
                bk.r = q__1.r;
                bk.i = q__1.i; // , expr subst
                i__2 = k - 1 + j * b_dim1;
                q__3.r = ak.r * bkm1.r - ak.i * bkm1.i;
                q__3.i = ak.r * bkm1.i + ak.i * bkm1.r; // , expr subst
                q__2.r = q__3.r - bk.r;
                q__2.i = q__3.i - bk.i; // , expr subst
                c_div(&q__1, &q__2, &denom);
                b[i__2].r = q__1.r;
                b[i__2].i = q__1.i; // , expr subst
                i__2 = k + j * b_dim1;
                q__3.r = akm1.r * bk.r - akm1.i * bk.i;
                q__3.i = akm1.r * bk.i + akm1.i * bk.r; // , expr subst
                q__2.r = q__3.r - bkm1.r;
                q__2.i = q__3.i - bkm1.i; // , expr subst
                c_div(&q__1, &q__2, &denom);
                b[i__2].r = q__1.r;
                b[i__2].i = q__1.i; // , expr subst
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
            q__1.r = -1.f;
            q__1.i = -0.f; // , expr subst
            cgemv_("Transpose", &i__1, nrhs, &q__1, &b[b_offset], ldb, &a[k * a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb) ;
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
            /* Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
            /* stored in columns K and K+1 of A. */
            i__1 = k - 1;
            q__1.r = -1.f;
            q__1.i = -0.f; // , expr subst
            cgemv_("Transpose", &i__1, nrhs, &q__1, &b[b_offset], ldb, &a[k * a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb) ;
            i__1 = k - 1;
            q__1.r = -1.f;
            q__1.i = -0.f; // , expr subst
            cgemv_("Transpose", &i__1, nrhs, &q__1, &b[b_offset], ldb, &a[(k + 1) * a_dim1 + 1], &c__1, &c_b1, &b[k + 1 + b_dim1], ldb);
            /* Interchange rows K and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k)
            {
                cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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
                cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(L(K)), where L(K) is the transformation */
            /* stored in column K of A. */
            if (k < *n)
            {
                i__1 = *n - k;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgeru_(&i__1, nrhs, &q__1, &a[k + 1 + k * a_dim1], &c__1, &b[ k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
            }
            /* Multiply by the inverse of the diagonal block. */
            c_div(&q__1, &c_b1, &a[k + k * a_dim1]);
            cscal_(nrhs, &q__1, &b[k + b_dim1], ldb);
            ++k;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Interchange rows K+1 and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k + 1)
            {
                cswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(L(K)), where L(K) is the transformation */
            /* stored in columns K and K+1 of A. */
            if (k < *n - 1)
            {
                i__1 = *n - k - 1;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgeru_(&i__1, nrhs, &q__1, &a[k + 2 + k * a_dim1], &c__1, &b[ k + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
                i__1 = *n - k - 1;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgeru_(&i__1, nrhs, &q__1, &a[k + 2 + (k + 1) * a_dim1], & c__1, &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
            }
            /* Multiply by the inverse of the diagonal block. */
            i__1 = k + 1 + k * a_dim1;
            akm1k.r = a[i__1].r;
            akm1k.i = a[i__1].i; // , expr subst
            c_div(&q__1, &a[k + k * a_dim1], &akm1k);
            akm1.r = q__1.r;
            akm1.i = q__1.i; // , expr subst
            c_div(&q__1, &a[k + 1 + (k + 1) * a_dim1], &akm1k);
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
                c_div(&q__1, &b[k + j * b_dim1], &akm1k);
                bkm1.r = q__1.r;
                bkm1.i = q__1.i; // , expr subst
                c_div(&q__1, &b[k + 1 + j * b_dim1], &akm1k);
                bk.r = q__1.r;
                bk.i = q__1.i; // , expr subst
                i__2 = k + j * b_dim1;
                q__3.r = ak.r * bkm1.r - ak.i * bkm1.i;
                q__3.i = ak.r * bkm1.i + ak.i * bkm1.r; // , expr subst
                q__2.r = q__3.r - bk.r;
                q__2.i = q__3.i - bk.i; // , expr subst
                c_div(&q__1, &q__2, &denom);
                b[i__2].r = q__1.r;
                b[i__2].i = q__1.i; // , expr subst
                i__2 = k + 1 + j * b_dim1;
                q__3.r = akm1.r * bk.r - akm1.i * bk.i;
                q__3.i = akm1.r * bk.i + akm1.i * bk.r; // , expr subst
                q__2.r = q__3.r - bkm1.r;
                q__2.i = q__3.i - bkm1.i; // , expr subst
                c_div(&q__1, &q__2, &denom);
                b[i__2].r = q__1.r;
                b[i__2].i = q__1.i; // , expr subst
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
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgemv_("Transpose", &i__1, nrhs, &q__1, &b[k + 1 + b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + b_dim1], ldb);
            }
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
            /* Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
            /* stored in columns K-1 and K of A. */
            if (k < *n)
            {
                i__1 = *n - k;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgemv_("Transpose", &i__1, nrhs, &q__1, &b[k + 1 + b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + b_dim1], ldb);
                i__1 = *n - k;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgemv_("Transpose", &i__1, nrhs, &q__1, &b[k + 1 + b_dim1], ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b1, &b[k - 1 + b_dim1], ldb);
            }
            /* Interchange rows K and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k)
            {
                cswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            k += -2;
        }
        goto L90;
L100:
        ;
    }
    return 0;
    /* End of CSYTRS */
}
/* csytrs_ */
