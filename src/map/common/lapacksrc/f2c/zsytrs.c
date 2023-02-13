/* ../netlib/zsytrs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    1.,0.
}
;
static integer c__1 = 1;
/* > \brief \b ZSYTRS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZSYTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrs. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrs. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrs. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LDB, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYTRS solves a system of linear equations A*X = B with a complex */
/* > symmetric matrix A using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by ZSYTRF. */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > The block diagonal matrix D and the multipliers used to */
/* > obtain the factor U or L as computed by ZSYTRF. */
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
/* > as determined by ZSYTRF. */
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
/* > \ingroup complex16SYcomputational */
/* ===================================================================== */
/* Subroutine */
int zsytrs_(char *uplo, integer *n, integer *nrhs, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    /* Local variables */
    integer j, k;
    doublecomplex ak, bk;
    integer kp;
    doublecomplex akm1, bkm1, akm1k;
    extern logical lsame_(char *, char *);
    doublecomplex denom;
    extern /* Subroutine */
    int zscal_(integer *, doublecomplex *, doublecomplex *, integer *), zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
    logical upper;
    extern /* Subroutine */
    int zgeru_(integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *), zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(char *, integer *);
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
        xerbla_("ZSYTRS", &i__1);
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
                zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(U(K)), where U(K) is the transformation */
            /* stored in column K of A. */
            i__1 = k - 1;
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + b_dim1], ldb, &b[b_dim1 + 1], ldb);
            /* Multiply by the inverse of the diagonal block. */
            z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
            zscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
            --k;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Interchange rows K-1 and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k - 1)
            {
                zswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(U(K)), where U(K) is the transformation */
            /* stored in columns K-1 and K of A. */
            i__1 = k - 2;
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgeru_(&i__1, nrhs, &z__1, &a[k * a_dim1 + 1], &c__1, &b[k + b_dim1], ldb, &b[b_dim1 + 1], ldb);
            i__1 = k - 2;
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgeru_(&i__1, nrhs, &z__1, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k - 1 + b_dim1], ldb, &b[b_dim1 + 1], ldb);
            /* Multiply by the inverse of the diagonal block. */
            i__1 = k - 1 + k * a_dim1;
            akm1k.r = a[i__1].r;
            akm1k.i = a[i__1].i; // , expr subst
            z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &akm1k);
            akm1.r = z__1.r;
            akm1.i = z__1.i; // , expr subst
            z_div(&z__1, &a[k + k * a_dim1], &akm1k);
            ak.r = z__1.r;
            ak.i = z__1.i; // , expr subst
            z__2.r = akm1.r * ak.r - akm1.i * ak.i;
            z__2.i = akm1.r * ak.i + akm1.i * ak.r; // , expr subst
            z__1.r = z__2.r - 1.;
            z__1.i = z__2.i - 0.; // , expr subst
            denom.r = z__1.r;
            denom.i = z__1.i; // , expr subst
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                z_div(&z__1, &b[k - 1 + j * b_dim1], &akm1k);
                bkm1.r = z__1.r;
                bkm1.i = z__1.i; // , expr subst
                z_div(&z__1, &b[k + j * b_dim1], &akm1k);
                bk.r = z__1.r;
                bk.i = z__1.i; // , expr subst
                i__2 = k - 1 + j * b_dim1;
                z__3.r = ak.r * bkm1.r - ak.i * bkm1.i;
                z__3.i = ak.r * bkm1.i + ak.i * bkm1.r; // , expr subst
                z__2.r = z__3.r - bk.r;
                z__2.i = z__3.i - bk.i; // , expr subst
                z_div(&z__1, &z__2, &denom);
                b[i__2].r = z__1.r;
                b[i__2].i = z__1.i; // , expr subst
                i__2 = k + j * b_dim1;
                z__3.r = akm1.r * bk.r - akm1.i * bk.i;
                z__3.i = akm1.r * bk.i + akm1.i * bk.r; // , expr subst
                z__2.r = z__3.r - bkm1.r;
                z__2.i = z__3.i - bkm1.i; // , expr subst
                z_div(&z__1, &z__2, &denom);
                b[i__2].r = z__1.r;
                b[i__2].i = z__1.i; // , expr subst
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
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[k * a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb) ;
            /* Interchange rows K and IPIV(K). */
            kp = ipiv[k];
            if (kp != k)
            {
                zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            ++k;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
            /* stored in columns K and K+1 of A. */
            i__1 = k - 1;
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[k * a_dim1 + 1], &c__1, &c_b1, &b[k + b_dim1], ldb) ;
            i__1 = k - 1;
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgemv_("Transpose", &i__1, nrhs, &z__1, &b[b_offset], ldb, &a[(k + 1) * a_dim1 + 1], &c__1, &c_b1, &b[k + 1 + b_dim1], ldb);
            /* Interchange rows K and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k)
            {
                zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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
                zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(L(K)), where L(K) is the transformation */
            /* stored in column K of A. */
            if (k < *n)
            {
                i__1 = *n - k;
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgeru_(&i__1, nrhs, &z__1, &a[k + 1 + k * a_dim1], &c__1, &b[ k + b_dim1], ldb, &b[k + 1 + b_dim1], ldb);
            }
            /* Multiply by the inverse of the diagonal block. */
            z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
            zscal_(nrhs, &z__1, &b[k + b_dim1], ldb);
            ++k;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Interchange rows K+1 and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k + 1)
            {
                zswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            /* Multiply by inv(L(K)), where L(K) is the transformation */
            /* stored in columns K and K+1 of A. */
            if (k < *n - 1)
            {
                i__1 = *n - k - 1;
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgeru_(&i__1, nrhs, &z__1, &a[k + 2 + k * a_dim1], &c__1, &b[ k + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
                i__1 = *n - k - 1;
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgeru_(&i__1, nrhs, &z__1, &a[k + 2 + (k + 1) * a_dim1], & c__1, &b[k + 1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
            }
            /* Multiply by the inverse of the diagonal block. */
            i__1 = k + 1 + k * a_dim1;
            akm1k.r = a[i__1].r;
            akm1k.i = a[i__1].i; // , expr subst
            z_div(&z__1, &a[k + k * a_dim1], &akm1k);
            akm1.r = z__1.r;
            akm1.i = z__1.i; // , expr subst
            z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &akm1k);
            ak.r = z__1.r;
            ak.i = z__1.i; // , expr subst
            z__2.r = akm1.r * ak.r - akm1.i * ak.i;
            z__2.i = akm1.r * ak.i + akm1.i * ak.r; // , expr subst
            z__1.r = z__2.r - 1.;
            z__1.i = z__2.i - 0.; // , expr subst
            denom.r = z__1.r;
            denom.i = z__1.i; // , expr subst
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                z_div(&z__1, &b[k + j * b_dim1], &akm1k);
                bkm1.r = z__1.r;
                bkm1.i = z__1.i; // , expr subst
                z_div(&z__1, &b[k + 1 + j * b_dim1], &akm1k);
                bk.r = z__1.r;
                bk.i = z__1.i; // , expr subst
                i__2 = k + j * b_dim1;
                z__3.r = ak.r * bkm1.r - ak.i * bkm1.i;
                z__3.i = ak.r * bkm1.i + ak.i * bkm1.r; // , expr subst
                z__2.r = z__3.r - bk.r;
                z__2.i = z__3.i - bk.i; // , expr subst
                z_div(&z__1, &z__2, &denom);
                b[i__2].r = z__1.r;
                b[i__2].i = z__1.i; // , expr subst
                i__2 = k + 1 + j * b_dim1;
                z__3.r = akm1.r * bk.r - akm1.i * bk.i;
                z__3.i = akm1.r * bk.i + akm1.i * bk.r; // , expr subst
                z__2.r = z__3.r - bkm1.r;
                z__2.i = z__3.i - bkm1.i; // , expr subst
                z_div(&z__1, &z__2, &denom);
                b[i__2].r = z__1.r;
                b[i__2].i = z__1.i; // , expr subst
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
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + b_dim1], ldb);
            }
            /* Interchange rows K and IPIV(K). */
            kp = ipiv[k];
            if (kp != k)
            {
                zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
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
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], ldb, &a[k + 1 + k * a_dim1], &c__1, &c_b1, &b[k + b_dim1], ldb);
                i__1 = *n - k;
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemv_("Transpose", &i__1, nrhs, &z__1, &b[k + 1 + b_dim1], ldb, &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b1, &b[k - 1 + b_dim1], ldb);
            }
            /* Interchange rows K and -IPIV(K). */
            kp = -ipiv[k];
            if (kp != k)
            {
                zswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            k += -2;
        }
        goto L90;
L100:
        ;
    }
    return 0;
    /* End of ZSYTRS */
}
/* zsytrs_ */
