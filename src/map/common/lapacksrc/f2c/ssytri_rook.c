/* ../netlib/ssytri_rook.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b11 = -1.f;
static real c_b13 = 0.f;
/* > \brief \b SSYTRI_ROOK */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSYTRI_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri_ rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri_ rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri_ rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSYTRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYTRI_ROOK computes the inverse of a real symmetric */
/* > matrix A using the factorization A = U*D*U**T or A = L*D*L**T */
/* > computed by SSYTRF_ROOK. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the block diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by SSYTRF_ROOK. */
/* > */
/* > On exit, if INFO = 0, the (symmetric) inverse of the original */
/* > matrix. If UPLO = 'U', the upper triangular part of the */
/* > inverse is formed and the part of A below the diagonal is not */
/* > referenced;
if UPLO = 'L' the lower triangular part of the */
/* > inverse is formed and the part of A above the diagonal is */
/* > not referenced. */
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
/* > as determined by SSYTRF_ROOK. */
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
/* > > 0: if INFO = i, D(i,i) = 0;
the matrix is singular and its */
/* > inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date April 2012 */
/* > \ingroup realSYcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > April 2012, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* > School of Mathematics, */
/* > University of Manchester */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int ssytri_rook_(char *uplo, integer *n, real *a, integer * lda, integer *ipiv, real *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    real r__1;
    /* Local variables */
    real d__;
    integer k;
    real t, ak;
    integer kp;
    real akp1, temp;
    extern real sdot_(integer *, real *, integer *, real *, integer *);
    real akkp1;
    extern logical lsame_(char *, char *);
    integer kstep;
    logical upper;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *), sswap_(integer *, real *, integer *, real *, integer * ), ssymv_(char *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *), xerbla_(char *, integer *);
    /* -- LAPACK computational routine (version 3.4.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* April 2012 */
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
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
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
    else if (*lda < max(1,*n))
    {
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSYTRI_ROOK", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Check that the diagonal matrix D is nonsingular. */
    if (upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for (*info = *n;
                *info >= 1;
                --(*info))
        {
            if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.f)
            {
                return 0;
            }
            /* L10: */
        }
    }
    else
    {
        /* Lower triangular storage: examine D from top to bottom. */
        i__1 = *n;
        for (*info = 1;
                *info <= i__1;
                ++(*info))
        {
            if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.f)
            {
                return 0;
            }
            /* L20: */
        }
    }
    *info = 0;
    if (upper)
    {
        /* Compute inv(A) from the factorization A = U*D*U**T. */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        k = 1;
L30: /* If K > N, exit from loop. */
        if (k > *n)
        {
            goto L40;
        }
        if (ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Invert the diagonal block. */
            a[k + k * a_dim1] = 1.f / a[k + k * a_dim1];
            /* Compute column K of the inverse. */
            if (k > 1)
            {
                i__1 = k - 1;
                scopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                ssymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], & c__1, &c_b13, &a[k * a_dim1 + 1], &c__1);
                i__1 = k - 1;
                a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
            }
            kstep = 1;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Invert the diagonal block. */
            t = (r__1 = a[k + (k + 1) * a_dim1], f2c_abs(r__1));
            ak = a[k + k * a_dim1] / t;
            akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
            akkp1 = a[k + (k + 1) * a_dim1] / t;
            d__ = t * (ak * akp1 - 1.f);
            a[k + k * a_dim1] = akp1 / d__;
            a[k + 1 + (k + 1) * a_dim1] = ak / d__;
            a[k + (k + 1) * a_dim1] = -akkp1 / d__;
            /* Compute columns K and K+1 of the inverse. */
            if (k > 1)
            {
                i__1 = k - 1;
                scopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                ssymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], & c__1, &c_b13, &a[k * a_dim1 + 1], &c__1);
                i__1 = k - 1;
                a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
                i__1 = k - 1;
                a[k + (k + 1) * a_dim1] -= sdot_(&i__1, &a[k * a_dim1 + 1], & c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
                i__1 = k - 1;
                scopy_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], & c__1);
                i__1 = k - 1;
                ssymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], & c__1, &c_b13, &a[(k + 1) * a_dim1 + 1], &c__1);
                i__1 = k - 1;
                a[k + 1 + (k + 1) * a_dim1] -= sdot_(&i__1, &work[1], &c__1, & a[(k + 1) * a_dim1 + 1], &c__1);
            }
            kstep = 2;
        }
        if (kstep == 1)
        {
            /* Interchange rows and columns K and IPIV(K) in the leading */
            /* submatrix A(1:k+1,1:k+1) */
            kp = ipiv[k];
            if (kp != k)
            {
                if (kp > 1)
                {
                    i__1 = kp - 1;
                    sswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
                }
                i__1 = k - kp - 1;
                sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1], lda);
                temp = a[k + k * a_dim1];
                a[k + k * a_dim1] = a[kp + kp * a_dim1];
                a[kp + kp * a_dim1] = temp;
            }
        }
        else
        {
            /* Interchange rows and columns K and K+1 with -IPIV(K) and */
            /* -IPIV(K+1)in the leading submatrix A(1:k+1,1:k+1) */
            kp = -ipiv[k];
            if (kp != k)
            {
                if (kp > 1)
                {
                    i__1 = kp - 1;
                    sswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
                }
                i__1 = k - kp - 1;
                sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1], lda);
                temp = a[k + k * a_dim1];
                a[k + k * a_dim1] = a[kp + kp * a_dim1];
                a[kp + kp * a_dim1] = temp;
                temp = a[k + (k + 1) * a_dim1];
                a[k + (k + 1) * a_dim1] = a[kp + (k + 1) * a_dim1];
                a[kp + (k + 1) * a_dim1] = temp;
            }
            ++k;
            kp = -ipiv[k];
            if (kp != k)
            {
                if (kp > 1)
                {
                    i__1 = kp - 1;
                    sswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
                }
                i__1 = k - kp - 1;
                sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1], lda);
                temp = a[k + k * a_dim1];
                a[k + k * a_dim1] = a[kp + kp * a_dim1];
                a[kp + kp * a_dim1] = temp;
            }
        }
        ++k;
        goto L30;
L40:
        ;
    }
    else
    {
        /* Compute inv(A) from the factorization A = L*D*L**T. */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        k = *n;
L50: /* If K < 1, exit from loop. */
        if (k < 1)
        {
            goto L60;
        }
        if (ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Invert the diagonal block. */
            a[k + k * a_dim1] = 1.f / a[k + k * a_dim1];
            /* Compute column K of the inverse. */
            if (k < *n)
            {
                i__1 = *n - k;
                scopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                ssymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1], &c__1, &c_b13, &a[k + 1 + k * a_dim1], & c__1);
                i__1 = *n - k;
                a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k + 1 + k * a_dim1], &c__1);
            }
            kstep = 1;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Invert the diagonal block. */
            t = (r__1 = a[k + (k - 1) * a_dim1], f2c_abs(r__1));
            ak = a[k - 1 + (k - 1) * a_dim1] / t;
            akp1 = a[k + k * a_dim1] / t;
            akkp1 = a[k + (k - 1) * a_dim1] / t;
            d__ = t * (ak * akp1 - 1.f);
            a[k - 1 + (k - 1) * a_dim1] = akp1 / d__;
            a[k + k * a_dim1] = ak / d__;
            a[k + (k - 1) * a_dim1] = -akkp1 / d__;
            /* Compute columns K-1 and K of the inverse. */
            if (k < *n)
            {
                i__1 = *n - k;
                scopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                ssymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1], &c__1, &c_b13, &a[k + 1 + k * a_dim1], & c__1);
                i__1 = *n - k;
                a[k + k * a_dim1] -= sdot_(&i__1, &work[1], &c__1, &a[k + 1 + k * a_dim1], &c__1);
                i__1 = *n - k;
                a[k + (k - 1) * a_dim1] -= sdot_(&i__1, &a[k + 1 + k * a_dim1] , &c__1, &a[k + 1 + (k - 1) * a_dim1], &c__1);
                i__1 = *n - k;
                scopy_(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], & c__1);
                i__1 = *n - k;
                ssymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1], &c__1, &c_b13, &a[k + 1 + (k - 1) * a_dim1] , &c__1);
                i__1 = *n - k;
                a[k - 1 + (k - 1) * a_dim1] -= sdot_(&i__1, &work[1], &c__1, & a[k + 1 + (k - 1) * a_dim1], &c__1);
            }
            kstep = 2;
        }
        if (kstep == 1)
        {
            /* Interchange rows and columns K and IPIV(K) in the trailing */
            /* submatrix A(k-1:n,k-1:n) */
            kp = ipiv[k];
            if (kp != k)
            {
                if (kp < *n)
                {
                    i__1 = *n - kp;
                    sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1], &c__1);
                }
                i__1 = kp - k - 1;
                sswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) * a_dim1], lda);
                temp = a[k + k * a_dim1];
                a[k + k * a_dim1] = a[kp + kp * a_dim1];
                a[kp + kp * a_dim1] = temp;
            }
        }
        else
        {
            /* Interchange rows and columns K and K-1 with -IPIV(K) and */
            /* -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n) */
            kp = -ipiv[k];
            if (kp != k)
            {
                if (kp < *n)
                {
                    i__1 = *n - kp;
                    sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1], &c__1);
                }
                i__1 = kp - k - 1;
                sswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) * a_dim1], lda);
                temp = a[k + k * a_dim1];
                a[k + k * a_dim1] = a[kp + kp * a_dim1];
                a[kp + kp * a_dim1] = temp;
                temp = a[k + (k - 1) * a_dim1];
                a[k + (k - 1) * a_dim1] = a[kp + (k - 1) * a_dim1];
                a[kp + (k - 1) * a_dim1] = temp;
            }
            --k;
            kp = -ipiv[k];
            if (kp != k)
            {
                if (kp < *n)
                {
                    i__1 = *n - kp;
                    sswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1], &c__1);
                }
                i__1 = kp - k - 1;
                sswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) * a_dim1], lda);
                temp = a[k + k * a_dim1];
                a[k + k * a_dim1] = a[kp + kp * a_dim1];
                a[kp + kp * a_dim1] = temp;
            }
        }
        --k;
        goto L50;
L60:
        ;
    }
    return 0;
    /* End of SSYTRI_ROOK */
}
/* ssytri_rook__ */
