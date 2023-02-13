/* ../netlib/dpstrf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b22 = -1.;
static doublereal c_b24 = 1.;
/* > \brief \b DPSTRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DPSTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpstrf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpstrf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpstrf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* DOUBLE PRECISION TOL */
/* INTEGER INFO, LDA, N, RANK */
/* CHARACTER UPLO */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), WORK( 2*N ) */
/* INTEGER PIV( N ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPSTRF computes the Cholesky factorization with complete */
/* > pivoting of a real symmetric positive semidefinite matrix A. */
/* > */
/* > The factorization has the form */
/* > P**T * A * P = U**T * U , if UPLO = 'U', */
/* > P**T * A * P = L * L**T, if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular, and */
/* > P is stored as vector PIV. */
/* > */
/* > This algorithm does not attempt to check that A is positive */
/* > semidefinite. This version of the algorithm calls level 3 BLAS. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > symmetric matrix A is stored. */
/* > = 'U': Upper triangular */
/* > = 'L': Lower triangular */
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
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the symmetric matrix A. If UPLO = 'U', the leading */
/* > n by n upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading n by n lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, if INFO = 0, the factor U or L from the Cholesky */
/* > factorization as above. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] PIV */
/* > \verbatim */
/* > PIV is INTEGER array, dimension (N) */
/* > PIV is such that the nonzero entries are P( PIV(K), K ) = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* > RANK is INTEGER */
/* > The rank of A given by the number of steps the algorithm */
/* > completed. */
/* > \endverbatim */
/* > */
/* > \param[in] TOL */
/* > \verbatim */
/* > TOL is DOUBLE PRECISION */
/* > User defined tolerance. If TOL < 0, then N*U*MAX( A(K,K) ) */
/* > will be used. The algorithm terminates at the (K-1)st step */
/* > if the pivot <= TOL. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (2*N) */
/* > Work space. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > < 0: If INFO = -K, the K-th argument had an illegal value, */
/* > = 0: algorithm completed successfully, and */
/* > > 0: the matrix A is either rank deficient with computed rank */
/* > as returned in RANK, or is indefinite. See Section 7 of */
/* > LAPACK Working Note #161 for further information. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup doubleOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int dpstrf_(char *uplo, integer *n, doublereal *a, integer * lda, integer *piv, integer *rank, doublereal *tol, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k, jb, nb;
    doublereal ajj;
    integer pvt;
    extern /* Subroutine */
    int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    doublereal dtemp;
    integer itemp;
    extern /* Subroutine */
    int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal dstop;
    logical upper;
    extern /* Subroutine */
    int dsyrk_(char *, char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, integer *), dpstf2_(char *, integer *, doublereal *, integer *, integer *, integer *, doublereal *, doublereal *, integer *);
    extern doublereal dlamch_(char *);
    extern logical disnan_(doublereal *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *), dmaxloc_(doublereal *, integer *);
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
    /* Test the input parameters. */
    /* Parameter adjustments */
    --work;
    --piv;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
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
        xerbla_("DPSTRF", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Get block size */
    nb = ilaenv_(&c__1, "DPOTRF", uplo, n, &c_n1, &c_n1, &c_n1);
    if (nb <= 1 || nb >= *n)
    {
        /* Use unblocked code */
        dpstf2_(uplo, n, &a[a_dim1 + 1], lda, &piv[1], rank, tol, &work[1], info);
        goto L200;
    }
    else
    {
        /* Initialize PIV */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            piv[i__] = i__;
            /* L100: */
        }
        /* Compute stopping value */
        pvt = 1;
        ajj = a[pvt + pvt * a_dim1];
        i__1 = *n;
        for (i__ = 2;
                i__ <= i__1;
                ++i__)
        {
            if (a[i__ + i__ * a_dim1] > ajj)
            {
                pvt = i__;
                ajj = a[pvt + pvt * a_dim1];
            }
        }
        if (ajj == 0. || disnan_(&ajj))
        {
            *rank = 0;
            *info = 1;
            goto L200;
        }
        /* Compute stopping value if not supplied */
        if (*tol < 0.)
        {
            dstop = *n * dlamch_("Epsilon") * ajj;
        }
        else
        {
            dstop = *tol;
        }
        if (upper)
        {
            /* Compute the Cholesky factorization P**T * A * P = U**T * U */
            i__1 = *n;
            i__2 = nb;
            for (k = 1;
                    i__2 < 0 ? k >= i__1 : k <= i__1;
                    k += i__2)
            {
                /* Account for last block not being NB wide */
                /* Computing MIN */
                i__3 = nb;
                i__4 = *n - k + 1; // , expr subst
                jb = min(i__3,i__4);
                /* Set relevant part of first half of WORK to zero, */
                /* holds dot products */
                i__3 = *n;
                for (i__ = k;
                        i__ <= i__3;
                        ++i__)
                {
                    work[i__] = 0.;
                    /* L110: */
                }
                i__3 = k + jb - 1;
                for (j = k;
                        j <= i__3;
                        ++j)
                {
                    /* Find pivot, test for exit, else swap rows and columns */
                    /* Update dot products, compute possible pivots which are */
                    /* stored in the second half of WORK */
                    i__4 = *n;
                    for (i__ = j;
                            i__ <= i__4;
                            ++i__)
                    {
                        if (j > k)
                        {
                            /* Computing 2nd power */
                            d__1 = a[j - 1 + i__ * a_dim1];
                            work[i__] += d__1 * d__1;
                        }
                        work[*n + i__] = a[i__ + i__ * a_dim1] - work[i__];
                        /* L120: */
                    }
                    if (j > 1)
                    {
                        i__4 = *n - j + 1;
                        itemp = dmaxloc_(&work[*n + j], &i__4);
                        pvt = itemp + j - 1;
                        ajj = work[*n + pvt];
                        if (ajj <= dstop || disnan_(&ajj))
                        {
                            a[j + j * a_dim1] = ajj;
                            goto L190;
                        }
                    }
                    if (j != pvt)
                    {
                        /* Pivot OK, so can now swap pivot rows and columns */
                        a[pvt + pvt * a_dim1] = a[j + j * a_dim1];
                        i__4 = j - 1;
                        dswap_(&i__4, &a[j * a_dim1 + 1], &c__1, &a[pvt * a_dim1 + 1], &c__1);
                        if (pvt < *n)
                        {
                            i__4 = *n - pvt;
                            dswap_(&i__4, &a[j + (pvt + 1) * a_dim1], lda, &a[ pvt + (pvt + 1) * a_dim1], lda);
                        }
                        i__4 = pvt - j - 1;
                        dswap_(&i__4, &a[j + (j + 1) * a_dim1], lda, &a[j + 1 + pvt * a_dim1], &c__1);
                        /* Swap dot products and PIV */
                        dtemp = work[j];
                        work[j] = work[pvt];
                        work[pvt] = dtemp;
                        itemp = piv[pvt];
                        piv[pvt] = piv[j];
                        piv[j] = itemp;
                    }
                    ajj = sqrt(ajj);
                    a[j + j * a_dim1] = ajj;
                    /* Compute elements J+1:N of row J. */
                    if (j < *n)
                    {
                        i__4 = j - k;
                        i__5 = *n - j;
                        dgemv_("Trans", &i__4, &i__5, &c_b22, &a[k + (j + 1) * a_dim1], lda, &a[k + j * a_dim1], &c__1, & c_b24, &a[j + (j + 1) * a_dim1], lda);
                        i__4 = *n - j;
                        d__1 = 1. / ajj;
                        dscal_(&i__4, &d__1, &a[j + (j + 1) * a_dim1], lda);
                    }
                    /* L130: */
                }
                /* Update trailing matrix, J already incremented */
                if (k + jb <= *n)
                {
                    i__3 = *n - j + 1;
                    dsyrk_("Upper", "Trans", &i__3, &jb, &c_b22, &a[k + j * a_dim1], lda, &c_b24, &a[j + j * a_dim1], lda);
                }
                /* L140: */
            }
        }
        else
        {
            /* Compute the Cholesky factorization P**T * A * P = L * L**T */
            i__2 = *n;
            i__1 = nb;
            for (k = 1;
                    i__1 < 0 ? k >= i__2 : k <= i__2;
                    k += i__1)
            {
                /* Account for last block not being NB wide */
                /* Computing MIN */
                i__3 = nb;
                i__4 = *n - k + 1; // , expr subst
                jb = min(i__3,i__4);
                /* Set relevant part of first half of WORK to zero, */
                /* holds dot products */
                i__3 = *n;
                for (i__ = k;
                        i__ <= i__3;
                        ++i__)
                {
                    work[i__] = 0.;
                    /* L150: */
                }
                i__3 = k + jb - 1;
                for (j = k;
                        j <= i__3;
                        ++j)
                {
                    /* Find pivot, test for exit, else swap rows and columns */
                    /* Update dot products, compute possible pivots which are */
                    /* stored in the second half of WORK */
                    i__4 = *n;
                    for (i__ = j;
                            i__ <= i__4;
                            ++i__)
                    {
                        if (j > k)
                        {
                            /* Computing 2nd power */
                            d__1 = a[i__ + (j - 1) * a_dim1];
                            work[i__] += d__1 * d__1;
                        }
                        work[*n + i__] = a[i__ + i__ * a_dim1] - work[i__];
                        /* L160: */
                    }
                    if (j > 1)
                    {
                        i__4 = *n - j + 1;
                        itemp = dmaxloc_(&work[*n + j], &i__4);
                        pvt = itemp + j - 1;
                        ajj = work[*n + pvt];
                        if (ajj <= dstop || disnan_(&ajj))
                        {
                            a[j + j * a_dim1] = ajj;
                            goto L190;
                        }
                    }
                    if (j != pvt)
                    {
                        /* Pivot OK, so can now swap pivot rows and columns */
                        a[pvt + pvt * a_dim1] = a[j + j * a_dim1];
                        i__4 = j - 1;
                        dswap_(&i__4, &a[j + a_dim1], lda, &a[pvt + a_dim1], lda);
                        if (pvt < *n)
                        {
                            i__4 = *n - pvt;
                            dswap_(&i__4, &a[pvt + 1 + j * a_dim1], &c__1, &a[ pvt + 1 + pvt * a_dim1], &c__1);
                        }
                        i__4 = pvt - j - 1;
                        dswap_(&i__4, &a[j + 1 + j * a_dim1], &c__1, &a[pvt + (j + 1) * a_dim1], lda);
                        /* Swap dot products and PIV */
                        dtemp = work[j];
                        work[j] = work[pvt];
                        work[pvt] = dtemp;
                        itemp = piv[pvt];
                        piv[pvt] = piv[j];
                        piv[j] = itemp;
                    }
                    ajj = sqrt(ajj);
                    a[j + j * a_dim1] = ajj;
                    /* Compute elements J+1:N of column J. */
                    if (j < *n)
                    {
                        i__4 = *n - j;
                        i__5 = j - k;
                        dgemv_("No Trans", &i__4, &i__5, &c_b22, &a[j + 1 + k * a_dim1], lda, &a[j + k * a_dim1], lda, & c_b24, &a[j + 1 + j * a_dim1], &c__1);
                        i__4 = *n - j;
                        d__1 = 1. / ajj;
                        dscal_(&i__4, &d__1, &a[j + 1 + j * a_dim1], &c__1);
                    }
                    /* L170: */
                }
                /* Update trailing matrix, J already incremented */
                if (k + jb <= *n)
                {
                    i__3 = *n - j + 1;
                    dsyrk_("Lower", "No Trans", &i__3, &jb, &c_b22, &a[j + k * a_dim1], lda, &c_b24, &a[j + j * a_dim1], lda);
                }
                /* L180: */
            }
        }
    }
    /* Ran to completion, A has full rank */
    *rank = *n;
    goto L200;
L190: /* Rank is the number of steps completed. Set INFO = 1 to signal */
    /* that the factorization cannot be used to solve a system. */
    *rank = j - 1;
    *info = 1;
L200:
    return 0;
    /* End of DPSTRF */
}
/* dpstrf_ */
