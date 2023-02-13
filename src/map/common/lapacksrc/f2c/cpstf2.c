/* ../netlib/cpstf2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    1.f,0.f
}
;
static integer c__1 = 1;
/* > \brief \b CPSTF2 computes the Cholesky factorization with complete pivoting of a real symmetric or comple x Hermitian positive semi-definite matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CPSTF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpstf2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpstf2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpstf2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* REAL TOL */
/* INTEGER INFO, LDA, N, RANK */
/* CHARACTER UPLO */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ) */
/* REAL WORK( 2*N ) */
/* INTEGER PIV( N ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPSTF2 computes the Cholesky factorization with complete */
/* > pivoting of a complex Hermitian positive semidefinite matrix A. */
/* > */
/* > The factorization has the form */
/* > P**T * A * P = U**H * U , if UPLO = 'U', */
/* > P**T * A * P = L * L**H, if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular, and */
/* > P is stored as vector PIV. */
/* > */
/* > This algorithm does not attempt to check that A is positive */
/* > semidefinite. This version of the algorithm calls level 2 BLAS. */
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
/* > A is COMPLEX array, dimension (LDA,N) */
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
/* > TOL is REAL */
/* > User defined tolerance. If TOL < 0, then N*U*MAX( A( K,K ) ) */
/* > will be used. The algorithm terminates at the (K-1)st step */
/* > if the pivot <= TOL. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (2*N) */
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
/* > \date September 2012 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int cpstf2_(char *uplo, integer *n, complex *a, integer *lda, integer *piv, integer *rank, real *tol, real *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;
    complex q__1, q__2;
    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j;
    real ajj;
    integer pvt;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int cgemv_(char *, integer *, integer *, complex * , complex *, integer *, complex *, integer *, complex *, complex * , integer *);
    complex ctemp;
    extern /* Subroutine */
    int cswap_(integer *, complex *, integer *, complex *, integer *);
    integer itemp;
    real stemp;
    logical upper;
    real sstop;
    extern /* Subroutine */
    int clacgv_(integer *, complex *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int csscal_(integer *, real *, complex *, integer *), xerbla_(char *, integer *);
    extern integer smaxloc_(real *, integer *);
    extern logical sisnan_(real *);
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* Test the input parameters */
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
        xerbla_("CPSTF2", &i__1);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
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
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__ + i__ * a_dim1;
        work[i__] = a[i__2].r;
        /* L110: */
    }
    pvt = smaxloc_(&work[1], n);
    i__1 = pvt + pvt * a_dim1;
    ajj = a[i__1].r;
    if (ajj == 0.f || sisnan_(&ajj))
    {
        *rank = 0;
        *info = 1;
        goto L200;
    }
    /* Compute stopping value if not supplied */
    if (*tol < 0.f)
    {
        sstop = *n * slamch_("Epsilon") * ajj;
    }
    else
    {
        sstop = *tol;
    }
    /* Set first half of WORK to zero, holds dot products */
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        work[i__] = 0.f;
        /* L120: */
    }
    if (upper)
    {
        /* Compute the Cholesky factorization P**T * A * P = U**H * U */
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Find pivot, test for exit, else swap rows and columns */
            /* Update dot products, compute possible pivots which are */
            /* stored in the second half of WORK */
            i__2 = *n;
            for (i__ = j;
                    i__ <= i__2;
                    ++i__)
            {
                if (j > 1)
                {
                    r_cnjg(&q__2, &a[j - 1 + i__ * a_dim1]);
                    i__3 = j - 1 + i__ * a_dim1;
                    q__1.r = q__2.r * a[i__3].r - q__2.i * a[i__3].i;
                    q__1.i = q__2.r * a[i__3].i + q__2.i * a[i__3].r; // , expr subst
                    work[i__] += q__1.r;
                }
                i__3 = i__ + i__ * a_dim1;
                work[*n + i__] = a[i__3].r - work[i__];
                /* L130: */
            }
            if (j > 1)
            {
                i__2 = *n - j + 1;
                itemp = smaxloc_(&work[*n + j], &i__2);
                pvt = itemp + j - 1;
                ajj = work[*n + pvt];
                if (ajj <= sstop || sisnan_(&ajj))
                {
                    i__2 = j + j * a_dim1;
                    a[i__2].r = ajj;
                    a[i__2].i = 0.f; // , expr subst
                    goto L190;
                }
            }
            if (j != pvt)
            {
                /* Pivot OK, so can now swap pivot rows and columns */
                i__2 = pvt + pvt * a_dim1;
                i__3 = j + j * a_dim1;
                a[i__2].r = a[i__3].r;
                a[i__2].i = a[i__3].i; // , expr subst
                i__2 = j - 1;
                cswap_(&i__2, &a[j * a_dim1 + 1], &c__1, &a[pvt * a_dim1 + 1], &c__1);
                if (pvt < *n)
                {
                    i__2 = *n - pvt;
                    cswap_(&i__2, &a[j + (pvt + 1) * a_dim1], lda, &a[pvt + ( pvt + 1) * a_dim1], lda);
                }
                i__2 = pvt - 1;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    r_cnjg(&q__1, &a[j + i__ * a_dim1]);
                    ctemp.r = q__1.r;
                    ctemp.i = q__1.i; // , expr subst
                    i__3 = j + i__ * a_dim1;
                    r_cnjg(&q__1, &a[i__ + pvt * a_dim1]);
                    a[i__3].r = q__1.r;
                    a[i__3].i = q__1.i; // , expr subst
                    i__3 = i__ + pvt * a_dim1;
                    a[i__3].r = ctemp.r;
                    a[i__3].i = ctemp.i; // , expr subst
                    /* L140: */
                }
                i__2 = j + pvt * a_dim1;
                r_cnjg(&q__1, &a[j + pvt * a_dim1]);
                a[i__2].r = q__1.r;
                a[i__2].i = q__1.i; // , expr subst
                /* Swap dot products and PIV */
                stemp = work[j];
                work[j] = work[pvt];
                work[pvt] = stemp;
                itemp = piv[pvt];
                piv[pvt] = piv[j];
                piv[j] = itemp;
            }
            ajj = sqrt(ajj);
            i__2 = j + j * a_dim1;
            a[i__2].r = ajj;
            a[i__2].i = 0.f; // , expr subst
            /* Compute elements J+1:N of row J */
            if (j < *n)
            {
                i__2 = j - 1;
                clacgv_(&i__2, &a[j * a_dim1 + 1], &c__1);
                i__2 = j - 1;
                i__3 = *n - j;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgemv_("Trans", &i__2, &i__3, &q__1, &a[(j + 1) * a_dim1 + 1], lda, &a[j * a_dim1 + 1], &c__1, &c_b1, &a[j + (j + 1) * a_dim1], lda);
                i__2 = j - 1;
                clacgv_(&i__2, &a[j * a_dim1 + 1], &c__1);
                i__2 = *n - j;
                r__1 = 1.f / ajj;
                csscal_(&i__2, &r__1, &a[j + (j + 1) * a_dim1], lda);
            }
            /* L150: */
        }
    }
    else
    {
        /* Compute the Cholesky factorization P**T * A * P = L * L**H */
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            /* Find pivot, test for exit, else swap rows and columns */
            /* Update dot products, compute possible pivots which are */
            /* stored in the second half of WORK */
            i__2 = *n;
            for (i__ = j;
                    i__ <= i__2;
                    ++i__)
            {
                if (j > 1)
                {
                    r_cnjg(&q__2, &a[i__ + (j - 1) * a_dim1]);
                    i__3 = i__ + (j - 1) * a_dim1;
                    q__1.r = q__2.r * a[i__3].r - q__2.i * a[i__3].i;
                    q__1.i = q__2.r * a[i__3].i + q__2.i * a[i__3].r; // , expr subst
                    work[i__] += q__1.r;
                }
                i__3 = i__ + i__ * a_dim1;
                work[*n + i__] = a[i__3].r - work[i__];
                /* L160: */
            }
            if (j > 1)
            {
                i__2 = *n - j + 1;
                itemp = smaxloc_(&work[*n + j], &i__2);
                pvt = itemp + j - 1;
                ajj = work[*n + pvt];
                if (ajj <= sstop || sisnan_(&ajj))
                {
                    i__2 = j + j * a_dim1;
                    a[i__2].r = ajj;
                    a[i__2].i = 0.f; // , expr subst
                    goto L190;
                }
            }
            if (j != pvt)
            {
                /* Pivot OK, so can now swap pivot rows and columns */
                i__2 = pvt + pvt * a_dim1;
                i__3 = j + j * a_dim1;
                a[i__2].r = a[i__3].r;
                a[i__2].i = a[i__3].i; // , expr subst
                i__2 = j - 1;
                cswap_(&i__2, &a[j + a_dim1], lda, &a[pvt + a_dim1], lda);
                if (pvt < *n)
                {
                    i__2 = *n - pvt;
                    cswap_(&i__2, &a[pvt + 1 + j * a_dim1], &c__1, &a[pvt + 1 + pvt * a_dim1], &c__1);
                }
                i__2 = pvt - 1;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    r_cnjg(&q__1, &a[i__ + j * a_dim1]);
                    ctemp.r = q__1.r;
                    ctemp.i = q__1.i; // , expr subst
                    i__3 = i__ + j * a_dim1;
                    r_cnjg(&q__1, &a[pvt + i__ * a_dim1]);
                    a[i__3].r = q__1.r;
                    a[i__3].i = q__1.i; // , expr subst
                    i__3 = pvt + i__ * a_dim1;
                    a[i__3].r = ctemp.r;
                    a[i__3].i = ctemp.i; // , expr subst
                    /* L170: */
                }
                i__2 = pvt + j * a_dim1;
                r_cnjg(&q__1, &a[pvt + j * a_dim1]);
                a[i__2].r = q__1.r;
                a[i__2].i = q__1.i; // , expr subst
                /* Swap dot products and PIV */
                stemp = work[j];
                work[j] = work[pvt];
                work[pvt] = stemp;
                itemp = piv[pvt];
                piv[pvt] = piv[j];
                piv[j] = itemp;
            }
            ajj = sqrt(ajj);
            i__2 = j + j * a_dim1;
            a[i__2].r = ajj;
            a[i__2].i = 0.f; // , expr subst
            /* Compute elements J+1:N of column J */
            if (j < *n)
            {
                i__2 = j - 1;
                clacgv_(&i__2, &a[j + a_dim1], lda);
                i__2 = *n - j;
                i__3 = j - 1;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgemv_("No Trans", &i__2, &i__3, &q__1, &a[j + 1 + a_dim1], lda, &a[j + a_dim1], lda, &c_b1, &a[j + 1 + j * a_dim1], &c__1);
                i__2 = j - 1;
                clacgv_(&i__2, &a[j + a_dim1], lda);
                i__2 = *n - j;
                r__1 = 1.f / ajj;
                csscal_(&i__2, &r__1, &a[j + 1 + j * a_dim1], &c__1);
            }
            /* L180: */
        }
    }
    /* Ran to completion, A has full rank */
    *rank = *n;
    goto L200;
L190: /* Rank is number of steps completed. Set INFO = 1 to signal */
    /* that the factorization cannot be used to solve a system. */
    *rank = j - 1;
    *info = 1;
L200:
    return 0;
    /* End of CPSTF2 */
}
/* cpstf2_ */
