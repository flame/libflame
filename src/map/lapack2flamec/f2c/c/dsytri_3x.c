/* ../netlib/v3.9.0/dsytri_3x.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b10 = 1.;
static doublereal c_b14 = 0.;
/* > \brief \b DSYTRI_3X */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSYTRI_3X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsytri_ 3x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsytri_ 3x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsytri_ 3x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSYTRI_3X( UPLO, N, A, LDA, E, IPIV, WORK, NB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N, NB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* DOUBLE PRECISION A( LDA, * ), E( * ), WORK( N+NB+1, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > DSYTRI_3X computes the inverse of a real symmetric indefinite */
/* > matrix A using the factorization computed by DSYTRF_RK or DSYTRF_BK: */
/* > */
/* > A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T), */
/* > */
/* > where U (or L) is unit upper (or lower) triangular matrix, */
/* > U**T (or L**T) is the transpose of U (or L), P is a permutation */
/* > matrix, P**T is the transpose of P, and D is symmetric and block */
/* > diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This is the blocked version of the algorithm, calling Level 3 BLAS. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are */
/* > stored as an upper or lower triangular matrix. */
/* > = 'U': Upper triangle of A is stored;
*/
/* > = 'L': Lower triangle of A is stored. */
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
/* > On entry, diagonal of the block diagonal matrix D and */
/* > factors U or L as computed by DSYTRF_RK and DSYTRF_BK: */
/* > a) ONLY diagonal elements of the symmetric block diagonal */
/* > matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
*/
/* > (superdiagonal (or subdiagonal) elements of D */
/* > should be provided on entry in array E), and */
/* > b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* > If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > */
/* > On exit, if INFO = 0, the symmetric inverse of the original */
/* > matrix. */
/* > If UPLO = 'U': the upper triangular part of the inverse */
/* > is formed and the part of A below the diagonal is not */
/* > referenced;
*/
/* > If UPLO = 'L': the lower triangular part of the inverse */
/* > is formed and the part of A above the diagonal is not */
/* > referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (N) */
/* > On entry, contains the superdiagonal (or subdiagonal) */
/* > elements of the symmetric block diagonal matrix D */
/* > with 1-by-1 or 2-by-2 diagonal blocks, where */
/* > If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) not referenced;
*/
/* > If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) not referenced. */
/* > */
/* > NOTE: For 1-by-1 diagonal block D(k), where */
/* > 1 <= k <= N, the element E(k) is not referenced in both */
/* > UPLO = 'U' or UPLO = 'L' cases. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by DSYTRF_RK or DSYTRF_BK. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (N+NB+1,NB+3). */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > Block size. */
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
/* > \date June 2017 */
/* > \ingroup doubleSYcomputational */
/* > \par Contributors: */
/* ================== */
/* > \verbatim */
/* > */
/* > June 2017, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int dsytri_3x_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *e, integer *ipiv, doublereal *work, integer *nb, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dsytri_3x inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS ", nb %" FLA_IS "",*uplo, *n, *lda, *nb);
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3;
    /* Local variables */
    extern /* Subroutine */
    int dsyswapr_(char *, integer *, doublereal *, integer *, integer *, integer *);
    doublereal d__;
    integer i__, j, k;
    doublereal t, ak;
    integer u11, ip, nnb, cut;
    doublereal akp1;
    integer invd;
    doublereal akkp1;
    extern /* Subroutine */
    int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    logical upper;
    doublereal u01_i_j__, u11_i_j__;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    integer icount;
    extern /* Subroutine */
    int dtrtri_(char *, char *, integer *, doublereal *, integer *, integer *);
    doublereal u01_ip1_j__, u11_ip1_j__;
    /* -- LAPACK computational routine (version 3.7.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2017 */
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
    --e;
    --ipiv;
    work_dim1 = *n + *nb + 1;
    work_offset = 1 + work_dim1;
    work -= work_offset;
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
    /* Quick return if possible */
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DSYTRI_3X", &i__1);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Workspace got Non-diag elements of D */
    i__1 = *n;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        work[k + work_dim1] = e[k];
    }
    /* Check that the diagonal matrix D is nonsingular. */
    if (upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for (*info = *n;
                *info >= 1;
                --(*info))
        {
            if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.)
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
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
            if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.)
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
        }
    }
    *info = 0;
    /* Splitting Workspace */
    /* U01 is a block ( N, NB+1 ) */
    /* The first element of U01 is in WORK( 1, 1 ) */
    /* U11 is a block ( NB+1, NB+1 ) */
    /* The first element of U11 is in WORK( N+1, 1 ) */
    u11 = *n;
    /* INVD is a block ( N, 2 ) */
    /* The first element of INVD is in WORK( 1, INVD ) */
    invd = *nb + 2;
    if (upper)
    {
        /* Begin Upper */
        /* invA = P * inv(U**T) * inv(D) * inv(U) * P**T. */
        dtrtri_(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D) * inv(U) */
        k = 1;
        while(k <= *n)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                work[k + invd * work_dim1] = 1. / a[k + k * a_dim1];
                work[k + (invd + 1) * work_dim1] = 0.;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                t = work[k + 1 + work_dim1];
                ak = a[k + k * a_dim1] / t;
                akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
                akkp1 = work[k + 1 + work_dim1] / t;
                d__ = t * (ak * akp1 - 1.);
                work[k + invd * work_dim1] = akp1 / d__;
                work[k + 1 + (invd + 1) * work_dim1] = ak / d__;
                work[k + (invd + 1) * work_dim1] = -akkp1 / d__;
                work[k + 1 + invd * work_dim1] = work[k + (invd + 1) * work_dim1];
                ++k;
            }
            ++k;
        }
        /* inv(U**T) = (inv(U))**T */
        /* inv(U**T) * inv(D) * inv(U) */
        cut = *n;
        while(cut > 0)
        {
            nnb = *nb;
            if (cut <= nnb)
            {
                nnb = cut;
            }
            else
            {
                icount = 0;
                /* count negative elements, */
                i__1 = cut;
                for (i__ = cut + 1 - nnb;
                        i__ <= i__1;
                        ++i__)
                {
                    if (ipiv[i__] < 0)
                    {
                        ++icount;
                    }
                }
                /* need a even number for a clear cut */
                if (icount % 2 == 1)
                {
                    ++nnb;
                }
            }
            cut -= nnb;
            /* U01 Block */
            i__1 = cut;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = nnb;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    work[i__ + j * work_dim1] = a[i__ + (cut + j) * a_dim1];
                }
            }
            /* U11 Block */
            i__1 = nnb;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                work[u11 + i__ + i__ * work_dim1] = 1.;
                i__2 = i__ - 1;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    work[u11 + i__ + j * work_dim1] = 0.;
                }
                i__2 = nnb;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    work[u11 + i__ + j * work_dim1] = a[cut + i__ + (cut + j) * a_dim1];
                }
            }
            /* invD * U01 */
            i__ = 1;
            while(i__ <= cut)
            {
                if (ipiv[i__] > 0)
                {
                    i__1 = nnb;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        work[i__ + j * work_dim1] = work[i__ + invd * work_dim1] * work[i__ + j * work_dim1];
                    }
                }
                else
                {
                    i__1 = nnb;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        u01_i_j__ = work[i__ + j * work_dim1];
                        u01_ip1_j__ = work[i__ + 1 + j * work_dim1];
                        work[i__ + j * work_dim1] = work[i__ + invd * work_dim1] * u01_i_j__ + work[i__ + (invd + 1) * work_dim1] * u01_ip1_j__;
                        work[i__ + 1 + j * work_dim1] = work[i__ + 1 + invd * work_dim1] * u01_i_j__ + work[i__ + 1 + (invd + 1) * work_dim1] * u01_ip1_j__;
                    }
                    ++i__;
                }
                ++i__;
            }
            /* invD1 * U11 */
            i__ = 1;
            while(i__ <= nnb)
            {
                if (ipiv[cut + i__] > 0)
                {
                    i__1 = nnb;
                    for (j = i__;
                            j <= i__1;
                            ++j)
                    {
                        work[u11 + i__ + j * work_dim1] = work[cut + i__ + invd * work_dim1] * work[u11 + i__ + j * work_dim1];
                    }
                }
                else
                {
                    i__1 = nnb;
                    for (j = i__;
                            j <= i__1;
                            ++j)
                    {
                        u11_i_j__ = work[u11 + i__ + j * work_dim1];
                        u11_ip1_j__ = work[u11 + i__ + 1 + j * work_dim1];
                        work[u11 + i__ + j * work_dim1] = work[cut + i__ + invd * work_dim1] * work[u11 + i__ + j * work_dim1] + work[cut + i__ + (invd + 1) * work_dim1] * work[u11 + i__ + 1 + j * work_dim1];
                        work[u11 + i__ + 1 + j * work_dim1] = work[cut + i__ + 1 + invd * work_dim1] * u11_i_j__ + work[ cut + i__ + 1 + (invd + 1) * work_dim1] * u11_ip1_j__;
                    }
                    ++i__;
                }
                ++i__;
            }
            /* U11**T * invD1 * U11 -> U11 */
            i__1 = *n + *nb + 1;
            dtrmm_("L", "U", "T", "U", &nnb, &nnb, &c_b10, &a[cut + 1 + (cut + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1);
            i__1 = nnb;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = nnb;
                for (j = i__;
                        j <= i__2;
                        ++j)
                {
                    a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + j * work_dim1];
                }
            }
            /* U01**T * invD * U01 -> A( CUT+I, CUT+J ) */
            i__1 = *n + *nb + 1;
            i__2 = *n + *nb + 1;
            dgemm_("T", "N", &nnb, &nnb, &cut, &c_b10, &a[(cut + 1) * a_dim1 + 1], lda, &work[work_offset], &i__1, &c_b14, &work[u11 + 1 + work_dim1], &i__2);
            /* U11 = U11**T * invD1 * U11 + U01**T * invD * U01 */
            i__1 = nnb;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = nnb;
                for (j = i__;
                        j <= i__2;
                        ++j)
                {
                    a[cut + i__ + (cut + j) * a_dim1] += work[u11 + i__ + j * work_dim1];
                }
            }
            /* U01 = U00**T * invD0 * U01 */
            i__1 = *n + *nb + 1;
            dtrmm_("L", uplo, "T", "U", &cut, &nnb, &c_b10, &a[a_offset], lda, &work[work_offset], &i__1);
            /* Update U01 */
            i__1 = cut;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = nnb;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    a[i__ + (cut + j) * a_dim1] = work[i__ + j * work_dim1];
                }
            }
            /* Next Block */
        }
        /* Apply PERMUTATIONS P and P**T: */
        /* P * inv(U**T) * inv(D) * inv(U) * P**T. */
        /* Interchange rows and columns I and IPIV(I) in reverse order */
        /* from the formation order of IPIV vector for Upper case. */
        /* ( We can use a loop over IPIV with increment 1, */
        /* since the ABS value of IPIV(I) represents the row (column) */
        /* index of the interchange with row (column) i in both 1x1 */
        /* and 2x2 pivot cases, i.e. we don't need separate code branches */
        /* for 1x1 and 2x2 pivot cases ) */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            ip = (i__2 = ipiv[i__], f2c_dabs(i__2));
            if (ip != i__)
            {
                if (i__ < ip)
                {
                    dsyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    dsyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
        }
    }
    else
    {
        /* Begin Lower */
        /* inv A = P * inv(L**T) * inv(D) * inv(L) * P**T. */
        dtrtri_(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D) * inv(L) */
        k = *n;
        while(k >= 1)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                work[k + invd * work_dim1] = 1. / a[k + k * a_dim1];
                work[k + (invd + 1) * work_dim1] = 0.;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                t = work[k - 1 + work_dim1];
                ak = a[k - 1 + (k - 1) * a_dim1] / t;
                akp1 = a[k + k * a_dim1] / t;
                akkp1 = work[k - 1 + work_dim1] / t;
                d__ = t * (ak * akp1 - 1.);
                work[k - 1 + invd * work_dim1] = akp1 / d__;
                work[k + invd * work_dim1] = ak / d__;
                work[k + (invd + 1) * work_dim1] = -akkp1 / d__;
                work[k - 1 + (invd + 1) * work_dim1] = work[k + (invd + 1) * work_dim1];
                --k;
            }
            --k;
        }
        /* inv(L**T) = (inv(L))**T */
        /* inv(L**T) * inv(D) * inv(L) */
        cut = 0;
        while(cut < *n)
        {
            nnb = *nb;
            if (cut + nnb > *n)
            {
                nnb = *n - cut;
            }
            else
            {
                icount = 0;
                /* count negative elements, */
                i__1 = cut + nnb;
                for (i__ = cut + 1;
                        i__ <= i__1;
                        ++i__)
                {
                    if (ipiv[i__] < 0)
                    {
                        ++icount;
                    }
                }
                /* need a even number for a clear cut */
                if (icount % 2 == 1)
                {
                    ++nnb;
                }
            }
            /* L21 Block */
            i__1 = *n - cut - nnb;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = nnb;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    work[i__ + j * work_dim1] = a[cut + nnb + i__ + (cut + j) * a_dim1];
                }
            }
            /* L11 Block */
            i__1 = nnb;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                work[u11 + i__ + i__ * work_dim1] = 1.;
                i__2 = nnb;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    work[u11 + i__ + j * work_dim1] = 0.;
                }
                i__2 = i__ - 1;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    work[u11 + i__ + j * work_dim1] = a[cut + i__ + (cut + j) * a_dim1];
                }
            }
            /* invD*L21 */
            i__ = *n - cut - nnb;
            while(i__ >= 1)
            {
                if (ipiv[cut + nnb + i__] > 0)
                {
                    i__1 = nnb;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        work[i__ + j * work_dim1] = work[cut + nnb + i__ + invd * work_dim1] * work[i__ + j * work_dim1];
                    }
                }
                else
                {
                    i__1 = nnb;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        u01_i_j__ = work[i__ + j * work_dim1];
                        u01_ip1_j__ = work[i__ - 1 + j * work_dim1];
                        work[i__ + j * work_dim1] = work[cut + nnb + i__ + invd * work_dim1] * u01_i_j__ + work[cut + nnb + i__ + (invd + 1) * work_dim1] * u01_ip1_j__;
                        work[i__ - 1 + j * work_dim1] = work[cut + nnb + i__ - 1 + (invd + 1) * work_dim1] * u01_i_j__ + work[cut + nnb + i__ - 1 + invd * work_dim1] * u01_ip1_j__;
                    }
                    --i__;
                }
                --i__;
            }
            /* invD1*L11 */
            i__ = nnb;
            while(i__ >= 1)
            {
                if (ipiv[cut + i__] > 0)
                {
                    i__1 = nnb;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        work[u11 + i__ + j * work_dim1] = work[cut + i__ + invd * work_dim1] * work[u11 + i__ + j * work_dim1];
                    }
                }
                else
                {
                    i__1 = nnb;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        u11_i_j__ = work[u11 + i__ + j * work_dim1];
                        u11_ip1_j__ = work[u11 + i__ - 1 + j * work_dim1];
                        work[u11 + i__ + j * work_dim1] = work[cut + i__ + invd * work_dim1] * work[u11 + i__ + j * work_dim1] + work[cut + i__ + (invd + 1) * work_dim1] * u11_ip1_j__;
                        work[u11 + i__ - 1 + j * work_dim1] = work[cut + i__ - 1 + (invd + 1) * work_dim1] * u11_i_j__ + work[cut + i__ - 1 + invd * work_dim1] * u11_ip1_j__;
                    }
                    --i__;
                }
                --i__;
            }
            /* L11**T * invD1 * L11 -> L11 */
            i__1 = *n + *nb + 1;
            dtrmm_("L", uplo, "T", "U", &nnb, &nnb, &c_b10, &a[cut + 1 + (cut + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1);
            i__1 = nnb;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = i__;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + j * work_dim1];
                }
            }
            if (cut + nnb < *n)
            {
                /* L21**T * invD2*L21 -> A( CUT+I, CUT+J ) */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                i__3 = *n + *nb + 1;
                dgemm_("T", "N", &nnb, &nnb, &i__1, &c_b10, &a[cut + nnb + 1 + (cut + 1) * a_dim1], lda, &work[work_offset], &i__2, &c_b14, &work[u11 + 1 + work_dim1], &i__3);
                /* L11 = L11**T * invD1 * L11 + U01**T * invD * U01 */
                i__1 = nnb;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = i__;
                    for (j = 1;
                            j <= i__2;
                            ++j)
                    {
                        a[cut + i__ + (cut + j) * a_dim1] += work[u11 + i__ + j * work_dim1];
                    }
                }
                /* L01 = L22**T * invD2 * L21 */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                dtrmm_("L", uplo, "T", "U", &i__1, &nnb, &c_b10, &a[cut + nnb + 1 + (cut + nnb + 1) * a_dim1], lda, &work[ work_offset], &i__2);
                /* Update L21 */
                i__1 = *n - cut - nnb;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = nnb;
                    for (j = 1;
                            j <= i__2;
                            ++j)
                    {
                        a[cut + nnb + i__ + (cut + j) * a_dim1] = work[i__ + j * work_dim1];
                    }
                }
            }
            else
            {
                /* L11 = L11**T * invD1 * L11 */
                i__1 = nnb;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = i__;
                    for (j = 1;
                            j <= i__2;
                            ++j)
                    {
                        a[cut + i__ + (cut + j) * a_dim1] = work[u11 + i__ + j * work_dim1];
                    }
                }
            }
            /* Next Block */
            cut += nnb;
        }
        /* Apply PERMUTATIONS P and P**T: */
        /* P * inv(L**T) * inv(D) * inv(L) * P**T. */
        /* Interchange rows and columns I and IPIV(I) in reverse order */
        /* from the formation order of IPIV vector for Lower case. */
        /* ( We can use a loop over IPIV with increment -1, */
        /* since the ABS value of IPIV(I) represents the row (column) */
        /* index of the interchange with row (column) i in both 1x1 */
        /* and 2x2 pivot cases, i.e. we don't need separate code branches */
        /* for 1x1 and 2x2 pivot cases ) */
        for (i__ = *n;
                i__ >= 1;
                --i__)
        {
            ip = (i__1 = ipiv[i__], f2c_dabs(i__1));
            if (ip != i__)
            {
                if (i__ < ip)
                {
                    dsyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    dsyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of DSYTRI_3X */
}
/* dsytri_3x__ */

