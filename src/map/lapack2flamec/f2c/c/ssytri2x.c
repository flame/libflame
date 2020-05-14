/* ../netlib/ssytri2x.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b11 = 1.f;
static real c_b15 = 0.f;
/* > \brief \b SSYTRI2X */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSYTRI2X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri2 x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri2 x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri2 x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N, NB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL A( LDA, * ), WORK( N+NB+1,* ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYTRI2X computes the inverse of a real symmetric indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > SSYTRF. */
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
/* > On entry, the NNB diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by SSYTRF. */
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
/* > Details of the interchanges and the NNB structure of D */
/* > as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (N+NNB+1,NNB+3) */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > Block size */
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
/* > \date November 2011 */
/* > \ingroup realSYcomputational */
/* ===================================================================== */
/* Subroutine */
int ssytri2x_(char *uplo, integer *n, real *a, integer *lda, integer *ipiv, real *work, integer *nb, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3;
    /* Local variables */
    real d__;
    integer i__, j, k;
    extern /* Subroutine */
    int ssyswapr_(char *, integer *, real *, integer * , integer *, integer *);
    real t, ak;
    integer u11, ip, nnb, cut;
    real akp1;
    integer invd;
    real akkp1;
    extern logical lsame_(char *, char *);
    integer iinfo;
    extern /* Subroutine */
    int sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *);
    integer count;
    logical upper;
    extern /* Subroutine */
    int strmm_(char *, char *, char *, char *, integer *, integer *, real *, real *, integer *, real *, integer * );
    real u01_i_j__, u11_i_j__;
    extern /* Subroutine */
    int xerbla_(char *, integer *), strtri_( char *, char *, integer *, real *, integer *, integer *);
    real u01_ip1_j__, u11_ip1_j__;
    extern /* Subroutine */
    int ssyconv_(char *, char *, integer *, real *, integer *, integer *, real *, integer *);
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
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
        xerbla_("SSYTRI2X", &i__1);
        return 0;
    }
    if (*n == 0)
    {
        return 0;
    }
    /* Convert A */
    /* Workspace got Non-diag elements of D */
    ssyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[work_offset], & iinfo);
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
        }
    }
    *info = 0;
    /* Splitting Workspace */
    /* U01 is a block (N,NB+1) */
    /* The first element of U01 is in WORK(1,1) */
    /* U11 is a block (NB+1,NB+1) */
    /* The first element of U11 is in WORK(N+1,1) */
    u11 = *n;
    /* INVD is a block (N,2) */
    /* The first element of INVD is in WORK(1,INVD) */
    invd = *nb + 2;
    if (upper)
    {
        /* invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */
        strtri_(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D)*inv(U) */
        k = 1;
        while(k <= *n)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                work[k + invd * work_dim1] = 1.f / a[k + k * a_dim1];
                work[k + (invd + 1) * work_dim1] = 0.f;
                ++k;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                t = work[k + 1 + work_dim1];
                ak = a[k + k * a_dim1] / t;
                akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
                akkp1 = work[k + 1 + work_dim1] / t;
                d__ = t * (ak * akp1 - 1.f);
                work[k + invd * work_dim1] = akp1 / d__;
                work[k + 1 + (invd + 1) * work_dim1] = ak / d__;
                work[k + (invd + 1) * work_dim1] = -akkp1 / d__;
                work[k + 1 + invd * work_dim1] = -akkp1 / d__;
                k += 2;
            }
        }
        /* inv(U**T) = (inv(U))**T */
        /* inv(U**T)*inv(D)*inv(U) */
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
                count = 0;
                /* count negative elements, */
                i__1 = cut;
                for (i__ = cut + 1 - nnb;
                        i__ <= i__1;
                        ++i__)
                {
                    if (ipiv[i__] < 0)
                    {
                        ++count;
                    }
                }
                /* need a even number for a clear cut */
                if (count % 2 == 1)
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
                work[u11 + i__ + i__ * work_dim1] = 1.f;
                i__2 = i__ - 1;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    work[u11 + i__ + j * work_dim1] = 0.f;
                }
                i__2 = nnb;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    work[u11 + i__ + j * work_dim1] = a[cut + i__ + (cut + j) * a_dim1];
                }
            }
            /* invD*U01 */
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
                    ++i__;
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
                    i__ += 2;
                }
            }
            /* invD1*U11 */
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
                    ++i__;
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
                    i__ += 2;
                }
            }
            /* U11**T*invD1*U11->U11 */
            i__1 = *n + *nb + 1;
            strmm_("L", "U", "T", "U", &nnb, &nnb, &c_b11, &a[cut + 1 + (cut + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1);
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
            /* U01**T*invD*U01->A(CUT+I,CUT+J) */
            i__1 = *n + *nb + 1;
            i__2 = *n + *nb + 1;
            sgemm_("T", "N", &nnb, &nnb, &cut, &c_b11, &a[(cut + 1) * a_dim1 + 1], lda, &work[work_offset], &i__1, &c_b15, &work[u11 + 1 + work_dim1], &i__2);
            /* U11 = U11**T*invD1*U11 + U01**T*invD*U01 */
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
            /* U01 = U00**T*invD0*U01 */
            i__1 = *n + *nb + 1;
            strmm_("L", uplo, "T", "U", &cut, &nnb, &c_b11, &a[a_offset], lda, &work[work_offset], &i__1);
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
        /* Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T */
        i__ = 1;
        while(i__ <= *n)
        {
            if (ipiv[i__] > 0)
            {
                ip = ipiv[i__];
                if (i__ < ip)
                {
                    ssyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    ssyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
            else
            {
                ip = -ipiv[i__];
                ++i__;
                if (i__ - 1 < ip)
                {
                    i__1 = i__ - 1;
                    ssyswapr_(uplo, n, &a[a_offset], lda, &i__1, &ip);
                }
                if (i__ - 1 > ip)
                {
                    i__1 = i__ - 1;
                    ssyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__1);
                }
            }
            ++i__;
        }
    }
    else
    {
        /* LOWER... */
        /* invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */
        strtri_(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D)*inv(U) */
        k = *n;
        while(k >= 1)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                work[k + invd * work_dim1] = 1.f / a[k + k * a_dim1];
                work[k + (invd + 1) * work_dim1] = 0.f;
                --k;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                t = work[k - 1 + work_dim1];
                ak = a[k - 1 + (k - 1) * a_dim1] / t;
                akp1 = a[k + k * a_dim1] / t;
                akkp1 = work[k - 1 + work_dim1] / t;
                d__ = t * (ak * akp1 - 1.f);
                work[k - 1 + invd * work_dim1] = akp1 / d__;
                work[k + invd * work_dim1] = ak / d__;
                work[k + (invd + 1) * work_dim1] = -akkp1 / d__;
                work[k - 1 + (invd + 1) * work_dim1] = -akkp1 / d__;
                k += -2;
            }
        }
        /* inv(U**T) = (inv(U))**T */
        /* inv(U**T)*inv(D)*inv(U) */
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
                count = 0;
                /* count negative elements, */
                i__1 = cut + nnb;
                for (i__ = cut + 1;
                        i__ <= i__1;
                        ++i__)
                {
                    if (ipiv[i__] < 0)
                    {
                        ++count;
                    }
                }
                /* need a even number for a clear cut */
                if (count % 2 == 1)
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
                work[u11 + i__ + i__ * work_dim1] = 1.f;
                i__2 = nnb;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    work[u11 + i__ + j * work_dim1] = 0.f;
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
                    --i__;
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
                    i__ += -2;
                }
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
                    --i__;
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
                    i__ += -2;
                }
            }
            /* L11**T*invD1*L11->L11 */
            i__1 = *n + *nb + 1;
            strmm_("L", uplo, "T", "U", &nnb, &nnb, &c_b11, &a[cut + 1 + (cut + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1);
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
                /* L21**T*invD2*L21->A(CUT+I,CUT+J) */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                i__3 = *n + *nb + 1;
                sgemm_("T", "N", &nnb, &nnb, &i__1, &c_b11, &a[cut + nnb + 1 + (cut + 1) * a_dim1], lda, &work[work_offset], &i__2, &c_b15, &work[u11 + 1 + work_dim1], &i__3);
                /* L11 = L11**T*invD1*L11 + U01**T*invD*U01 */
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
                /* L01 = L22**T*invD2*L21 */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                strmm_("L", uplo, "T", "U", &i__1, &nnb, &c_b11, &a[cut + nnb + 1 + (cut + nnb + 1) * a_dim1], lda, &work[ work_offset], &i__2);
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
                /* L11 = L11**T*invD1*L11 */
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
        /* Apply PERMUTATIONS P and P**T: P * inv(U**T)*inv(D)*inv(U) *P**T */
        i__ = *n;
        while(i__ >= 1)
        {
            if (ipiv[i__] > 0)
            {
                ip = ipiv[i__];
                if (i__ < ip)
                {
                    ssyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    ssyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
            else
            {
                ip = -ipiv[i__];
                if (i__ < ip)
                {
                    ssyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    ssyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
                --i__;
            }
            --i__;
        }
    }
    return 0;
    /* End of SSYTRI2X */
}
/* ssytri2x_ */
