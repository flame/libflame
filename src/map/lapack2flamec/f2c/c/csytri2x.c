/* ../netlib/csytri2x.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    1.f,0.f
}
;
static complex c_b2 =
{
    0.f,0.f
}
;
/* > \brief \b CSYTRI2X */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSYTRI2X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytri2 x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytri2 x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytri2 x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N, NB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), WORK( N+NB+1,* ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYTRI2X computes the inverse of a real symmetric indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > CSYTRF. */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the NNB diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by CSYTRF. */
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
/* > as determined by CSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N+NNB+1,NNB+3) */
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
/* > \ingroup complexSYcomputational */
/* ===================================================================== */
/* Subroutine */
int csytri2x_(char *uplo, integer *n, complex *a, integer * lda, integer *ipiv, complex *work, integer *nb, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    complex q__1, q__2, q__3;
    /* Builtin functions */
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    extern /* Subroutine */
    int csyswapr_(char *, integer *, complex *, integer *, integer *, integer *);
    complex d__;
    integer i__, j, k;
    complex t, ak;
    integer u11, ip, nnb, cut;
    complex akp1;
    integer invd;
    complex akkp1;
    extern /* Subroutine */
    int cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *);
    integer iinfo;
    extern /* Subroutine */
    int ctrmm_(char *, char *, char *, char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *);
    integer count;
    logical upper;
    complex u01_i_j__, u11_i_j__;
    extern /* Subroutine */
    int xerbla_(char *, integer *), ctrtri_( char *, char *, integer *, complex *, integer *, integer *), csyconv_(char *, char *, integer *, complex *, integer *, integer *, complex *, integer *);
    complex u01_ip1_j__, u11_ip1_j__;
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
        xerbla_("CSYTRI2X", &i__1);
        return 0;
    }
    if (*n == 0)
    {
        return 0;
    }
    /* Convert A */
    /* Workspace got Non-diag elements of D */
    csyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[work_offset], & iinfo);
    /* Check that the diagonal matrix D is nonsingular. */
    if (upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for (*info = *n;
                *info >= 1;
                --(*info))
        {
            i__1 = *info + *info * a_dim1;
            if (ipiv[*info] > 0 && (a[i__1].r == 0.f && a[i__1].i == 0.f))
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
            i__2 = *info + *info * a_dim1;
            if (ipiv[*info] > 0 && (a[i__2].r == 0.f && a[i__2].i == 0.f))
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
        ctrtri_(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D)*inv(U) */
        k = 1;
        while(k <= *n)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                i__1 = k + invd * work_dim1;
                c_div(&q__1, &c_b1, &a[k + k * a_dim1]);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                work[i__1].r = 0.f;
                work[i__1].i = 0.f; // , expr subst
                ++k;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                i__1 = k + 1 + work_dim1;
                t.r = work[i__1].r;
                t.i = work[i__1].i; // , expr subst
                c_div(&q__1, &a[k + k * a_dim1], &t);
                ak.r = q__1.r;
                ak.i = q__1.i; // , expr subst
                c_div(&q__1, &a[k + 1 + (k + 1) * a_dim1], &t);
                akp1.r = q__1.r;
                akp1.i = q__1.i; // , expr subst
                c_div(&q__1, &work[k + 1 + work_dim1], &t);
                akkp1.r = q__1.r;
                akkp1.i = q__1.i; // , expr subst
                q__3.r = ak.r * akp1.r - ak.i * akp1.i;
                q__3.i = ak.r * akp1.i + ak.i * akp1.r; // , expr subst
                q__2.r = q__3.r - 1.f;
                q__2.i = q__3.i - 0.f; // , expr subst
                q__1.r = t.r * q__2.r - t.i * q__2.i;
                q__1.i = t.r * q__2.i + t.i * q__2.r; // , expr subst
                d__.r = q__1.r;
                d__.i = q__1.i; // , expr subst
                i__1 = k + invd * work_dim1;
                c_div(&q__1, &akp1, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + 1 + (invd + 1) * work_dim1;
                c_div(&q__1, &ak, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                q__2.r = -akkp1.r;
                q__2.i = -akkp1.i; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + 1 + invd * work_dim1;
                q__2.r = -akkp1.r;
                q__2.i = -akkp1.i; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
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
                    i__3 = i__ + j * work_dim1;
                    i__4 = i__ + (cut + j) * a_dim1;
                    work[i__3].r = a[i__4].r;
                    work[i__3].i = a[i__4].i; // , expr subst
                }
            }
            /* U11 Block */
            i__1 = nnb;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = u11 + i__ + i__ * work_dim1;
                work[i__2].r = 1.f;
                work[i__2].i = 0.f; // , expr subst
                i__2 = i__ - 1;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    work[i__3].r = 0.f;
                    work[i__3].i = 0.f; // , expr subst
                }
                i__2 = nnb;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    i__4 = cut + i__ + (cut + j) * a_dim1;
                    work[i__3].r = a[i__4].r;
                    work[i__3].i = a[i__4].i; // , expr subst
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
                        i__2 = i__ + j * work_dim1;
                        i__3 = i__ + invd * work_dim1;
                        i__4 = i__ + j * work_dim1;
                        q__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__1.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
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
                        i__2 = i__ + j * work_dim1;
                        u01_i_j__.r = work[i__2].r;
                        u01_i_j__.i = work[i__2] .i; // , expr subst
                        i__2 = i__ + 1 + j * work_dim1;
                        u01_ip1_j__.r = work[i__2].r;
                        u01_ip1_j__.i = work[ i__2].i; // , expr subst
                        i__2 = i__ + j * work_dim1;
                        i__3 = i__ + invd * work_dim1;
                        q__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        q__2.i = work[i__3].r * u01_i_j__.i + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = i__ + (invd + 1) * work_dim1;
                        q__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        q__3.i = work[i__4].r * u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                        i__2 = i__ + 1 + j * work_dim1;
                        i__3 = i__ + 1 + invd * work_dim1;
                        q__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        q__2.i = work[i__3].r * u01_i_j__.i + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = i__ + 1 + (invd + 1) * work_dim1;
                        q__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        q__3.i = work[i__4].r * u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
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
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__1.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
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
                        i__2 = u11 + i__ + j * work_dim1;
                        u11_i_j__.r = work[i__2].r;
                        u11_i_j__.i = work[i__2] .i; // , expr subst
                        i__2 = u11 + i__ + 1 + j * work_dim1;
                        u11_ip1_j__.r = work[i__2].r;
                        u11_ip1_j__.i = work[ i__2].i; // , expr subst
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__2.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__2.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        i__5 = cut + i__ + (invd + 1) * work_dim1;
                        i__6 = u11 + i__ + 1 + j * work_dim1;
                        q__3.r = work[i__5].r * work[i__6].r - work[i__5].i * work[i__6].i;
                        q__3.i = work[i__5].r * work[ i__6].i + work[i__5].i * work[i__6].r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                        i__2 = u11 + i__ + 1 + j * work_dim1;
                        i__3 = cut + i__ + 1 + invd * work_dim1;
                        q__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * u11_i_j__.i;
                        q__2.i = work[i__3].r * u11_i_j__.i + work[i__3].i * u11_i_j__.r; // , expr subst
                        i__4 = cut + i__ + 1 + (invd + 1) * work_dim1;
                        q__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i * u11_ip1_j__.i;
                        q__3.i = work[i__4].r * u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                    }
                    i__ += 2;
                }
            }
            /* U11**T*invD1*U11->U11 */
            i__1 = *n + *nb + 1;
            ctrmm_("L", "U", "T", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1);
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
                    i__3 = cut + i__ + (cut + j) * a_dim1;
                    i__4 = u11 + i__ + j * work_dim1;
                    a[i__3].r = work[i__4].r;
                    a[i__3].i = work[i__4].i; // , expr subst
                }
            }
            /* U01**T*invD*U01->A(CUT+I,CUT+J) */
            i__1 = *n + *nb + 1;
            i__2 = *n + *nb + 1;
            cgemm_("T", "N", &nnb, &nnb, &cut, &c_b1, &a[(cut + 1) * a_dim1 + 1], lda, &work[work_offset], &i__1, &c_b2, &work[u11 + 1 + work_dim1], &i__2);
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
                    i__3 = cut + i__ + (cut + j) * a_dim1;
                    i__4 = cut + i__ + (cut + j) * a_dim1;
                    i__5 = u11 + i__ + j * work_dim1;
                    q__1.r = a[i__4].r + work[i__5].r;
                    q__1.i = a[i__4].i + work[i__5].i; // , expr subst
                    a[i__3].r = q__1.r;
                    a[i__3].i = q__1.i; // , expr subst
                }
            }
            /* U01 = U00**T*invD0*U01 */
            i__1 = *n + *nb + 1;
            ctrmm_("L", uplo, "T", "U", &cut, &nnb, &c_b1, &a[a_offset], lda, &work[work_offset], &i__1);
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
                    i__3 = i__ + (cut + j) * a_dim1;
                    i__4 = i__ + j * work_dim1;
                    a[i__3].r = work[i__4].r;
                    a[i__3].i = work[i__4].i; // , expr subst
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
                    csyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    csyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
            else
            {
                ip = -ipiv[i__];
                ++i__;
                if (i__ - 1 < ip)
                {
                    i__1 = i__ - 1;
                    csyswapr_(uplo, n, &a[a_offset], lda, &i__1, &ip);
                }
                if (i__ - 1 > ip)
                {
                    i__1 = i__ - 1;
                    csyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__1);
                }
            }
            ++i__;
        }
    }
    else
    {
        /* LOWER... */
        /* invA = P * inv(U**T)*inv(D)*inv(U)*P**T. */
        ctrtri_(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D)*inv(U) */
        k = *n;
        while(k >= 1)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                i__1 = k + invd * work_dim1;
                c_div(&q__1, &c_b1, &a[k + k * a_dim1]);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                work[i__1].r = 0.f;
                work[i__1].i = 0.f; // , expr subst
                --k;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                i__1 = k - 1 + work_dim1;
                t.r = work[i__1].r;
                t.i = work[i__1].i; // , expr subst
                c_div(&q__1, &a[k - 1 + (k - 1) * a_dim1], &t);
                ak.r = q__1.r;
                ak.i = q__1.i; // , expr subst
                c_div(&q__1, &a[k + k * a_dim1], &t);
                akp1.r = q__1.r;
                akp1.i = q__1.i; // , expr subst
                c_div(&q__1, &work[k - 1 + work_dim1], &t);
                akkp1.r = q__1.r;
                akkp1.i = q__1.i; // , expr subst
                q__3.r = ak.r * akp1.r - ak.i * akp1.i;
                q__3.i = ak.r * akp1.i + ak.i * akp1.r; // , expr subst
                q__2.r = q__3.r - 1.f;
                q__2.i = q__3.i - 0.f; // , expr subst
                q__1.r = t.r * q__2.r - t.i * q__2.i;
                q__1.i = t.r * q__2.i + t.i * q__2.r; // , expr subst
                d__.r = q__1.r;
                d__.i = q__1.i; // , expr subst
                i__1 = k - 1 + invd * work_dim1;
                c_div(&q__1, &akp1, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + invd * work_dim1;
                c_div(&q__1, &ak, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                q__2.r = -akkp1.r;
                q__2.i = -akkp1.i; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                i__1 = k - 1 + (invd + 1) * work_dim1;
                q__2.r = -akkp1.r;
                q__2.i = -akkp1.i; // , expr subst
                c_div(&q__1, &q__2, &d__);
                work[i__1].r = q__1.r;
                work[i__1].i = q__1.i; // , expr subst
                k += -2;
            }
        }
        /* inv(U**T) = (inv(U))**T */
        /* inv(U**T)*inv(D)*inv(U) */
        cut = 0;
        while(cut < *n)
        {
            nnb = *nb;
            if (cut + nnb >= *n)
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
                    i__3 = i__ + j * work_dim1;
                    i__4 = cut + nnb + i__ + (cut + j) * a_dim1;
                    work[i__3].r = a[i__4].r;
                    work[i__3].i = a[i__4].i; // , expr subst
                }
            }
            /* L11 Block */
            i__1 = nnb;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = u11 + i__ + i__ * work_dim1;
                work[i__2].r = 1.f;
                work[i__2].i = 0.f; // , expr subst
                i__2 = nnb;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    work[i__3].r = 0.f;
                    work[i__3].i = 0.f; // , expr subst
                }
                i__2 = i__ - 1;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    i__4 = cut + i__ + (cut + j) * a_dim1;
                    work[i__3].r = a[i__4].r;
                    work[i__3].i = a[i__4].i; // , expr subst
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
                        i__2 = i__ + j * work_dim1;
                        i__3 = cut + nnb + i__ + invd * work_dim1;
                        i__4 = i__ + j * work_dim1;
                        q__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__1.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
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
                        i__2 = i__ + j * work_dim1;
                        u01_i_j__.r = work[i__2].r;
                        u01_i_j__.i = work[i__2] .i; // , expr subst
                        i__2 = i__ - 1 + j * work_dim1;
                        u01_ip1_j__.r = work[i__2].r;
                        u01_ip1_j__.i = work[ i__2].i; // , expr subst
                        i__2 = i__ + j * work_dim1;
                        i__3 = cut + nnb + i__ + invd * work_dim1;
                        q__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        q__2.i = work[i__3].r * u01_i_j__.i + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = cut + nnb + i__ + (invd + 1) * work_dim1;
                        q__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        q__3.i = work[i__4].r * u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                        i__2 = i__ - 1 + j * work_dim1;
                        i__3 = cut + nnb + i__ - 1 + (invd + 1) * work_dim1;
                        q__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        q__2.i = work[i__3].r * u01_i_j__.i + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = cut + nnb + i__ - 1 + invd * work_dim1;
                        q__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        q__3.i = work[i__4].r * u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
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
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__1.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
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
                        i__2 = u11 + i__ + j * work_dim1;
                        u11_i_j__.r = work[i__2].r;
                        u11_i_j__.i = work[i__2] .i; // , expr subst
                        i__2 = u11 + i__ - 1 + j * work_dim1;
                        u11_ip1_j__.r = work[i__2].r;
                        u11_ip1_j__.i = work[ i__2].i; // , expr subst
                        i__2 = u11 + i__ + j * work_dim1;
                        i__3 = cut + i__ + invd * work_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        q__2.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        q__2.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        i__5 = cut + i__ + (invd + 1) * work_dim1;
                        q__3.r = work[i__5].r * u11_ip1_j__.r - work[i__5].i * u11_ip1_j__.i;
                        q__3.i = work[i__5].r * u11_ip1_j__.i + work[i__5].i * u11_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                        i__2 = u11 + i__ - 1 + j * work_dim1;
                        i__3 = cut + i__ - 1 + (invd + 1) * work_dim1;
                        q__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * u11_i_j__.i;
                        q__2.i = work[i__3].r * u11_i_j__.i + work[i__3].i * u11_i_j__.r; // , expr subst
                        i__4 = cut + i__ - 1 + invd * work_dim1;
                        q__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i * u11_ip1_j__.i;
                        q__3.i = work[i__4].r * u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r; // , expr subst
                        q__1.r = q__2.r + q__3.r;
                        q__1.i = q__2.i + q__3.i; // , expr subst
                        work[i__2].r = q__1.r;
                        work[i__2].i = q__1.i; // , expr subst
                    }
                    i__ += -2;
                }
            }
            /* L11**T*invD1*L11->L11 */
            i__1 = *n + *nb + 1;
            ctrmm_("L", uplo, "T", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1);
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
                    i__3 = cut + i__ + (cut + j) * a_dim1;
                    i__4 = u11 + i__ + j * work_dim1;
                    a[i__3].r = work[i__4].r;
                    a[i__3].i = work[i__4].i; // , expr subst
                }
            }
            if (cut + nnb < *n)
            {
                /* L21**T*invD2*L21->A(CUT+I,CUT+J) */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                i__3 = *n + *nb + 1;
                cgemm_("T", "N", &nnb, &nnb, &i__1, &c_b1, &a[cut + nnb + 1 + (cut + 1) * a_dim1], lda, &work[work_offset], &i__2, & c_b2, &work[u11 + 1 + work_dim1], &i__3);
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
                        i__3 = cut + i__ + (cut + j) * a_dim1;
                        i__4 = cut + i__ + (cut + j) * a_dim1;
                        i__5 = u11 + i__ + j * work_dim1;
                        q__1.r = a[i__4].r + work[i__5].r;
                        q__1.i = a[i__4].i + work[i__5].i; // , expr subst
                        a[i__3].r = q__1.r;
                        a[i__3].i = q__1.i; // , expr subst
                    }
                }
                /* L01 = L22**T*invD2*L21 */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                ctrmm_("L", uplo, "T", "U", &i__1, &nnb, &c_b1, &a[cut + nnb + 1 + (cut + nnb + 1) * a_dim1], lda, &work[ work_offset], &i__2);
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
                        i__3 = cut + nnb + i__ + (cut + j) * a_dim1;
                        i__4 = i__ + j * work_dim1;
                        a[i__3].r = work[i__4].r;
                        a[i__3].i = work[i__4].i; // , expr subst
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
                        i__3 = cut + i__ + (cut + j) * a_dim1;
                        i__4 = u11 + i__ + j * work_dim1;
                        a[i__3].r = work[i__4].r;
                        a[i__3].i = work[i__4].i; // , expr subst
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
                    csyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    csyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
            else
            {
                ip = -ipiv[i__];
                if (i__ < ip)
                {
                    csyswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    csyswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
                --i__;
            }
            --i__;
        }
    }
    return 0;
    /* End of CSYTRI2X */
}
/* csytri2x_ */
