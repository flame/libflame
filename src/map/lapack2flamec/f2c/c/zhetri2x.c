/* ../netlib/zhetri2x.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    1.,0.
}
;
static doublecomplex c_b2 =
{
    0.,0.
}
;
/* > \brief \b ZHETRI2X */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHETRI2X + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri2 x.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri2 x.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri2 x.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHETRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N, NB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), WORK( N+NB+1,* ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRI2X computes the inverse of a COMPLEX*16 Hermitian indefinite matrix */
/* > A using the factorization A = U*D*U**H or A = L*D*L**H computed by */
/* > ZHETRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**H;
*/
/* > = 'L': Lower triangular, form is A = L*D*L**H. */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the NNB diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by ZHETRF. */
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
/* > as determined by ZHETRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N+NNB+1,NNB+3) */
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
/* > \ingroup complex16HEcomputational */
/* ===================================================================== */
/* Subroutine */
int zhetri2x_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *work, integer *nb, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, work_dim1, work_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    double z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg( doublecomplex *, doublecomplex *);
    /* Local variables */
    extern /* Subroutine */
    int zheswapr_(char *, integer *, doublecomplex *, integer *, integer *, integer *);
    doublecomplex d__;
    integer i__, j, k;
    doublecomplex t, ak;
    integer u11, ip, nnb, cut;
    doublecomplex akp1;
    integer invd;
    doublecomplex akkp1;
    extern logical lsame_(char *, char *);
    integer iinfo;
    extern /* Subroutine */
    int zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
    integer count;
    logical upper;
    extern /* Subroutine */
    int ztrmm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *);
    doublecomplex u01_i_j__, u11_i_j__;
    extern /* Subroutine */
    int xerbla_(char *, integer *), ztrtri_( char *, char *, integer *, doublecomplex *, integer *, integer *);
    doublecomplex u01_ip1_j__, u11_ip1_j__;
    extern /* Subroutine */
    int zsyconv_(char *, char *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *);
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
        xerbla_("ZHETRI2X", &i__1);
        return 0;
    }
    if (*n == 0)
    {
        return 0;
    }
    /* Convert A */
    /* Workspace got Non-diag elements of D */
    zsyconv_(uplo, "C", n, &a[a_offset], lda, &ipiv[1], &work[work_offset], & iinfo);
    /* Check that the diagonal matrix D is nonsingular. */
    if (upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for (*info = *n;
                *info >= 1;
                --(*info))
        {
            i__1 = *info + *info * a_dim1;
            if (ipiv[*info] > 0 && (a[i__1].r == 0. && a[i__1].i == 0.))
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
            if (ipiv[*info] > 0 && (a[i__2].r == 0. && a[i__2].i == 0.))
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
        /* invA = P * inv(U**H)*inv(D)*inv(U)*P**H. */
        ztrtri_(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D)*inv(U) */
        k = 1;
        while(k <= *n)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                i__1 = k + invd * work_dim1;
                i__2 = k + k * a_dim1;
                d__1 = 1.f / a[i__2].r;
                work[i__1].r = d__1;
                work[i__1].i = 0.; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                work[i__1].r = 0.;
                work[i__1].i = 0.; // , expr subst
                ++k;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                d__1 = z_abs(&work[k + 1 + work_dim1]);
                t.r = d__1;
                t.i = 0.; // , expr subst
                i__1 = k + k * a_dim1;
                d__1 = a[i__1].r;
                z__2.r = d__1;
                z__2.i = 0.; // , expr subst
                z_div(&z__1, &z__2, &t);
                ak.r = z__1.r;
                ak.i = z__1.i; // , expr subst
                i__1 = k + 1 + (k + 1) * a_dim1;
                d__1 = a[i__1].r;
                z__2.r = d__1;
                z__2.i = 0.; // , expr subst
                z_div(&z__1, &z__2, &t);
                akp1.r = z__1.r;
                akp1.i = z__1.i; // , expr subst
                z_div(&z__1, &work[k + 1 + work_dim1], &t);
                akkp1.r = z__1.r;
                akkp1.i = z__1.i; // , expr subst
                z__3.r = ak.r * akp1.r - ak.i * akp1.i;
                z__3.i = ak.r * akp1.i + ak.i * akp1.r; // , expr subst
                z__2.r = z__3.r - 1.f;
                z__2.i = z__3.i; // , expr subst
                z__1.r = t.r * z__2.r - t.i * z__2.i;
                z__1.i = t.r * z__2.i + t.i * z__2.r; // , expr subst
                d__.r = z__1.r;
                d__.i = z__1.i; // , expr subst
                i__1 = k + invd * work_dim1;
                z_div(&z__1, &akp1, &d__);
                work[i__1].r = z__1.r;
                work[i__1].i = z__1.i; // , expr subst
                i__1 = k + 1 + (invd + 1) * work_dim1;
                z_div(&z__1, &ak, &d__);
                work[i__1].r = z__1.r;
                work[i__1].i = z__1.i; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                z__2.r = -akkp1.r;
                z__2.i = -akkp1.i; // , expr subst
                z_div(&z__1, &z__2, &d__);
                work[i__1].r = z__1.r;
                work[i__1].i = z__1.i; // , expr subst
                i__1 = k + 1 + invd * work_dim1;
                d_cnjg(&z__1, &work[k + (invd + 1) * work_dim1]);
                work[i__1].r = z__1.r;
                work[i__1].i = z__1.i; // , expr subst
                k += 2;
            }
        }
        /* inv(U**H) = (inv(U))**H */
        /* inv(U**H)*inv(D)*inv(U) */
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
                work[i__2].r = 1.;
                work[i__2].i = 0.; // , expr subst
                i__2 = i__ - 1;
                for (j = 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    work[i__3].r = 0.;
                    work[i__3].i = 0.; // , expr subst
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
                        z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        z__1.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
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
                        z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        z__2.i = work[i__3].r * u01_i_j__.i + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = i__ + (invd + 1) * work_dim1;
                        z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        z__3.i = work[i__4].r * u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        z__1.r = z__2.r + z__3.r;
                        z__1.i = z__2.i + z__3.i; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
                        i__2 = i__ + 1 + j * work_dim1;
                        i__3 = i__ + 1 + invd * work_dim1;
                        z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        z__2.i = work[i__3].r * u01_i_j__.i + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = i__ + 1 + (invd + 1) * work_dim1;
                        z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        z__3.i = work[i__4].r * u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        z__1.r = z__2.r + z__3.r;
                        z__1.i = z__2.i + z__3.i; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
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
                        z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        z__1.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
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
                        z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        z__2.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        i__5 = cut + i__ + (invd + 1) * work_dim1;
                        i__6 = u11 + i__ + 1 + j * work_dim1;
                        z__3.r = work[i__5].r * work[i__6].r - work[i__5].i * work[i__6].i;
                        z__3.i = work[i__5].r * work[ i__6].i + work[i__5].i * work[i__6].r; // , expr subst
                        z__1.r = z__2.r + z__3.r;
                        z__1.i = z__2.i + z__3.i; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
                        i__2 = u11 + i__ + 1 + j * work_dim1;
                        i__3 = cut + i__ + 1 + invd * work_dim1;
                        z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * u11_i_j__.i;
                        z__2.i = work[i__3].r * u11_i_j__.i + work[i__3].i * u11_i_j__.r; // , expr subst
                        i__4 = cut + i__ + 1 + (invd + 1) * work_dim1;
                        z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i * u11_ip1_j__.i;
                        z__3.i = work[i__4].r * u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r; // , expr subst
                        z__1.r = z__2.r + z__3.r;
                        z__1.i = z__2.i + z__3.i; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
                    }
                    i__ += 2;
                }
            }
            /* U11**H*invD1*U11->U11 */
            i__1 = *n + *nb + 1;
            ztrmm_("L", "U", "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1);
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
            /* U01**H*invD*U01->A(CUT+I,CUT+J) */
            i__1 = *n + *nb + 1;
            i__2 = *n + *nb + 1;
            zgemm_("C", "N", &nnb, &nnb, &cut, &c_b1, &a[(cut + 1) * a_dim1 + 1], lda, &work[work_offset], &i__1, &c_b2, &work[u11 + 1 + work_dim1], &i__2);
            /* U11 = U11**H*invD1*U11 + U01**H*invD*U01 */
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
                    z__1.r = a[i__4].r + work[i__5].r;
                    z__1.i = a[i__4].i + work[i__5].i; // , expr subst
                    a[i__3].r = z__1.r;
                    a[i__3].i = z__1.i; // , expr subst
                }
            }
            /* U01 = U00**H*invD0*U01 */
            i__1 = *n + *nb + 1;
            ztrmm_("L", uplo, "C", "U", &cut, &nnb, &c_b1, &a[a_offset], lda, &work[work_offset], &i__1);
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
        /* Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H */
        i__ = 1;
        while(i__ <= *n)
        {
            if (ipiv[i__] > 0)
            {
                ip = ipiv[i__];
                if (i__ < ip)
                {
                    zheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    zheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
            else
            {
                ip = -ipiv[i__];
                ++i__;
                if (i__ - 1 < ip)
                {
                    i__1 = i__ - 1;
                    zheswapr_(uplo, n, &a[a_offset], lda, &i__1, &ip);
                }
                if (i__ - 1 > ip)
                {
                    i__1 = i__ - 1;
                    zheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__1);
                }
            }
            ++i__;
        }
    }
    else
    {
        /* LOWER... */
        /* invA = P * inv(U**H)*inv(D)*inv(U)*P**H. */
        ztrtri_(uplo, "U", n, &a[a_offset], lda, info);
        /* inv(D) and inv(D)*inv(U) */
        k = *n;
        while(k >= 1)
        {
            if (ipiv[k] > 0)
            {
                /* 1 x 1 diagonal NNB */
                i__1 = k + invd * work_dim1;
                i__2 = k + k * a_dim1;
                d__1 = 1.f / a[i__2].r;
                work[i__1].r = d__1;
                work[i__1].i = 0.; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                work[i__1].r = 0.;
                work[i__1].i = 0.; // , expr subst
                --k;
            }
            else
            {
                /* 2 x 2 diagonal NNB */
                d__1 = z_abs(&work[k - 1 + work_dim1]);
                t.r = d__1;
                t.i = 0.; // , expr subst
                i__1 = k - 1 + (k - 1) * a_dim1;
                d__1 = a[i__1].r;
                z__2.r = d__1;
                z__2.i = 0.; // , expr subst
                z_div(&z__1, &z__2, &t);
                ak.r = z__1.r;
                ak.i = z__1.i; // , expr subst
                i__1 = k + k * a_dim1;
                d__1 = a[i__1].r;
                z__2.r = d__1;
                z__2.i = 0.; // , expr subst
                z_div(&z__1, &z__2, &t);
                akp1.r = z__1.r;
                akp1.i = z__1.i; // , expr subst
                z_div(&z__1, &work[k - 1 + work_dim1], &t);
                akkp1.r = z__1.r;
                akkp1.i = z__1.i; // , expr subst
                z__3.r = ak.r * akp1.r - ak.i * akp1.i;
                z__3.i = ak.r * akp1.i + ak.i * akp1.r; // , expr subst
                z__2.r = z__3.r - 1.f;
                z__2.i = z__3.i; // , expr subst
                z__1.r = t.r * z__2.r - t.i * z__2.i;
                z__1.i = t.r * z__2.i + t.i * z__2.r; // , expr subst
                d__.r = z__1.r;
                d__.i = z__1.i; // , expr subst
                i__1 = k - 1 + invd * work_dim1;
                z_div(&z__1, &akp1, &d__);
                work[i__1].r = z__1.r;
                work[i__1].i = z__1.i; // , expr subst
                i__1 = k + invd * work_dim1;
                z_div(&z__1, &ak, &d__);
                work[i__1].r = z__1.r;
                work[i__1].i = z__1.i; // , expr subst
                i__1 = k + (invd + 1) * work_dim1;
                z__2.r = -akkp1.r;
                z__2.i = -akkp1.i; // , expr subst
                z_div(&z__1, &z__2, &d__);
                work[i__1].r = z__1.r;
                work[i__1].i = z__1.i; // , expr subst
                i__1 = k - 1 + (invd + 1) * work_dim1;
                d_cnjg(&z__1, &work[k + (invd + 1) * work_dim1]);
                work[i__1].r = z__1.r;
                work[i__1].i = z__1.i; // , expr subst
                k += -2;
            }
        }
        /* inv(U**H) = (inv(U))**H */
        /* inv(U**H)*inv(D)*inv(U) */
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
                work[i__2].r = 1.;
                work[i__2].i = 0.; // , expr subst
                i__2 = nnb;
                for (j = i__ + 1;
                        j <= i__2;
                        ++j)
                {
                    i__3 = u11 + i__ + j * work_dim1;
                    work[i__3].r = 0.;
                    work[i__3].i = 0.; // , expr subst
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
                        z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        z__1.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
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
                        z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        z__2.i = work[i__3].r * u01_i_j__.i + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = cut + nnb + i__ + (invd + 1) * work_dim1;
                        z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        z__3.i = work[i__4].r * u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        z__1.r = z__2.r + z__3.r;
                        z__1.i = z__2.i + z__3.i; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
                        i__2 = i__ - 1 + j * work_dim1;
                        i__3 = cut + nnb + i__ - 1 + (invd + 1) * work_dim1;
                        z__2.r = work[i__3].r * u01_i_j__.r - work[i__3].i * u01_i_j__.i;
                        z__2.i = work[i__3].r * u01_i_j__.i + work[i__3].i * u01_i_j__.r; // , expr subst
                        i__4 = cut + nnb + i__ - 1 + invd * work_dim1;
                        z__3.r = work[i__4].r * u01_ip1_j__.r - work[i__4].i * u01_ip1_j__.i;
                        z__3.i = work[i__4].r * u01_ip1_j__.i + work[i__4].i * u01_ip1_j__.r; // , expr subst
                        z__1.r = z__2.r + z__3.r;
                        z__1.i = z__2.i + z__3.i; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
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
                        z__1.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        z__1.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
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
                        z__2.r = work[i__3].r * work[i__4].r - work[i__3].i * work[i__4].i;
                        z__2.i = work[i__3].r * work[ i__4].i + work[i__3].i * work[i__4].r; // , expr subst
                        i__5 = cut + i__ + (invd + 1) * work_dim1;
                        z__3.r = work[i__5].r * u11_ip1_j__.r - work[i__5].i * u11_ip1_j__.i;
                        z__3.i = work[i__5].r * u11_ip1_j__.i + work[i__5].i * u11_ip1_j__.r; // , expr subst
                        z__1.r = z__2.r + z__3.r;
                        z__1.i = z__2.i + z__3.i; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
                        i__2 = u11 + i__ - 1 + j * work_dim1;
                        i__3 = cut + i__ - 1 + (invd + 1) * work_dim1;
                        z__2.r = work[i__3].r * u11_i_j__.r - work[i__3].i * u11_i_j__.i;
                        z__2.i = work[i__3].r * u11_i_j__.i + work[i__3].i * u11_i_j__.r; // , expr subst
                        i__4 = cut + i__ - 1 + invd * work_dim1;
                        z__3.r = work[i__4].r * u11_ip1_j__.r - work[i__4].i * u11_ip1_j__.i;
                        z__3.i = work[i__4].r * u11_ip1_j__.i + work[i__4].i * u11_ip1_j__.r; // , expr subst
                        z__1.r = z__2.r + z__3.r;
                        z__1.i = z__2.i + z__3.i; // , expr subst
                        work[i__2].r = z__1.r;
                        work[i__2].i = z__1.i; // , expr subst
                    }
                    i__ += -2;
                }
            }
            /* L11**H*invD1*L11->L11 */
            i__1 = *n + *nb + 1;
            ztrmm_("L", uplo, "C", "U", &nnb, &nnb, &c_b1, &a[cut + 1 + (cut + 1) * a_dim1], lda, &work[u11 + 1 + work_dim1], &i__1);
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
                /* L21**H*invD2*L21->A(CUT+I,CUT+J) */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                i__3 = *n + *nb + 1;
                zgemm_("C", "N", &nnb, &nnb, &i__1, &c_b1, &a[cut + nnb + 1 + (cut + 1) * a_dim1], lda, &work[work_offset], &i__2, & c_b2, &work[u11 + 1 + work_dim1], &i__3);
                /* L11 = L11**H*invD1*L11 + U01**H*invD*U01 */
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
                        z__1.r = a[i__4].r + work[i__5].r;
                        z__1.i = a[i__4].i + work[i__5].i; // , expr subst
                        a[i__3].r = z__1.r;
                        a[i__3].i = z__1.i; // , expr subst
                    }
                }
                /* L01 = L22**H*invD2*L21 */
                i__1 = *n - nnb - cut;
                i__2 = *n + *nb + 1;
                ztrmm_("L", uplo, "C", "U", &i__1, &nnb, &c_b1, &a[cut + nnb + 1 + (cut + nnb + 1) * a_dim1], lda, &work[ work_offset], &i__2);
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
                /* L11 = L11**H*invD1*L11 */
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
        /* Apply PERMUTATIONS P and P**H: P * inv(U**H)*inv(D)*inv(U) *P**H */
        i__ = *n;
        while(i__ >= 1)
        {
            if (ipiv[i__] > 0)
            {
                ip = ipiv[i__];
                if (i__ < ip)
                {
                    zheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    zheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
            }
            else
            {
                ip = -ipiv[i__];
                if (i__ < ip)
                {
                    zheswapr_(uplo, n, &a[a_offset], lda, &i__, &ip);
                }
                if (i__ > ip)
                {
                    zheswapr_(uplo, n, &a[a_offset], lda, &ip, &i__);
                }
                --i__;
            }
            --i__;
        }
    }
    return 0;
    /* End of ZHETRI2X */
}
/* zhetri2x_ */
