/* ../netlib/v3.9.0/ctrevc3.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    0.f,0.f
}
;
static complex c_b2 =
{
    1.f,0.f
}
;
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b CTREVC3 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CTREVC3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrevc3 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrevc3 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrevc3 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, */
/* LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO) */
/* .. Scalar Arguments .. */
/* CHARACTER HOWMNY, SIDE */
/* INTEGER INFO, LDT, LDVL, LDVR, LWORK, M, MM, N */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL SELECT( * ) */
/* REAL RWORK( * ) */
/* COMPLEX T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTREVC3 computes some or all of the right and/or left eigenvectors of */
/* > a complex upper triangular matrix T. */
/* > Matrices of this type are produced by the Schur factorization of */
/* > a complex general matrix: A = Q*T*Q**H, as computed by CHSEQR. */
/* > */
/* > The right eigenvector x and the left eigenvector y of T corresponding */
/* > to an eigenvalue w are defined by: */
/* > */
/* > T*x = w*x, (y**H)*T = w*(y**H) */
/* > */
/* > where y**H denotes the conjugate transpose of the vector y. */
/* > The eigenvalues are not input to this routine, but are read directly */
/* > from the diagonal of T. */
/* > */
/* > This routine returns the matrices X and/or Y of right and left */
/* > eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an */
/* > input matrix. If Q is the unitary factor that reduces a matrix A to */
/* > Schur form T, then Q*X and Q*Y are the matrices of right and left */
/* > eigenvectors of A. */
/* > */
/* > This uses a Level 3 BLAS version of the back transformation. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'R': compute right eigenvectors only;
*/
/* > = 'L': compute left eigenvectors only;
*/
/* > = 'B': compute both right and left eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] HOWMNY */
/* > \verbatim */
/* > HOWMNY is CHARACTER*1 */
/* > = 'A': compute all right and/or left eigenvectors;
*/
/* > = 'B': compute all right and/or left eigenvectors, */
/* > backtransformed using the matrices supplied in */
/* > VR and/or VL;
*/
/* > = 'S': compute selected right and/or left eigenvectors, */
/* > as indicated by the logical array SELECT. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* > SELECT is LOGICAL array, dimension (N) */
/* > If HOWMNY = 'S', SELECT specifies the eigenvectors to be */
/* > computed. */
/* > The eigenvector corresponding to the j-th eigenvalue is */
/* > computed if SELECT(j) = .TRUE.. */
/* > Not referenced if HOWMNY = 'A' or 'B'. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix T. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] T */
/* > \verbatim */
/* > T is COMPLEX array, dimension (LDT,N) */
/* > The upper triangular matrix T. T is modified, but restored */
/* > on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* > VL is COMPLEX array, dimension (LDVL,MM) */
/* > On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/* > contain an N-by-N matrix Q (usually the unitary matrix Q of */
/* > Schur vectors returned by CHSEQR). */
/* > On exit, if SIDE = 'L' or 'B', VL contains: */
/* > if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
*/
/* > if HOWMNY = 'B', the matrix Q*Y;
*/
/* > if HOWMNY = 'S', the left eigenvectors of T specified by */
/* > SELECT, stored consecutively in the columns */
/* > of VL, in the same order as their */
/* > eigenvalues. */
/* > Not referenced if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* > LDVL is INTEGER */
/* > The leading dimension of the array VL. */
/* > LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* > VR is COMPLEX array, dimension (LDVR,MM) */
/* > On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/* > contain an N-by-N matrix Q (usually the unitary matrix Q of */
/* > Schur vectors returned by CHSEQR). */
/* > On exit, if SIDE = 'R' or 'B', VR contains: */
/* > if HOWMNY = 'A', the matrix X of right eigenvectors of T;
*/
/* > if HOWMNY = 'B', the matrix Q*X;
*/
/* > if HOWMNY = 'S', the right eigenvectors of T specified by */
/* > SELECT, stored consecutively in the columns */
/* > of VR, in the same order as their */
/* > eigenvalues. */
/* > Not referenced if SIDE = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* > LDVR is INTEGER */
/* > The leading dimension of the array VR. */
/* > LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] MM */
/* > \verbatim */
/* > MM is INTEGER */
/* > The number of columns in the arrays VL and/or VR. MM >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of columns in the arrays VL and/or VR actually */
/* > used to store the eigenvectors. */
/* > If HOWMNY = 'A' or 'B', M is set to N. */
/* > Each selected eigenvector occupies one column. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of array WORK. LWORK >= fla_max(1,2*N). */
/* > For optimum performance, LWORK >= N + 2*N*NB, where NB is */
/* > the optimal blocksize. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (LRWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* > LRWORK is INTEGER */
/* > The dimension of array RWORK. LRWORK >= fla_max(1,N). */
/* > */
/* > If LRWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the RWORK array, returns */
/* > this value as the first entry of the RWORK array, and no error */
/* > message related to LRWORK is issued by XERBLA. */
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
/* > \date November 2017 */
/* @generated from ztrevc3.f, fortran z -> c, Tue Apr 19 01:47:44 2016 */
/* > \ingroup complexOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The algorithm used in this program is basically backward (forward) */
/* > substitution, with scaling to make the the code robust against */
/* > possible overflow. */
/* > */
/* > Each eigenvector is normalized so that the element of largest */
/* > magnitude has magnitude 1;
here the magnitude of a complex number */
/* > (x,y) is taken to be |x| + |y|. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int ctrevc3_(char *side, char *howmny, logical *select, integer *n, complex *t, integer *ldt, complex *vl, integer *ldvl, complex *vr, integer *ldvr, integer *mm, integer *m, complex *work, integer *lwork, real *rwork, integer *lrwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"ctrevc3 inputs: side %c, howmny %c, n %lld, ldt %lld, ldvl %lld, ldvr %lld, mm %lld, lwork %lld, lrwork %lld",*side, *howmny, *n, *ldt, *ldvl, *ldvr, *mm, *lwork, *lrwork);
#else
    snprintf(buffer, 256,"ctrevc3 inputs: side %c, howmny %c, n %d, ldt %d, ldvl %d, ldvr %d, mm %d, lwork %d, lrwork %d",*side, *howmny, *n, *ldt, *ldvl, *ldvr, *mm, *lwork, *lrwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__3, i__4, i__5, i__6;
    real r__1, r__2, r__3;
    complex q__1, q__2;
    char ch__1[2];
    /* Local variables */
    integer i__, j, k, nb, ii, ki, is, iv;
    real ulp;
    logical allv;
    real unfl, ovfl, smin;
    logical over;
    real scale;
    extern /* Subroutine */
    int cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int cgemv_(char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *, complex *, complex *, integer *);
    real remax;
    extern /* Subroutine */
    int ccopy_(integer *, complex *, integer *, complex *, integer *);
    logical leftv, bothv, somev;
    extern /* Subroutine */
    int slabad_(real *, real *);
    extern integer icamax_(integer *, complex *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int csscal_(integer *, real *, complex *, integer *), claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *), clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int clatrs_(char *, char *, char *, char *, integer *, complex *, integer *, complex *, real *, real *, integer *);
    extern real scasum_(integer *, complex *, integer *);
    logical rightv;
    integer maxwrk;
    real smlnum;
    logical lquery;
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Decode and test the input parameters */
    /* Parameter adjustments */
    --select;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    --rwork;
    /* Function Body */
    bothv = lsame_(side, "B");
    rightv = lsame_(side, "R") || bothv;
    leftv = lsame_(side, "L") || bothv;
    allv = lsame_(howmny, "A");
    over = lsame_(howmny, "B");
    somev = lsame_(howmny, "S");
    /* Set M to the number of columns required to store the selected */
    /* eigenvectors. */
    if (somev)
    {
        *m = 0;
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            if (select[j])
            {
                ++(*m);
            }
            /* L10: */
        }
    }
    else
    {
        *m = *n;
    }
    *info = 0;
    nb = ilaenv_(&c__1, "CTREVC", ch__1, n, &c_n1, &c_n1, &c_n1);
    maxwrk = *n + (*n << 1) * nb;
    work[1].r = (real) maxwrk;
    work[1].i = 0.f; // , expr subst
    rwork[1] = (real) (*n);
    lquery = *lwork == -1 || *lrwork == -1;
    if (! rightv && ! leftv)
    {
        *info = -1;
    }
    else if (! allv && ! over && ! somev)
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*ldt < fla_max(1,*n))
    {
        *info = -6;
    }
    else if (*ldvl < 1 || leftv && *ldvl < *n)
    {
        *info = -8;
    }
    else if (*ldvr < 1 || rightv && *ldvr < *n)
    {
        *info = -10;
    }
    else if (*mm < *m)
    {
        *info = -11;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__3 = *n << 1; // , expr subst
        if (*lwork < fla_max(i__1,i__3) && ! lquery)
        {
            *info = -14;
        }
        else if (*lrwork < fla_max(1,*n) && ! lquery)
        {
            *info = -16;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CTREVC3", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Use blocked version of back-transformation if sufficient workspace. */
    /* Zero-out the workspace to avoid potential NaN propagation. */
    if (over && *lwork >= *n + (*n << 4))
    {
        nb = (*lwork - *n) / (*n << 1);
        nb = fla_min(nb,128);
        i__1 = (nb << 1) + 1;
        claset_("F", n, &i__1, &c_b1, &c_b1, &work[1], n);
    }
    else
    {
        nb = 1;
    }
    /* Set the constants to control overflow. */
    unfl = slamch_("Safe minimum");
    ovfl = 1.f / unfl;
    slabad_(&unfl, &ovfl);
    ulp = slamch_("Precision");
    smlnum = unfl * (*n / ulp);
    /* Store the diagonal elements of T in working array WORK. */
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__3 = i__;
        i__4 = i__ + i__ * t_dim1;
        work[i__3].r = t[i__4].r;
        work[i__3].i = t[i__4].i; // , expr subst
        /* L20: */
    }
    /* Compute 1-norm of each column of strictly upper triangular */
    /* part of T to control overflow in triangular solver. */
    rwork[1] = 0.f;
    i__1 = *n;
    for (j = 2;
            j <= i__1;
            ++j)
    {
        i__3 = j - 1;
        rwork[j] = scasum_(&i__3, &t[j * t_dim1 + 1], &c__1);
        /* L30: */
    }
    if (rightv)
    {
        /* ============================================================ */
        /* Compute right eigenvectors. */
        /* IV is index of column in current block. */
        /* Non-blocked version always uses IV=NB=1;
        */
        /* blocked version starts with IV=NB, goes down to 1. */
        /* (Note the "0-th" column is used to store the original diagonal.) */
        iv = nb;
        is = *m;
        for (ki = *n;
                ki >= 1;
                --ki)
        {
            if (somev)
            {
                if (! select[ki])
                {
                    goto L80;
                }
            }
            /* Computing MAX */
            i__1 = ki + ki * t_dim1;
            r__3 = ulp * ((r__1 = t[i__1].r, f2c_abs(r__1)) + (r__2 = t[ ki + ki * t_dim1].i, f2c_abs(r__2)));
            smin = fla_max(r__3,smlnum);
            /* -------------------------------------------------------- */
            /* Complex right eigenvector */
            i__1 = ki + iv * *n;
            work[i__1].r = 1.f;
            work[i__1].i = 0.f; // , expr subst
            /* Form right-hand side. */
            i__1 = ki - 1;
            for (k = 1;
                    k <= i__1;
                    ++k)
            {
                i__3 = k + iv * *n;
                i__4 = k + ki * t_dim1;
                q__1.r = -t[i__4].r;
                q__1.i = -t[i__4].i; // , expr subst
                work[i__3].r = q__1.r;
                work[i__3].i = q__1.i; // , expr subst
                /* L40: */
            }
            /* Solve upper triangular system: */
            /* [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE*WORK. */
            i__1 = ki - 1;
            for (k = 1;
                    k <= i__1;
                    ++k)
            {
                i__3 = k + k * t_dim1;
                i__4 = k + k * t_dim1;
                i__5 = ki + ki * t_dim1;
                q__1.r = t[i__4].r - t[i__5].r;
                q__1.i = t[i__4].i - t[i__5] .i; // , expr subst
                t[i__3].r = q__1.r;
                t[i__3].i = q__1.i; // , expr subst
                i__3 = k + k * t_dim1;
                if ((r__1 = t[i__3].r, f2c_abs(r__1)) + (r__2 = t[k + k * t_dim1].i, f2c_abs(r__2)) < smin)
                {
                    i__4 = k + k * t_dim1;
                    t[i__4].r = smin;
                    t[i__4].i = 0.f; // , expr subst
                }
                /* L50: */
            }
            if (ki > 1)
            {
                i__1 = ki - 1;
                clatrs_("Upper", "No transpose", "Non-unit", "Y", &i__1, &t[ t_offset], ldt, &work[iv * *n + 1], &scale, &rwork[1], info);
                i__1 = ki + iv * *n;
                work[i__1].r = scale;
                work[i__1].i = 0.f; // , expr subst
            }
            /* Copy the vector x or Q*x to VR and normalize. */
            if (! over)
            {
                /* ------------------------------ */
                /* no back-transform: copy x to VR and normalize. */
                ccopy_(&ki, &work[iv * *n + 1], &c__1, &vr[is * vr_dim1 + 1], &c__1);
                ii = icamax_(&ki, &vr[is * vr_dim1 + 1], &c__1);
                i__1 = ii + is * vr_dim1;
                remax = 1.f / ((r__1 = vr[i__1].r, f2c_abs(r__1)) + (r__2 = vr[ii + is * vr_dim1].i, f2c_abs(r__2)));
                csscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);
                i__1 = *n;
                for (k = ki + 1;
                        k <= i__1;
                        ++k)
                {
                    i__3 = k + is * vr_dim1;
                    vr[i__3].r = 0.f;
                    vr[i__3].i = 0.f; // , expr subst
                    /* L60: */
                }
            }
            else if (nb == 1)
            {
                /* ------------------------------ */
                /* version 1: back-transform each vector with GEMV, Q*x. */
                if (ki > 1)
                {
                    i__1 = ki - 1;
                    q__1.r = scale;
                    q__1.i = 0.f; // , expr subst
                    cgemv_("N", n, &i__1, &c_b2, &vr[vr_offset], ldvr, &work[ iv * *n + 1], &c__1, &q__1, &vr[ki * vr_dim1 + 1], &c__1);
                }
                ii = icamax_(n, &vr[ki * vr_dim1 + 1], &c__1);
                i__1 = ii + ki * vr_dim1;
                remax = 1.f / ((r__1 = vr[i__1].r, f2c_abs(r__1)) + (r__2 = vr[ii + ki * vr_dim1].i, f2c_abs(r__2)));
                csscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
            }
            else
            {
                /* ------------------------------ */
                /* version 2: back-transform block of vectors with GEMM */
                /* zero out below vector */
                i__1 = *n;
                for (k = ki + 1;
                        k <= i__1;
                        ++k)
                {
                    i__3 = k + iv * *n;
                    work[i__3].r = 0.f;
                    work[i__3].i = 0.f; // , expr subst
                }
                /* Columns IV:NB of work are valid vectors. */
                /* When the number of vectors stored reaches NB, */
                /* or if this was last vector, do the GEMM */
                if (iv == 1 || ki == 1)
                {
                    i__1 = nb - iv + 1;
                    i__3 = ki + nb - iv;
                    cgemm_("N", "N", n, &i__1, &i__3, &c_b2, &vr[vr_offset], ldvr, &work[iv * *n + 1], n, &c_b1, &work[(nb + iv) * *n + 1], n);
                    /* normalize vectors */
                    i__1 = nb;
                    for (k = iv;
                            k <= i__1;
                            ++k)
                    {
                        ii = icamax_(n, &work[(nb + k) * *n + 1], &c__1);
                        i__3 = ii + (nb + k) * *n;
                        remax = 1.f / ((r__1 = work[i__3].r, f2c_abs(r__1)) + ( r__2 = work[ii + (nb + k) * *n].i, f2c_abs( r__2)));
                        csscal_(n, &remax, &work[(nb + k) * *n + 1], &c__1);
                    }
                    i__1 = nb - iv + 1;
                    clacpy_("F", n, &i__1, &work[(nb + iv) * *n + 1], n, &vr[ ki * vr_dim1 + 1], ldvr);
                    iv = nb;
                }
                else
                {
                    --iv;
                }
            }
            /* Restore the original diagonal elements of T. */
            i__1 = ki - 1;
            for (k = 1;
                    k <= i__1;
                    ++k)
            {
                i__3 = k + k * t_dim1;
                i__4 = k;
                t[i__3].r = work[i__4].r;
                t[i__3].i = work[i__4].i; // , expr subst
                /* L70: */
            }
            --is;
L80:
            ;
        }
    }
    if (leftv)
    {
        /* ============================================================ */
        /* Compute left eigenvectors. */
        /* IV is index of column in current block. */
        /* Non-blocked version always uses IV=1;
        */
        /* blocked version starts with IV=1, goes up to NB. */
        /* (Note the "0-th" column is used to store the original diagonal.) */
        iv = 1;
        is = 1;
        i__1 = *n;
        for (ki = 1;
                ki <= i__1;
                ++ki)
        {
            if (somev)
            {
                if (! select[ki])
                {
                    goto L130;
                }
            }
            /* Computing MAX */
            i__3 = ki + ki * t_dim1;
            r__3 = ulp * ((r__1 = t[i__3].r, f2c_abs(r__1)) + (r__2 = t[ ki + ki * t_dim1].i, f2c_abs(r__2)));
            smin = fla_max(r__3,smlnum);
            /* -------------------------------------------------------- */
            /* Complex left eigenvector */
            i__3 = ki + iv * *n;
            work[i__3].r = 1.f;
            work[i__3].i = 0.f; // , expr subst
            /* Form right-hand side. */
            i__3 = *n;
            for (k = ki + 1;
                    k <= i__3;
                    ++k)
            {
                i__4 = k + iv * *n;
                q__2.r = t[ki + k * t_dim1].r;
                q__2.i = -t[ki + k * t_dim1].i;
                q__1.r = -q__2.r;
                q__1.i = -q__2.i; // , expr subst
                work[i__4].r = q__1.r;
                work[i__4].i = q__1.i; // , expr subst
                /* L90: */
            }
            /* Solve conjugate-transposed triangular system: */
            /* [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE*WORK. */
            i__3 = *n;
            for (k = ki + 1;
                    k <= i__3;
                    ++k)
            {
                i__4 = k + k * t_dim1;
                i__5 = k + k * t_dim1;
                i__6 = ki + ki * t_dim1;
                q__1.r = t[i__5].r - t[i__6].r;
                q__1.i = t[i__5].i - t[i__6] .i; // , expr subst
                t[i__4].r = q__1.r;
                t[i__4].i = q__1.i; // , expr subst
                i__4 = k + k * t_dim1;
                if ((r__1 = t[i__4].r, f2c_abs(r__1)) + (r__2 = t[k + k * t_dim1].i, f2c_abs(r__2)) < smin)
                {
                    i__5 = k + k * t_dim1;
                    t[i__5].r = smin;
                    t[i__5].i = 0.f; // , expr subst
                }
                /* L100: */
            }
            if (ki < *n)
            {
                i__3 = *n - ki;
                clatrs_("Upper", "Conjugate transpose", "Non-unit", "Y", & i__3, &t[ki + 1 + (ki + 1) * t_dim1], ldt, &work[ki + 1 + iv * *n], &scale, &rwork[1], info);
                i__3 = ki + iv * *n;
                work[i__3].r = scale;
                work[i__3].i = 0.f; // , expr subst
            }
            /* Copy the vector x or Q*x to VL and normalize. */
            if (! over)
            {
                /* ------------------------------ */
                /* no back-transform: copy x to VL and normalize. */
                i__3 = *n - ki + 1;
                ccopy_(&i__3, &work[ki + iv * *n], &c__1, &vl[ki + is * vl_dim1], &c__1);
                i__3 = *n - ki + 1;
                ii = icamax_(&i__3, &vl[ki + is * vl_dim1], &c__1) + ki - 1;
                i__3 = ii + is * vl_dim1;
                remax = 1.f / ((r__1 = vl[i__3].r, f2c_abs(r__1)) + (r__2 = vl[ii + is * vl_dim1].i, f2c_abs(r__2)));
                i__3 = *n - ki + 1;
                csscal_(&i__3, &remax, &vl[ki + is * vl_dim1], &c__1);
                i__3 = ki - 1;
                for (k = 1;
                        k <= i__3;
                        ++k)
                {
                    i__4 = k + is * vl_dim1;
                    vl[i__4].r = 0.f;
                    vl[i__4].i = 0.f; // , expr subst
                    /* L110: */
                }
            }
            else if (nb == 1)
            {
                /* ------------------------------ */
                /* version 1: back-transform each vector with GEMV, Q*x. */
                if (ki < *n)
                {
                    i__3 = *n - ki;
                    q__1.r = scale;
                    q__1.i = 0.f; // , expr subst
                    cgemv_("N", n, &i__3, &c_b2, &vl[(ki + 1) * vl_dim1 + 1], ldvl, &work[ki + 1 + iv * *n], &c__1, &q__1, &vl[ ki * vl_dim1 + 1], &c__1);
                }
                ii = icamax_(n, &vl[ki * vl_dim1 + 1], &c__1);
                i__3 = ii + ki * vl_dim1;
                remax = 1.f / ((r__1 = vl[i__3].r, f2c_abs(r__1)) + (r__2 = vl[ii + ki * vl_dim1].i, f2c_abs(r__2)));
                csscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
            }
            else
            {
                /* ------------------------------ */
                /* version 2: back-transform block of vectors with GEMM */
                /* zero out above vector */
                /* could go from KI-NV+1 to KI-1 */
                i__3 = ki - 1;
                for (k = 1;
                        k <= i__3;
                        ++k)
                {
                    i__4 = k + iv * *n;
                    work[i__4].r = 0.f;
                    work[i__4].i = 0.f; // , expr subst
                }
                /* Columns 1:IV of work are valid vectors. */
                /* When the number of vectors stored reaches NB, */
                /* or if this was last vector, do the GEMM */
                if (iv == nb || ki == *n)
                {
                    i__3 = *n - ki + iv;
                    cgemm_("N", "N", n, &iv, &i__3, &c_b2, &vl[(ki - iv + 1) * vl_dim1 + 1], ldvl, &work[ki - iv + 1 + *n], n, & c_b1, &work[(nb + 1) * *n + 1], n);
                    /* normalize vectors */
                    i__3 = iv;
                    for (k = 1;
                            k <= i__3;
                            ++k)
                    {
                        ii = icamax_(n, &work[(nb + k) * *n + 1], &c__1);
                        i__4 = ii + (nb + k) * *n;
                        remax = 1.f / ((r__1 = work[i__4].r, f2c_abs(r__1)) + ( r__2 = work[ii + (nb + k) * *n].i, f2c_abs( r__2)));
                        csscal_(n, &remax, &work[(nb + k) * *n + 1], &c__1);
                    }
                    clacpy_("F", n, &iv, &work[(nb + 1) * *n + 1], n, &vl[(ki - iv + 1) * vl_dim1 + 1], ldvl);
                    iv = 1;
                }
                else
                {
                    ++iv;
                }
            }
            /* Restore the original diagonal elements of T. */
            i__3 = *n;
            for (k = ki + 1;
                    k <= i__3;
                    ++k)
            {
                i__4 = k + k * t_dim1;
                i__5 = k;
                t[i__4].r = work[i__5].r;
                t[i__4].i = work[i__5].i; // , expr subst
                /* L120: */
            }
            ++is;
L130:
            ;
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CTREVC3 */
}
/* ctrevc3_ */

