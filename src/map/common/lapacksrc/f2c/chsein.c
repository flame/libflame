/* ../netlib/chsein.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static logical c_false = FALSE_;
static logical c_true = TRUE_;
/* > \brief \b CHSEIN */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHSEIN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chsein. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chsein. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chsein. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL, */
/* LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL, */
/* IFAILR, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER EIGSRC, INITV, SIDE */
/* INTEGER INFO, LDH, LDVL, LDVR, M, MM, N */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL SELECT( * ) */
/* INTEGER IFAILL( * ), IFAILR( * ) */
/* REAL RWORK( * ) */
/* COMPLEX H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ), */
/* $ W( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHSEIN uses inverse iteration to find specified right and/or left */
/* > eigenvectors of a complex upper Hessenberg matrix H. */
/* > */
/* > The right eigenvector x and the left eigenvector y of the matrix H */
/* > corresponding to an eigenvalue w are defined by: */
/* > */
/* > H * x = w * x, y**h * H = w * y**h */
/* > */
/* > where y**h denotes the conjugate transpose of the vector y. */
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
/* > \param[in] EIGSRC */
/* > \verbatim */
/* > EIGSRC is CHARACTER*1 */
/* > Specifies the source of eigenvalues supplied in W: */
/* > = 'Q': the eigenvalues were found using CHSEQR;
thus, if */
/* > H has zero subdiagonal elements, and so is */
/* > block-triangular, then the j-th eigenvalue can be */
/* > assumed to be an eigenvalue of the block containing */
/* > the j-th row/column. This property allows CHSEIN to */
/* > perform inverse iteration on just one diagonal block. */
/* > = 'N': no assumptions are made on the correspondence */
/* > between eigenvalues and diagonal blocks. In this */
/* > case, CHSEIN must always perform inverse iteration */
/* > using the whole matrix H. */
/* > \endverbatim */
/* > */
/* > \param[in] INITV */
/* > \verbatim */
/* > INITV is CHARACTER*1 */
/* > = 'N': no initial vectors are supplied;
*/
/* > = 'U': user-supplied initial vectors are stored in the arrays */
/* > VL and/or VR. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* > SELECT is LOGICAL array, dimension (N) */
/* > Specifies the eigenvectors to be computed. To select the */
/* > eigenvector corresponding to the eigenvalue W(j), */
/* > SELECT(j) must be set to .TRUE.. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix H. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* > H is COMPLEX array, dimension (LDH,N) */
/* > The upper Hessenberg matrix H. */
/* > If a NaN is detected in H, the routine will return with INFO=-6. */
/* > \endverbatim */
/* > */
/* > \param[in] LDH */
/* > \verbatim */
/* > LDH is INTEGER */
/* > The leading dimension of the array H. LDH >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] W */
/* > \verbatim */
/* > W is COMPLEX array, dimension (N) */
/* > On entry, the eigenvalues of H. */
/* > On exit, the real parts of W may have been altered since */
/* > close eigenvalues are perturbed slightly in searching for */
/* > independent eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VL */
/* > \verbatim */
/* > VL is COMPLEX array, dimension (LDVL,MM) */
/* > On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must */
/* > contain starting vectors for the inverse iteration for the */
/* > left eigenvectors;
the starting vector for each eigenvector */
/* > must be in the same column in which the eigenvector will be */
/* > stored. */
/* > On exit, if SIDE = 'L' or 'B', the left eigenvectors */
/* > specified by SELECT will be stored consecutively in the */
/* > columns of VL, in the same order as their eigenvalues. */
/* > If SIDE = 'R', VL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* > LDVL is INTEGER */
/* > The leading dimension of the array VL. */
/* > LDVL >= max(1,N) if SIDE = 'L' or 'B';
LDVL >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VR */
/* > \verbatim */
/* > VR is COMPLEX array, dimension (LDVR,MM) */
/* > On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must */
/* > contain starting vectors for the inverse iteration for the */
/* > right eigenvectors;
the starting vector for each eigenvector */
/* > must be in the same column in which the eigenvector will be */
/* > stored. */
/* > On exit, if SIDE = 'R' or 'B', the right eigenvectors */
/* > specified by SELECT will be stored consecutively in the */
/* > columns of VR, in the same order as their eigenvalues. */
/* > If SIDE = 'L', VR is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* > LDVR is INTEGER */
/* > The leading dimension of the array VR. */
/* > LDVR >= max(1,N) if SIDE = 'R' or 'B';
LDVR >= 1 otherwise. */
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
/* > The number of columns in the arrays VL and/or VR required to */
/* > store the eigenvectors (= the number of .TRUE. elements in */
/* > SELECT). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] IFAILL */
/* > \verbatim */
/* > IFAILL is INTEGER array, dimension (MM) */
/* > If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left */
/* > eigenvector in the i-th column of VL (corresponding to the */
/* > eigenvalue w(j)) failed to converge;
IFAILL(i) = 0 if the */
/* > eigenvector converged satisfactorily. */
/* > If SIDE = 'R', IFAILL is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] IFAILR */
/* > \verbatim */
/* > IFAILR is INTEGER array, dimension (MM) */
/* > If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right */
/* > eigenvector in the i-th column of VR (corresponding to the */
/* > eigenvalue w(j)) failed to converge;
IFAILR(i) = 0 if the */
/* > eigenvector converged satisfactorily. */
/* > If SIDE = 'L', IFAILR is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, i is the number of eigenvectors which */
/* > failed to converge;
see IFAILL and IFAILR for further */
/* > details. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup complexOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Each eigenvector is normalized so that the element of largest */
/* > magnitude has magnitude 1;
here the magnitude of a complex number */
/* > (x,y) is taken to be |x|+|y|. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int chsein_(char *side, char *eigsrc, char *initv, logical * select, integer *n, complex *h__, integer *ldh, complex *w, complex * vl, integer *ldvl, complex *vr, integer *ldvr, integer *mm, integer * m, complex *work, real *rwork, integer *ifaill, integer *ifailr, integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1, q__2;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    integer i__, k, kl, kr, ks;
    complex wk;
    integer kln;
    real ulp, eps3, unfl;
    extern logical lsame_(char *, char *);
    integer iinfo;
    logical leftv, bothv;
    real hnorm;
    extern /* Subroutine */
    int claein_(logical *, logical *, integer *, complex *, integer *, complex *, complex *, complex *, integer *, real *, real *, real *, integer *);
    extern real slamch_(char *), clanhs_(char *, integer *, complex *, integer *, real *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern logical sisnan_(real *);
    logical noinit;
    integer ldwork;
    logical rightv, fromqr;
    real smlnum;
    /* -- LAPACK computational routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
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
    /* Decode and test the input parameters. */
    /* Parameter adjustments */
    --select;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --w;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    --rwork;
    --ifaill;
    --ifailr;
    /* Function Body */
    bothv = lsame_(side, "B");
    rightv = lsame_(side, "R") || bothv;
    leftv = lsame_(side, "L") || bothv;
    fromqr = lsame_(eigsrc, "Q");
    noinit = lsame_(initv, "N");
    /* Set M to the number of columns required to store the selected */
    /* eigenvectors. */
    *m = 0;
    i__1 = *n;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        if (select[k])
        {
            ++(*m);
        }
        /* L10: */
    }
    *info = 0;
    if (! rightv && ! leftv)
    {
        *info = -1;
    }
    else if (! fromqr && ! lsame_(eigsrc, "N"))
    {
        *info = -2;
    }
    else if (! noinit && ! lsame_(initv, "U"))
    {
        *info = -3;
    }
    else if (*n < 0)
    {
        *info = -5;
    }
    else if (*ldh < max(1,*n))
    {
        *info = -7;
    }
    else if (*ldvl < 1 || leftv && *ldvl < *n)
    {
        *info = -10;
    }
    else if (*ldvr < 1 || rightv && *ldvr < *n)
    {
        *info = -12;
    }
    else if (*mm < *m)
    {
        *info = -13;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHSEIN", &i__1);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0)
    {
        return 0;
    }
    /* Set machine-dependent constants. */
    unfl = slamch_("Safe minimum");
    ulp = slamch_("Precision");
    smlnum = unfl * (*n / ulp);
    ldwork = *n;
    kl = 1;
    kln = 0;
    if (fromqr)
    {
        kr = 0;
    }
    else
    {
        kr = *n;
    }
    ks = 1;
    i__1 = *n;
    for (k = 1;
            k <= i__1;
            ++k)
    {
        if (select[k])
        {
            /* Compute eigenvector(s) corresponding to W(K). */
            if (fromqr)
            {
                /* If affiliation of eigenvalues is known, check whether */
                /* the matrix splits. */
                /* Determine KL and KR such that 1 <= KL <= K <= KR <= N */
                /* and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or */
                /* KR = N). */
                /* Then inverse iteration can be performed with the */
                /* submatrix H(KL:N,KL:N) for a left eigenvector, and with */
                /* the submatrix H(1:KR,1:KR) for a right eigenvector. */
                i__2 = kl + 1;
                for (i__ = k;
                        i__ >= i__2;
                        --i__)
                {
                    i__3 = i__ + (i__ - 1) * h_dim1;
                    if (h__[i__3].r == 0.f && h__[i__3].i == 0.f)
                    {
                        goto L30;
                    }
                    /* L20: */
                }
L30:
                kl = i__;
                if (k > kr)
                {
                    i__2 = *n - 1;
                    for (i__ = k;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + 1 + i__ * h_dim1;
                        if (h__[i__3].r == 0.f && h__[i__3].i == 0.f)
                        {
                            goto L50;
                        }
                        /* L40: */
                    }
L50:
                    kr = i__;
                }
            }
            if (kl != kln)
            {
                kln = kl;
                /* Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it */
                /* has not ben computed before. */
                i__2 = kr - kl + 1;
                hnorm = clanhs_("I", &i__2, &h__[kl + kl * h_dim1], ldh, & rwork[1]);
                if (sisnan_(&hnorm))
                {
                    *info = -6;
                    return 0;
                }
                else if (hnorm > 0.f)
                {
                    eps3 = hnorm * ulp;
                }
                else
                {
                    eps3 = smlnum;
                }
            }
            /* Perturb eigenvalue if it is close to any previous */
            /* selected eigenvalues affiliated to the submatrix */
            /* H(KL:KR,KL:KR). Close roots are modified by EPS3. */
            i__2 = k;
            wk.r = w[i__2].r;
            wk.i = w[i__2].i; // , expr subst
L60:
            i__2 = kl;
            for (i__ = k - 1;
                    i__ >= i__2;
                    --i__)
            {
                i__3 = i__;
                q__2.r = w[i__3].r - wk.r;
                q__2.i = w[i__3].i - wk.i; // , expr subst
                q__1.r = q__2.r;
                q__1.i = q__2.i; // , expr subst
                if (select[i__] && (r__1 = q__1.r, f2c_abs(r__1)) + (r__2 = r_imag(&q__1), f2c_abs(r__2)) < eps3)
                {
                    q__1.r = wk.r + eps3;
                    q__1.i = wk.i; // , expr subst
                    wk.r = q__1.r;
                    wk.i = q__1.i; // , expr subst
                    goto L60;
                }
                /* L70: */
            }
            i__2 = k;
            w[i__2].r = wk.r;
            w[i__2].i = wk.i; // , expr subst
            if (leftv)
            {
                /* Compute left eigenvector. */
                i__2 = *n - kl + 1;
                claein_(&c_false, &noinit, &i__2, &h__[kl + kl * h_dim1], ldh, &wk, &vl[kl + ks * vl_dim1], &work[1], &ldwork, & rwork[1], &eps3, &smlnum, &iinfo);
                if (iinfo > 0)
                {
                    ++(*info);
                    ifaill[ks] = k;
                }
                else
                {
                    ifaill[ks] = 0;
                }
                i__2 = kl - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + ks * vl_dim1;
                    vl[i__3].r = 0.f;
                    vl[i__3].i = 0.f; // , expr subst
                    /* L80: */
                }
            }
            if (rightv)
            {
                /* Compute right eigenvector. */
                claein_(&c_true, &noinit, &kr, &h__[h_offset], ldh, &wk, &vr[ ks * vr_dim1 + 1], &work[1], &ldwork, &rwork[1], & eps3, &smlnum, &iinfo);
                if (iinfo > 0)
                {
                    ++(*info);
                    ifailr[ks] = k;
                }
                else
                {
                    ifailr[ks] = 0;
                }
                i__2 = *n;
                for (i__ = kr + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + ks * vr_dim1;
                    vr[i__3].r = 0.f;
                    vr[i__3].i = 0.f; // , expr subst
                    /* L90: */
                }
            }
            ++ks;
        }
        /* L100: */
    }
    return 0;
    /* End of CHSEIN */
}
/* chsein_ */
