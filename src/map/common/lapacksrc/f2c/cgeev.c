/* ../netlib/cgeev.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
/* > \brief <b> CGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matr ices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGEEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeev.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeev.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeev.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, */
/* WORK, LWORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVL, JOBVR */
/* INTEGER INFO, LDA, LDVL, LDVR, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* REAL RWORK( * ) */
/* COMPLEX A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), */
/* $ W( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEEV computes for an N-by-N complex nonsymmetric matrix A, the */
/* > eigenvalues and, optionally, the left and/or right eigenvectors. */
/* > */
/* > The right eigenvector v(j) of A satisfies */
/* > A * v(j) = lambda(j) * v(j) */
/* > where lambda(j) is its eigenvalue. */
/* > The left eigenvector u(j) of A satisfies */
/* > u(j)**H * A = lambda(j) * u(j)**H */
/* > where u(j)**H denotes the conjugate transpose of u(j). */
/* > */
/* > The computed eigenvectors are normalized to have Euclidean norm */
/* > equal to 1 and largest component real. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVL */
/* > \verbatim */
/* > JOBVL is CHARACTER*1 */
/* > = 'N': left eigenvectors of A are not computed;
*/
/* > = 'V': left eigenvectors of are computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVR */
/* > \verbatim */
/* > JOBVR is CHARACTER*1 */
/* > = 'N': right eigenvectors of A are not computed;
*/
/* > = 'V': right eigenvectors of A are computed. */
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
/* > On entry, the N-by-N matrix A. */
/* > On exit, A has been overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is COMPLEX array, dimension (N) */
/* > W contains the computed eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* > VL is COMPLEX array, dimension (LDVL,N) */
/* > If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/* > after another in the columns of VL, in the same order */
/* > as their eigenvalues. */
/* > If JOBVL = 'N', VL is not referenced. */
/* > u(j) = VL(:,j), the j-th column of VL. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* > LDVL is INTEGER */
/* > The leading dimension of the array VL. LDVL >= 1;
if */
/* > JOBVL = 'V', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VR */
/* > \verbatim */
/* > VR is COMPLEX array, dimension (LDVR,N) */
/* > If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/* > after another in the columns of VR, in the same order */
/* > as their eigenvalues. */
/* > If JOBVR = 'N', VR is not referenced. */
/* > v(j) = VR(:,j), the j-th column of VR. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* > LDVR is INTEGER */
/* > The leading dimension of the array VR. LDVR >= 1;
if */
/* > JOBVR = 'V', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= max(1,2*N). */
/* > For good performance, LWORK must generally be larger. */
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
/* > RWORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = i, the QR algorithm failed to compute all the */
/* > eigenvalues, and no eigenvectors have been computed;
*/
/* > elements and i+1:N of W contain eigenvalues which have */
/* > converged. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexGEeigen */
/* ===================================================================== */
/* Subroutine */
int cgeev_(char *jobvl, char *jobvr, integer *n, complex *a, integer *lda, complex *w, complex *vl, integer *ldvl, complex *vr, integer *ldvr, complex *work, integer *lwork, real *rwork, integer * info)
{
    /* System generated locals */
    integer a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1, q__2;
    /* Builtin functions */
    double sqrt(doublereal), r_imag(complex *);
    void r_cnjg(complex *, complex *);
    /* Local variables */
    integer i__, k, ihi;
    real scl;
    integer ilo;
    real dum[1], eps;
    complex tmp;
    integer ibal;
    char side[1];
    real anrm;
    integer ierr, itau, iwrk, nout;
    extern /* Subroutine */
    int cscal_(integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *);
    extern real scnrm2_(integer *, complex *, integer *);
    extern /* Subroutine */
    int cgebak_(char *, char *, integer *, integer *, integer *, real *, integer *, complex *, integer *, integer *), cgebal_(char *, integer *, complex *, integer *, integer *, integer *, real *, integer *), slabad_(real *, real *);
    logical scalea;
    extern real clange_(char *, integer *, integer *, complex *, integer *, real *);
    real cscale;
    extern /* Subroutine */
    int cgehrd_(integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), clascl_(char *, integer *, integer *, real *, real *, integer *, integer *, complex *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int csscal_(integer *, real *, complex *, integer *), clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    logical select[1];
    real bignum;
    extern integer isamax_(integer *, real *, integer *);
    extern /* Subroutine */
    int chseqr_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *), ctrevc_(char *, char *, logical *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, integer *, integer *, complex *, real *, integer *), cunghr_(integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *);
    integer minwrk, maxwrk;
    logical wantvl;
    real smlnum;
    integer hswork, irwork;
    logical lquery, wantvr;
    /* -- LAPACK driver routine (version 3.4.0) -- */
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    --rwork;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    wantvl = lsame_(jobvl, "V");
    wantvr = lsame_(jobvr, "V");
    if (! wantvl && ! lsame_(jobvl, "N"))
    {
        *info = -1;
    }
    else if (! wantvr && ! lsame_(jobvr, "N"))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    else if (*ldvl < 1 || wantvl && *ldvl < *n)
    {
        *info = -8;
    }
    else if (*ldvr < 1 || wantvr && *ldvr < *n)
    {
        *info = -10;
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* CWorkspace refers to complex workspace, and RWorkspace to real */
    /* workspace. NB refers to the optimal block size for the */
    /* immediately following subroutine, as returned by ILAENV. */
    /* HSWORK refers to the workspace preferred by CHSEQR, as */
    /* calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
    /* the worst case.) */
    if (*info == 0)
    {
        if (*n == 0)
        {
            minwrk = 1;
            maxwrk = 1;
        }
        else
        {
            maxwrk = *n + *n * ilaenv_(&c__1, "CGEHRD", " ", n, &c__1, n, & c__0);
            minwrk = *n << 1;
            if (wantvl)
            {
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n + (*n - 1) * ilaenv_(&c__1, "CUNGHR", " ", n, &c__1, n, &c_n1); // , expr subst
                maxwrk = max(i__1,i__2);
                chseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vl[ vl_offset], ldvl, &work[1], &c_n1, info);
            }
            else if (wantvr)
            {
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n + (*n - 1) * ilaenv_(&c__1, "CUNGHR", " ", n, &c__1, n, &c_n1); // , expr subst
                maxwrk = max(i__1,i__2);
                chseqr_("S", "V", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[ vr_offset], ldvr, &work[1], &c_n1, info);
            }
            else
            {
                chseqr_("E", "N", n, &c__1, n, &a[a_offset], lda, &w[1], &vr[ vr_offset], ldvr, &work[1], &c_n1, info);
            }
            hswork = work[1].r;
            /* Computing MAX */
            i__1 = max(maxwrk,hswork);
            maxwrk = max(i__1,minwrk);
        }
        work[1].r = (real) maxwrk;
        work[1].i = 0.f; // , expr subst
        if (*lwork < minwrk && ! lquery)
        {
            *info = -12;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGEEV ", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        return 0;
    }
    /* Get machine constants */
    eps = slamch_("P");
    smlnum = slamch_("S");
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = 1.f / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = clange_("M", n, n, &a[a_offset], lda, dum);
    scalea = FALSE_;
    if (anrm > 0.f && anrm < smlnum)
    {
        scalea = TRUE_;
        cscale = smlnum;
    }
    else if (anrm > bignum)
    {
        scalea = TRUE_;
        cscale = bignum;
    }
    if (scalea)
    {
        clascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, & ierr);
    }
    /* Balance the matrix */
    /* (CWorkspace: none) */
    /* (RWorkspace: need N) */
    ibal = 1;
    cgebal_("B", n, &a[a_offset], lda, &ilo, &ihi, &rwork[ibal], &ierr);
    /* Reduce to upper Hessenberg form */
    /* (CWorkspace: need 2*N, prefer N+N*NB) */
    /* (RWorkspace: none) */
    itau = 1;
    iwrk = itau + *n;
    i__1 = *lwork - iwrk + 1;
    cgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &ierr);
    if (wantvl)
    {
        /* Want left eigenvectors */
        /* Copy Householder vectors to VL */
        *(unsigned char *)side = 'L';
        clacpy_("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl) ;
        /* Generate unitary matrix in VL */
        /* (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
        /* (RWorkspace: none) */
        i__1 = *lwork - iwrk + 1;
        cunghr_(n, &ilo, &ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk], &i__1, &ierr);
        /* Perform QR iteration, accumulating Schur vectors in VL */
        /* (CWorkspace: need 1, prefer HSWORK (see comments) ) */
        /* (RWorkspace: none) */
        iwrk = itau;
        i__1 = *lwork - iwrk + 1;
        chseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vl[ vl_offset], ldvl, &work[iwrk], &i__1, info);
        if (wantvr)
        {
            /* Want left and right eigenvectors */
            /* Copy Schur vectors to VR */
            *(unsigned char *)side = 'B';
            clacpy_("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr);
        }
    }
    else if (wantvr)
    {
        /* Want right eigenvectors */
        /* Copy Householder vectors to VR */
        *(unsigned char *)side = 'R';
        clacpy_("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr) ;
        /* Generate unitary matrix in VR */
        /* (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
        /* (RWorkspace: none) */
        i__1 = *lwork - iwrk + 1;
        cunghr_(n, &ilo, &ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk], &i__1, &ierr);
        /* Perform QR iteration, accumulating Schur vectors in VR */
        /* (CWorkspace: need 1, prefer HSWORK (see comments) ) */
        /* (RWorkspace: none) */
        iwrk = itau;
        i__1 = *lwork - iwrk + 1;
        chseqr_("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vr[ vr_offset], ldvr, &work[iwrk], &i__1, info);
    }
    else
    {
        /* Compute eigenvalues only */
        /* (CWorkspace: need 1, prefer HSWORK (see comments) ) */
        /* (RWorkspace: none) */
        iwrk = itau;
        i__1 = *lwork - iwrk + 1;
        chseqr_("E", "N", n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vr[ vr_offset], ldvr, &work[iwrk], &i__1, info);
    }
    /* If INFO > 0 from CHSEQR, then quit */
    if (*info > 0)
    {
        goto L50;
    }
    if (wantvl || wantvr)
    {
        /* Compute left and/or right eigenvectors */
        /* (CWorkspace: need 2*N) */
        /* (RWorkspace: need 2*N) */
        irwork = ibal + *n;
        ctrevc_(side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &rwork[irwork], &ierr);
    }
    if (wantvl)
    {
        /* Undo balancing of left eigenvectors */
        /* (CWorkspace: none) */
        /* (RWorkspace: need N) */
        cgebak_("B", "L", n, &ilo, &ihi, &rwork[ibal], n, &vl[vl_offset], ldvl, &ierr);
        /* Normalize left eigenvectors and make largest component real */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            scl = 1.f / scnrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
            csscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
            i__2 = *n;
            for (k = 1;
                    k <= i__2;
                    ++k)
            {
                i__3 = k + i__ * vl_dim1;
                /* Computing 2nd power */
                r__1 = vl[i__3].r;
                /* Computing 2nd power */
                r__2 = r_imag(&vl[k + i__ * vl_dim1]);
                rwork[irwork + k - 1] = r__1 * r__1 + r__2 * r__2;
                /* L10: */
            }
            k = isamax_(n, &rwork[irwork], &c__1);
            r_cnjg(&q__2, &vl[k + i__ * vl_dim1]);
            r__1 = sqrt(rwork[irwork + k - 1]);
            q__1.r = q__2.r / r__1;
            q__1.i = q__2.i / r__1; // , expr subst
            tmp.r = q__1.r;
            tmp.i = q__1.i; // , expr subst
            cscal_(n, &tmp, &vl[i__ * vl_dim1 + 1], &c__1);
            i__2 = k + i__ * vl_dim1;
            i__3 = k + i__ * vl_dim1;
            r__1 = vl[i__3].r;
            q__1.r = r__1;
            q__1.i = 0.f; // , expr subst
            vl[i__2].r = q__1.r;
            vl[i__2].i = q__1.i; // , expr subst
            /* L20: */
        }
    }
    if (wantvr)
    {
        /* Undo balancing of right eigenvectors */
        /* (CWorkspace: none) */
        /* (RWorkspace: need N) */
        cgebak_("B", "R", n, &ilo, &ihi, &rwork[ibal], n, &vr[vr_offset], ldvr, &ierr);
        /* Normalize right eigenvectors and make largest component real */
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            scl = 1.f / scnrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
            csscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
            i__2 = *n;
            for (k = 1;
                    k <= i__2;
                    ++k)
            {
                i__3 = k + i__ * vr_dim1;
                /* Computing 2nd power */
                r__1 = vr[i__3].r;
                /* Computing 2nd power */
                r__2 = r_imag(&vr[k + i__ * vr_dim1]);
                rwork[irwork + k - 1] = r__1 * r__1 + r__2 * r__2;
                /* L30: */
            }
            k = isamax_(n, &rwork[irwork], &c__1);
            r_cnjg(&q__2, &vr[k + i__ * vr_dim1]);
            r__1 = sqrt(rwork[irwork + k - 1]);
            q__1.r = q__2.r / r__1;
            q__1.i = q__2.i / r__1; // , expr subst
            tmp.r = q__1.r;
            tmp.i = q__1.i; // , expr subst
            cscal_(n, &tmp, &vr[i__ * vr_dim1 + 1], &c__1);
            i__2 = k + i__ * vr_dim1;
            i__3 = k + i__ * vr_dim1;
            r__1 = vr[i__3].r;
            q__1.r = r__1;
            q__1.i = 0.f; // , expr subst
            vr[i__2].r = q__1.r;
            vr[i__2].i = q__1.i; // , expr subst
            /* L40: */
        }
    }
    /* Undo scaling if necessary */
L50:
    if (scalea)
    {
        i__1 = *n - *info;
        /* Computing MAX */
        i__3 = *n - *info;
        i__2 = max(i__3,1);
        clascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[*info + 1] , &i__2, &ierr);
        if (*info > 0)
        {
            i__1 = ilo - 1;
            clascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &w[1], n, &ierr);
        }
    }
    work[1].r = (real) maxwrk;
    work[1].i = 0.f; // , expr subst
    return 0;
    /* End of CGEEV */
}
/* cgeev_ */
