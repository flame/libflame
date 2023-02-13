/* ../netlib/cgees.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
/* > \brief <b> CGEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors f or GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGEES + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgees.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgees.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgees.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS, */
/* LDVS, WORK, LWORK, RWORK, BWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVS, SORT */
/* INTEGER INFO, LDA, LDVS, LWORK, N, SDIM */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL BWORK( * ) */
/* REAL RWORK( * ) */
/* COMPLEX A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * ) */
/* .. */
/* .. Function Arguments .. */
/* LOGICAL SELECT */
/* EXTERNAL SELECT */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEES computes for an N-by-N complex nonsymmetric matrix A, the */
/* > eigenvalues, the Schur form T, and, optionally, the matrix of Schur */
/* > vectors Z. This gives the Schur factorization A = Z*T*(Z**H). */
/* > */
/* > Optionally, it also orders the eigenvalues on the diagonal of the */
/* > Schur form so that selected eigenvalues are at the top left. */
/* > The leading columns of Z then form an orthonormal basis for the */
/* > invariant subspace corresponding to the selected eigenvalues. */
/* > */
/* > A complex matrix is in Schur form if it is upper triangular. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVS */
/* > \verbatim */
/* > JOBVS is CHARACTER*1 */
/* > = 'N': Schur vectors are not computed;
*/
/* > = 'V': Schur vectors are computed. */
/* > \endverbatim */
/* > */
/* > \param[in] SORT */
/* > \verbatim */
/* > SORT is CHARACTER*1 */
/* > Specifies whether or not to order the eigenvalues on the */
/* > diagonal of the Schur form. */
/* > = 'N': Eigenvalues are not ordered: */
/* > = 'S': Eigenvalues are ordered (see SELECT). */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* > SELECT is a LOGICAL FUNCTION of one COMPLEX argument */
/* > SELECT must be declared EXTERNAL in the calling subroutine. */
/* > If SORT = 'S', SELECT is used to select eigenvalues to order */
/* > to the top left of the Schur form. */
/* > IF SORT = 'N', SELECT is not referenced. */
/* > The eigenvalue W(j) is selected if SELECT(W(j)) is true. */
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
/* > On exit, A has been overwritten by its Schur form T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] SDIM */
/* > \verbatim */
/* > SDIM is INTEGER */
/* > If SORT = 'N', SDIM = 0. */
/* > If SORT = 'S', SDIM = number of eigenvalues for which */
/* > SELECT is true. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is COMPLEX array, dimension (N) */
/* > W contains the computed eigenvalues, in the same order that */
/* > they appear on the diagonal of the output Schur form T. */
/* > \endverbatim */
/* > */
/* > \param[out] VS */
/* > \verbatim */
/* > VS is COMPLEX array, dimension (LDVS,N) */
/* > If JOBVS = 'V', VS contains the unitary matrix Z of Schur */
/* > vectors. */
/* > If JOBVS = 'N', VS is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVS */
/* > \verbatim */
/* > LDVS is INTEGER */
/* > The leading dimension of the array VS. LDVS >= 1;
if */
/* > JOBVS = 'V', LDVS >= N. */
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
/* > RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BWORK */
/* > \verbatim */
/* > BWORK is LOGICAL array, dimension (N) */
/* > Not referenced if SORT = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = i, and i is */
/* > <= N: the QR algorithm failed to compute all the */
/* > eigenvalues;
elements 1:ILO-1 and i+1:N of W */
/* > contain those eigenvalues which have converged;
*/
/* > if JOBVS = 'V', VS contains the matrix which */
/* > reduces A to its partially converged Schur form. */
/* > = N+1: the eigenvalues could not be reordered because */
/* > some eigenvalues were too close to separate (the */
/* > problem is very ill-conditioned);
*/
/* > = N+2: after reordering, roundoff changed values of */
/* > some complex eigenvalues so that leading */
/* > eigenvalues in the Schur form no longer satisfy */
/* > SELECT = .TRUE.. This could also be caused by */
/* > underflow due to scaling. */
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
int cgees_(char *jobvs, char *sort, L_fp select, integer *n, complex *a, integer *lda, integer *sdim, complex *w, complex *vs, integer *ldvs, complex *work, integer *lwork, real *rwork, logical * bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, vs_dim1, vs_offset, i__1, i__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__;
    real s;
    integer ihi, ilo;
    real dum[1], eps, sep;
    integer ibal;
    real anrm;
    integer ierr, itau, iwrk, icond, ieval;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int ccopy_(integer *, complex *, integer *, complex *, integer *), cgebak_(char *, char *, integer *, integer *, integer *, real *, integer *, complex *, integer *, integer *), cgebal_(char *, integer *, complex *, integer *, integer *, integer *, real *, integer *), slabad_(real *, real *);
    logical scalea;
    extern real clange_(char *, integer *, integer *, complex *, integer *, real *);
    real cscale;
    extern /* Subroutine */
    int cgehrd_(integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), clascl_(char *, integer *, integer *, real *, real *, integer *, integer *, complex *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    real bignum;
    extern /* Subroutine */
    int chseqr_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *), cunghr_(integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), ctrsen_(char *, char *, logical *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, real *, real *, complex *, integer *, integer *);
    integer minwrk, maxwrk;
    real smlnum;
    integer hswork;
    logical wantst, lquery, wantvs;
    /* -- LAPACK driver routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* .. Function Arguments .. */
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
    vs_dim1 = *ldvs;
    vs_offset = 1 + vs_dim1;
    vs -= vs_offset;
    --work;
    --rwork;
    --bwork;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    wantvs = lsame_(jobvs, "V");
    wantst = lsame_(sort, "S");
    if (! wantvs && ! lsame_(jobvs, "N"))
    {
        *info = -1;
    }
    else if (! wantst && ! lsame_(sort, "N"))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -4;
    }
    else if (*lda < max(1,*n))
    {
        *info = -6;
    }
    else if (*ldvs < 1 || wantvs && *ldvs < *n)
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
            chseqr_("S", jobvs, n, &c__1, n, &a[a_offset], lda, &w[1], &vs[ vs_offset], ldvs, &work[1], &c_n1, &ieval);
            hswork = work[1].r;
            if (! wantvs)
            {
                maxwrk = max(maxwrk,hswork);
            }
            else
            {
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n + (*n - 1) * ilaenv_(&c__1, "CUNGHR", " ", n, &c__1, n, &c_n1); // , expr subst
                maxwrk = max(i__1,i__2);
                maxwrk = max(maxwrk,hswork);
            }
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
        xerbla_("CGEES ", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        *sdim = 0;
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
    /* Permute the matrix to make it more nearly triangular */
    /* (CWorkspace: none) */
    /* (RWorkspace: need N) */
    ibal = 1;
    cgebal_("P", n, &a[a_offset], lda, &ilo, &ihi, &rwork[ibal], &ierr);
    /* Reduce to upper Hessenberg form */
    /* (CWorkspace: need 2*N, prefer N+N*NB) */
    /* (RWorkspace: none) */
    itau = 1;
    iwrk = *n + itau;
    i__1 = *lwork - iwrk + 1;
    cgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &ierr);
    if (wantvs)
    {
        /* Copy Householder vectors to VS */
        clacpy_("L", n, n, &a[a_offset], lda, &vs[vs_offset], ldvs) ;
        /* Generate unitary matrix in VS */
        /* (CWorkspace: need 2*N-1, prefer N+(N-1)*NB) */
        /* (RWorkspace: none) */
        i__1 = *lwork - iwrk + 1;
        cunghr_(n, &ilo, &ihi, &vs[vs_offset], ldvs, &work[itau], &work[iwrk], &i__1, &ierr);
    }
    *sdim = 0;
    /* Perform QR iteration, accumulating Schur vectors in VS if desired */
    /* (CWorkspace: need 1, prefer HSWORK (see comments) ) */
    /* (RWorkspace: none) */
    iwrk = itau;
    i__1 = *lwork - iwrk + 1;
    chseqr_("S", jobvs, n, &ilo, &ihi, &a[a_offset], lda, &w[1], &vs[ vs_offset], ldvs, &work[iwrk], &i__1, &ieval);
    if (ieval > 0)
    {
        *info = ieval;
    }
    /* Sort eigenvalues if desired */
    if (wantst && *info == 0)
    {
        if (scalea)
        {
            clascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &w[1], n, & ierr);
        }
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            bwork[i__] = (*select)(&w[i__]);
            /* L10: */
        }
        /* Reorder eigenvalues and transform Schur vectors */
        /* (CWorkspace: none) */
        /* (RWorkspace: none) */
        i__1 = *lwork - iwrk + 1;
        ctrsen_("N", jobvs, &bwork[1], n, &a[a_offset], lda, &vs[vs_offset], ldvs, &w[1], sdim, &s, &sep, &work[iwrk], &i__1, &icond);
    }
    if (wantvs)
    {
        /* Undo balancing */
        /* (CWorkspace: none) */
        /* (RWorkspace: need N) */
        cgebak_("P", "R", n, &ilo, &ihi, &rwork[ibal], n, &vs[vs_offset], ldvs, &ierr);
    }
    if (scalea)
    {
        /* Undo scaling for the Schur form of A */
        clascl_("U", &c__0, &c__0, &cscale, &anrm, n, n, &a[a_offset], lda, & ierr);
        i__1 = *lda + 1;
        ccopy_(n, &a[a_offset], &i__1, &w[1], &c__1);
    }
    work[1].r = (real) maxwrk;
    work[1].i = 0.f; // , expr subst
    return 0;
    /* End of CGEES */
}
/* cgees_ */
