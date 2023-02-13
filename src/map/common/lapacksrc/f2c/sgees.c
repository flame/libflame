/* ../netlib/sgees.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
/* > \brief <b> SGEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors f or GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGEES + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgees.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgees.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgees.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI, */
/* VS, LDVS, WORK, LWORK, BWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVS, SORT */
/* INTEGER INFO, LDA, LDVS, LWORK, N, SDIM */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL BWORK( * ) */
/* REAL A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ), */
/* $ WR( * ) */
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
/* > SGEES computes for an N-by-N real nonsymmetric matrix A, the */
/* > eigenvalues, the real Schur form T, and, optionally, the matrix of */
/* > Schur vectors Z. This gives the Schur factorization A = Z*T*(Z**T). */
/* > */
/* > Optionally, it also orders the eigenvalues on the diagonal of the */
/* > real Schur form so that selected eigenvalues are at the top left. */
/* > The leading columns of Z then form an orthonormal basis for the */
/* > invariant subspace corresponding to the selected eigenvalues. */
/* > */
/* > A matrix is in real Schur form if it is upper quasi-triangular with */
/* > 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the */
/* > form */
/* > [ a b ] */
/* > [ c a ] */
/* > */
/* > where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc). */
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
/* > = 'N': Eigenvalues are not ordered;
*/
/* > = 'S': Eigenvalues are ordered (see SELECT). */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* > SELECT is LOGICAL FUNCTION of two REAL arguments */
/* > SELECT must be declared EXTERNAL in the calling subroutine. */
/* > If SORT = 'S', SELECT is used to select eigenvalues to sort */
/* > to the top left of the Schur form. */
/* > If SORT = 'N', SELECT is not referenced. */
/* > An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if */
/* > SELECT(WR(j),WI(j)) is true;
i.e., if either one of a complex */
/* > conjugate pair of eigenvalues is selected, then both complex */
/* > eigenvalues are selected. */
/* > Note that a selected complex eigenvalue may no longer */
/* > satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since */
/* > ordering may change the value of complex eigenvalues */
/* > (especially if the eigenvalue is ill-conditioned);
in this */
/* > case INFO is set to N+2 (see INFO below). */
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
/* > On entry, the N-by-N matrix A. */
/* > On exit, A has been overwritten by its real Schur form T. */
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
/* > If SORT = 'S', SDIM = number of eigenvalues (after sorting) */
/* > for which SELECT is true. (Complex conjugate */
/* > pairs for which SELECT is true for either */
/* > eigenvalue count as 2.) */
/* > \endverbatim */
/* > */
/* > \param[out] WR */
/* > \verbatim */
/* > WR is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* > WI is REAL array, dimension (N) */
/* > WR and WI contain the real and imaginary parts, */
/* > respectively, of the computed eigenvalues in the same order */
/* > that they appear on the diagonal of the output Schur form T. */
/* > Complex conjugate pairs of eigenvalues will appear */
/* > consecutively with the eigenvalue having the positive */
/* > imaginary part first. */
/* > \endverbatim */
/* > */
/* > \param[out] VS */
/* > \verbatim */
/* > VS is REAL array, dimension (LDVS,N) */
/* > If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur */
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
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) contains the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= max(1,3*N). */
/* > For good performance, LWORK must generally be larger. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
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
elements 1:ILO-1 and i+1:N of WR and WI */
/* > contain those eigenvalues which have converged;
if */
/* > JOBVS = 'V', VS contains the matrix which reduces A */
/* > to its partially converged Schur form. */
/* > = N+1: the eigenvalues could not be reordered because some */
/* > eigenvalues were too close to separate (the problem */
/* > is very ill-conditioned);
*/
/* > = N+2: after reordering, roundoff changed values of some */
/* > complex eigenvalues so that leading eigenvalues in */
/* > the Schur form no longer satisfy SELECT=.TRUE. This */
/* > could also be caused by underflow due to scaling. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realGEeigen */
/* ===================================================================== */
/* Subroutine */
int sgees_(char *jobvs, char *sort, L_fp select, integer *n, real *a, integer *lda, integer *sdim, real *wr, real *wi, real *vs, integer *ldvs, real *work, integer *lwork, logical *bwork, integer * info)
{
    /* System generated locals */
    integer a_dim1, a_offset, vs_dim1, vs_offset, i__1, i__2, i__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__;
    real s;
    integer i1, i2, ip, ihi, ilo;
    real dum[1], eps, sep;
    integer ibal;
    real anrm;
    integer idum[1], ierr, itau, iwrk, inxt, icond, ieval;
    extern logical lsame_(char *, char *);
    logical cursl;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *), sswap_(integer *, real *, integer *, real *, integer * );
    logical lst2sl;
    extern /* Subroutine */
    int slabad_(real *, real *);
    logical scalea;
    real cscale;
    extern /* Subroutine */
    int sgebak_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, integer *, integer *), sgebal_(char *, integer *, real *, integer *, integer *, integer *, real *, integer *);
    extern real slamch_(char *), slange_(char *, integer *, integer *, real *, integer *, real *);
    extern /* Subroutine */
    int sgehrd_(integer *, integer *, integer *, real *, integer *, real *, real *, integer *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    real bignum;
    extern /* Subroutine */
    int slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *, integer *, integer *), slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *);
    logical lastsl;
    extern /* Subroutine */
    int sorghr_(integer *, integer *, integer *, real *, integer *, real *, real *, integer *, integer *), shseqr_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, real *, integer *, real *, integer *, integer *);
    integer minwrk, maxwrk;
    real smlnum;
    integer hswork;
    extern /* Subroutine */
    int strsen_(char *, char *, logical *, integer *, real *, integer *, real *, integer *, real *, real *, integer *, real *, real *, real *, integer *, integer *, integer *, integer * );
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
    --wr;
    --wi;
    vs_dim1 = *ldvs;
    vs_offset = 1 + vs_dim1;
    vs -= vs_offset;
    --work;
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
        *info = -11;
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* NB refers to the optimal block size for the immediately */
    /* following subroutine, as returned by ILAENV. */
    /* HSWORK refers to the workspace preferred by SHSEQR, as */
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
            maxwrk = (*n << 1) + *n * ilaenv_(&c__1, "SGEHRD", " ", n, &c__1, n, &c__0);
            minwrk = *n * 3;
            shseqr_("S", jobvs, n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[1] , &vs[vs_offset], ldvs, &work[1], &c_n1, &ieval);
            hswork = work[1];
            if (! wantvs)
            {
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n + hswork; // , expr subst
                maxwrk = max(i__1,i__2);
            }
            else
            {
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, "SORGHR", " ", n, &c__1, n, &c_n1); // , expr subst
                maxwrk = max(i__1,i__2);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n + hswork; // , expr subst
                maxwrk = max(i__1,i__2);
            }
        }
        work[1] = (real) maxwrk;
        if (*lwork < minwrk && ! lquery)
        {
            *info = -13;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGEES ", &i__1);
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
    anrm = slange_("M", n, n, &a[a_offset], lda, dum);
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
        slascl_("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, & ierr);
    }
    /* Permute the matrix to make it more nearly triangular */
    /* (Workspace: need N) */
    ibal = 1;
    sgebal_("P", n, &a[a_offset], lda, &ilo, &ihi, &work[ibal], &ierr);
    /* Reduce to upper Hessenberg form */
    /* (Workspace: need 3*N, prefer 2*N+N*NB) */
    itau = *n + ibal;
    iwrk = *n + itau;
    i__1 = *lwork - iwrk + 1;
    sgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &ierr);
    if (wantvs)
    {
        /* Copy Householder vectors to VS */
        slacpy_("L", n, n, &a[a_offset], lda, &vs[vs_offset], ldvs) ;
        /* Generate orthogonal matrix in VS */
        /* (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */
        i__1 = *lwork - iwrk + 1;
        sorghr_(n, &ilo, &ihi, &vs[vs_offset], ldvs, &work[itau], &work[iwrk], &i__1, &ierr);
    }
    *sdim = 0;
    /* Perform QR iteration, accumulating Schur vectors in VS if desired */
    /* (Workspace: need N+1, prefer N+HSWORK (see comments) ) */
    iwrk = itau;
    i__1 = *lwork - iwrk + 1;
    shseqr_("S", jobvs, n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &vs[ vs_offset], ldvs, &work[iwrk], &i__1, &ieval);
    if (ieval > 0)
    {
        *info = ieval;
    }
    /* Sort eigenvalues if desired */
    if (wantst && *info == 0)
    {
        if (scalea)
        {
            slascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &wr[1], n, & ierr);
            slascl_("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &wi[1], n, & ierr);
        }
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            bwork[i__] = (*select)(&wr[i__], &wi[i__]);
            /* L10: */
        }
        /* Reorder eigenvalues and transform Schur vectors */
        /* (Workspace: none needed) */
        i__1 = *lwork - iwrk + 1;
        strsen_("N", jobvs, &bwork[1], n, &a[a_offset], lda, &vs[vs_offset], ldvs, &wr[1], &wi[1], sdim, &s, &sep, &work[iwrk], &i__1, idum, &c__1, &icond);
        if (icond > 0)
        {
            *info = *n + icond;
        }
    }
    if (wantvs)
    {
        /* Undo balancing */
        /* (Workspace: need N) */
        sgebak_("P", "R", n, &ilo, &ihi, &work[ibal], n, &vs[vs_offset], ldvs, &ierr);
    }
    if (scalea)
    {
        /* Undo scaling for the Schur form of A */
        slascl_("H", &c__0, &c__0, &cscale, &anrm, n, n, &a[a_offset], lda, & ierr);
        i__1 = *lda + 1;
        scopy_(n, &a[a_offset], &i__1, &wr[1], &c__1);
        if (cscale == smlnum)
        {
            /* If scaling back towards underflow, adjust WI if an */
            /* offdiagonal element of a 2-by-2 block in the Schur form */
            /* underflows. */
            if (ieval > 0)
            {
                i1 = ieval + 1;
                i2 = ihi - 1;
                i__1 = ilo - 1;
                /* Computing MAX */
                i__3 = ilo - 1;
                i__2 = max(i__3,1);
                slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[ 1], &i__2, &ierr);
            }
            else if (wantst)
            {
                i1 = 1;
                i2 = *n - 1;
            }
            else
            {
                i1 = ilo;
                i2 = ihi - 1;
            }
            inxt = i1 - 1;
            i__1 = i2;
            for (i__ = i1;
                    i__ <= i__1;
                    ++i__)
            {
                if (i__ < inxt)
                {
                    goto L20;
                }
                if (wi[i__] == 0.f)
                {
                    inxt = i__ + 1;
                }
                else
                {
                    if (a[i__ + 1 + i__ * a_dim1] == 0.f)
                    {
                        wi[i__] = 0.f;
                        wi[i__ + 1] = 0.f;
                    }
                    else if (a[i__ + 1 + i__ * a_dim1] != 0.f && a[i__ + ( i__ + 1) * a_dim1] == 0.f)
                    {
                        wi[i__] = 0.f;
                        wi[i__ + 1] = 0.f;
                        if (i__ > 1)
                        {
                            i__2 = i__ - 1;
                            sswap_(&i__2, &a[i__ * a_dim1 + 1], &c__1, &a[( i__ + 1) * a_dim1 + 1], &c__1);
                        }
                        if (*n > i__ + 1)
                        {
                            i__2 = *n - i__ - 1;
                            sswap_(&i__2, &a[i__ + (i__ + 2) * a_dim1], lda, & a[i__ + 1 + (i__ + 2) * a_dim1], lda);
                        }
                        if (wantvs)
                        {
                            sswap_(n, &vs[i__ * vs_dim1 + 1], &c__1, &vs[(i__ + 1) * vs_dim1 + 1], &c__1);
                        }
                        a[i__ + (i__ + 1) * a_dim1] = a[i__ + 1 + i__ * a_dim1];
                        a[i__ + 1 + i__ * a_dim1] = 0.f;
                    }
                    inxt = i__ + 2;
                }
L20:
                ;
            }
        }
        /* Undo scaling for the imaginary part of the eigenvalues */
        i__1 = *n - ieval;
        /* Computing MAX */
        i__3 = *n - ieval;
        i__2 = max(i__3,1);
        slascl_("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[ieval + 1], &i__2, &ierr);
    }
    if (wantst && *info == 0)
    {
        /* Check if reordering successful */
        lastsl = TRUE_;
        lst2sl = TRUE_;
        *sdim = 0;
        ip = 0;
        i__1 = *n;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            cursl = (*select)(&wr[i__], &wi[i__]);
            if (wi[i__] == 0.f)
            {
                if (cursl)
                {
                    ++(*sdim);
                }
                ip = 0;
                if (cursl && ! lastsl)
                {
                    *info = *n + 2;
                }
            }
            else
            {
                if (ip == 1)
                {
                    /* Last eigenvalue of conjugate pair */
                    cursl = cursl || lastsl;
                    lastsl = cursl;
                    if (cursl)
                    {
                        *sdim += 2;
                    }
                    ip = -1;
                    if (cursl && ! lst2sl)
                    {
                        *info = *n + 2;
                    }
                }
                else
                {
                    /* First eigenvalue of conjugate pair */
                    ip = 1;
                }
            }
            lst2sl = lastsl;
            lastsl = cursl;
            /* L30: */
        }
    }
    work[1] = (real) maxwrk;
    return 0;
    /* End of SGEES */
}
/* sgees_ */
