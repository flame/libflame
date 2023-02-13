/* ../netlib/dggev.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
static doublereal c_b36 = 0.;
static doublereal c_b37 = 1.;
/* > \brief <b> DGGEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matr ices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGGEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggev.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggev.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggev.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, */
/* BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVL, JOBVR */
/* INTEGER INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/* $ B( LDB, * ), BETA( * ), VL( LDVL, * ), */
/* $ VR( LDVR, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B) */
/* > the generalized eigenvalues, and optionally, the left and/or right */
/* > generalized eigenvectors. */
/* > */
/* > A generalized eigenvalue for a pair of matrices (A,B) is a scalar */
/* > lambda or a ratio alpha/beta = lambda, such that A - lambda*B is */
/* > singular. It is usually represented as the pair (alpha,beta), as */
/* > there is a reasonable interpretation for beta=0, and even for both */
/* > being zero. */
/* > */
/* > The right eigenvector v(j) corresponding to the eigenvalue lambda(j) */
/* > of (A,B) satisfies */
/* > */
/* > A * v(j) = lambda(j) * B * v(j). */
/* > */
/* > The left eigenvector u(j) corresponding to the eigenvalue lambda(j) */
/* > of (A,B) satisfies */
/* > */
/* > u(j)**H * A = lambda(j) * u(j)**H * B . */
/* > */
/* > where u(j)**H is the conjugate-transpose of u(j). */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVL */
/* > \verbatim */
/* > JOBVL is CHARACTER*1 */
/* > = 'N': do not compute the left generalized eigenvectors;
*/
/* > = 'V': compute the left generalized eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVR */
/* > \verbatim */
/* > JOBVR is CHARACTER*1 */
/* > = 'N': do not compute the right generalized eigenvectors;
*/
/* > = 'V': compute the right generalized eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A, B, VL, and VR. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA, N) */
/* > On entry, the matrix A in the pair (A,B). */
/* > On exit, A has been overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB, N) */
/* > On entry, the matrix B in the pair (A,B). */
/* > On exit, B has been overwritten. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAR */
/* > \verbatim */
/* > ALPHAR is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* > ALPHAI is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is DOUBLE PRECISION array, dimension (N) */
/* > On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will */
/* > be the generalized eigenvalues. If ALPHAI(j) is zero, then */
/* > the j-th eigenvalue is real;
if positive, then the j-th and */
/* > (j+1)-st eigenvalues are a complex conjugate pair, with */
/* > ALPHAI(j+1) negative. */
/* > */
/* > Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j) */
/* > may easily over- or underflow, and BETA(j) may even be zero. */
/* > Thus, the user should avoid naively computing the ratio */
/* > alpha/beta. However, ALPHAR and ALPHAI will be always less */
/* > than and usually comparable with norm(A) in magnitude, and */
/* > BETA always less than and usually comparable with norm(B). */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* > VL is DOUBLE PRECISION array, dimension (LDVL,N) */
/* > If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/* > after another in the columns of VL, in the same order as */
/* > their eigenvalues. If the j-th eigenvalue is real, then */
/* > u(j) = VL(:,j), the j-th column of VL. If the j-th and */
/* > (j+1)-th eigenvalues form a complex conjugate pair, then */
/* > u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1). */
/* > Each eigenvector is scaled so the largest component has */
/* > f2c_abs(real part)+f2c_abs(imag. part)=1. */
/* > Not referenced if JOBVL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* > LDVL is INTEGER */
/* > The leading dimension of the matrix VL. LDVL >= 1, and */
/* > if JOBVL = 'V', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VR */
/* > \verbatim */
/* > VR is DOUBLE PRECISION array, dimension (LDVR,N) */
/* > If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/* > after another in the columns of VR, in the same order as */
/* > their eigenvalues. If the j-th eigenvalue is real, then */
/* > v(j) = VR(:,j), the j-th column of VR. If the j-th and */
/* > (j+1)-th eigenvalues form a complex conjugate pair, then */
/* > v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1). */
/* > Each eigenvector is scaled so the largest component has */
/* > f2c_abs(real part)+f2c_abs(imag. part)=1. */
/* > Not referenced if JOBVR = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* > LDVR is INTEGER */
/* > The leading dimension of the matrix VR. LDVR >= 1, and */
/* > if JOBVR = 'V', LDVR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= max(1,8*N). */
/* > For good performance, LWORK must generally be larger. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > = 1,...,N: */
/* > The QZ iteration failed. No eigenvectors have been */
/* > calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) */
/* > should be correct for j=INFO+1,...,N. */
/* > > N: =N+1: other than QZ iteration failed in DHGEQZ. */
/* > =N+2: error return from DTGEVC. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date April 2012 */
/* > \ingroup doubleGEeigen */
/* ===================================================================== */
/* Subroutine */
int dggev_(char *jobvl, char *jobvr, integer *n, doublereal * a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer jc, in, jr, ihi, ilo;
    doublereal eps;
    logical ilv;
    doublereal anrm, bnrm;
    integer ierr, itau;
    doublereal temp;
    logical ilvl, ilvr;
    integer iwrk;
    extern logical lsame_(char *, char *);
    integer ileft, icols, irows;
    extern /* Subroutine */
    int dlabad_(doublereal *, doublereal *), dggbak_( char *, char *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *), dggbal_(char *, integer *, doublereal *, integer *, doublereal *, integer *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *);
    extern doublereal dlamch_(char *), dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *);
    extern /* Subroutine */
    int dgghrd_(char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *), dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *);
    logical ilascl, ilbscl;
    extern /* Subroutine */
    int dgeqrf_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, integer *), dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *), dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *), dtgevc_(char *, char *, logical *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *, integer *, doublereal *, integer *);
    logical ldumma[1];
    char chtemp[1];
    doublereal bignum;
    extern /* Subroutine */
    int dhgeqz_(char *, char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer ijobvl, iright, ijobvr;
    extern /* Subroutine */
    int dorgqr_(integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, integer *);
    doublereal anrmto, bnrmto;
    extern /* Subroutine */
    int dormqr_(char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *);
    integer minwrk, maxwrk;
    doublereal smlnum;
    logical lquery;
    /* -- LAPACK driver routine (version 3.4.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* April 2012 */
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
    /* Decode the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alphar;
    --alphai;
    --beta;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    /* Function Body */
    if (lsame_(jobvl, "N"))
    {
        ijobvl = 1;
        ilvl = FALSE_;
    }
    else if (lsame_(jobvl, "V"))
    {
        ijobvl = 2;
        ilvl = TRUE_;
    }
    else
    {
        ijobvl = -1;
        ilvl = FALSE_;
    }
    if (lsame_(jobvr, "N"))
    {
        ijobvr = 1;
        ilvr = FALSE_;
    }
    else if (lsame_(jobvr, "V"))
    {
        ijobvr = 2;
        ilvr = TRUE_;
    }
    else
    {
        ijobvr = -1;
        ilvr = FALSE_;
    }
    ilv = ilvl || ilvr;
    /* Test the input arguments */
    *info = 0;
    lquery = *lwork == -1;
    if (ijobvl <= 0)
    {
        *info = -1;
    }
    else if (ijobvr <= 0)
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
    else if (*ldb < max(1,*n))
    {
        *info = -7;
    }
    else if (*ldvl < 1 || ilvl && *ldvl < *n)
    {
        *info = -12;
    }
    else if (*ldvr < 1 || ilvr && *ldvr < *n)
    {
        *info = -14;
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* NB refers to the optimal block size for the immediately */
    /* following subroutine, as returned by ILAENV. The workspace is */
    /* computed assuming ILO = 1 and IHI = N, the worst case.) */
    if (*info == 0)
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n << 3; // , expr subst
        minwrk = max(i__1,i__2);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n * (ilaenv_(&c__1, "DGEQRF", " ", n, &c__1, n, & c__0) + 7); // , expr subst
        maxwrk = max(i__1,i__2);
        /* Computing MAX */
        i__1 = maxwrk;
        i__2 = *n * (ilaenv_(&c__1, "DORMQR", " ", n, &c__1, n, &c__0) + 7); // , expr subst
        maxwrk = max(i__1,i__2);
        if (ilvl)
        {
            /* Computing MAX */
            i__1 = maxwrk;
            i__2 = *n * (ilaenv_(&c__1, "DORGQR", " ", n, & c__1, n, &c_n1) + 7); // , expr subst
            maxwrk = max(i__1,i__2);
        }
        work[1] = (doublereal) maxwrk;
        if (*lwork < minwrk && ! lquery)
        {
            *info = -16;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGGEV ", &i__1);
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
    eps = dlamch_("P");
    smlnum = dlamch_("S");
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = 1. / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1]);
    ilascl = FALSE_;
    if (anrm > 0. && anrm < smlnum)
    {
        anrmto = smlnum;
        ilascl = TRUE_;
    }
    else if (anrm > bignum)
    {
        anrmto = bignum;
        ilascl = TRUE_;
    }
    if (ilascl)
    {
        dlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, & ierr);
    }
    /* Scale B if max element outside range [SMLNUM,BIGNUM] */
    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1]);
    ilbscl = FALSE_;
    if (bnrm > 0. && bnrm < smlnum)
    {
        bnrmto = smlnum;
        ilbscl = TRUE_;
    }
    else if (bnrm > bignum)
    {
        bnrmto = bignum;
        ilbscl = TRUE_;
    }
    if (ilbscl)
    {
        dlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, & ierr);
    }
    /* Permute the matrices A, B to isolate eigenvalues if possible */
    /* (Workspace: need 6*N) */
    ileft = 1;
    iright = *n + 1;
    iwrk = iright + *n;
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[ ileft], &work[iright], &work[iwrk], &ierr);
    /* Reduce B to triangular form (QR decomposition of B) */
    /* (Workspace: need N, prefer N*NB) */
    irows = ihi + 1 - ilo;
    if (ilv)
    {
        icols = *n + 1 - ilo;
    }
    else
    {
        icols = irows;
    }
    itau = iwrk;
    iwrk = itau + irows;
    i__1 = *lwork + 1 - iwrk;
    dgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[ iwrk], &i__1, &ierr);
    /* Apply the orthogonal transformation to matrix A */
    /* (Workspace: need N, prefer N*NB) */
    i__1 = *lwork + 1 - iwrk;
    dormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, & work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, & ierr);
    /* Initialize VL */
    /* (Workspace: need N, prefer N*NB) */
    if (ilvl)
    {
        dlaset_("Full", n, n, &c_b36, &c_b37, &vl[vl_offset], ldvl) ;
        if (irows > 1)
        {
            i__1 = irows - 1;
            i__2 = irows - 1;
            dlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vl[ ilo + 1 + ilo * vl_dim1], ldvl);
        }
        i__1 = *lwork + 1 - iwrk;
        dorgqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[ itau], &work[iwrk], &i__1, &ierr);
    }
    /* Initialize VR */
    if (ilvr)
    {
        dlaset_("Full", n, n, &c_b36, &c_b37, &vr[vr_offset], ldvr) ;
    }
    /* Reduce to generalized Hessenberg form */
    /* (Workspace: none needed) */
    if (ilv)
    {
        /* Eigenvectors requested -- work on whole matrix. */
        dgghrd_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &ierr);
    }
    else
    {
        dgghrd_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda, &b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[ vr_offset], ldvr, &ierr);
    }
    /* Perform QZ algorithm (Compute eigenvalues, and optionally, the */
    /* Schur forms and Schur vectors) */
    /* (Workspace: need N) */
    iwrk = itau;
    if (ilv)
    {
        *(unsigned char *)chtemp = 'S';
    }
    else
    {
        *(unsigned char *)chtemp = 'E';
    }
    i__1 = *lwork + 1 - iwrk;
    dhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[ b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[iwrk], &i__1, &ierr);
    if (ierr != 0)
    {
        if (ierr > 0 && ierr <= *n)
        {
            *info = ierr;
        }
        else if (ierr > *n && ierr <= *n << 1)
        {
            *info = ierr - *n;
        }
        else
        {
            *info = *n + 1;
        }
        goto L110;
    }
    /* Compute Eigenvectors */
    /* (Workspace: need 6*N) */
    if (ilv)
    {
        if (ilvl)
        {
            if (ilvr)
            {
                *(unsigned char *)chtemp = 'B';
            }
            else
            {
                *(unsigned char *)chtemp = 'L';
            }
        }
        else
        {
            *(unsigned char *)chtemp = 'R';
        }
        dtgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[ iwrk], &ierr);
        if (ierr != 0)
        {
            *info = *n + 2;
            goto L110;
        }
        /* Undo balancing on VL and VR and normalization */
        /* (Workspace: none needed) */
        if (ilvl)
        {
            dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, & vl[vl_offset], ldvl, &ierr);
            i__1 = *n;
            for (jc = 1;
                    jc <= i__1;
                    ++jc)
            {
                if (alphai[jc] < 0.)
                {
                    goto L50;
                }
                temp = 0.;
                if (alphai[jc] == 0.)
                {
                    i__2 = *n;
                    for (jr = 1;
                            jr <= i__2;
                            ++jr)
                    {
                        /* Computing MAX */
                        d__2 = temp;
                        d__3 = (d__1 = vl[jr + jc * vl_dim1], f2c_abs(d__1)); // , expr subst
                        temp = max(d__2,d__3);
                        /* L10: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for (jr = 1;
                            jr <= i__2;
                            ++jr)
                    {
                        /* Computing MAX */
                        d__3 = temp;
                        d__4 = (d__1 = vl[jr + jc * vl_dim1], f2c_abs(d__1)) + (d__2 = vl[jr + (jc + 1) * vl_dim1], f2c_abs(d__2)); // , expr subst
                        temp = max(d__3,d__4);
                        /* L20: */
                    }
                }
                if (temp < smlnum)
                {
                    goto L50;
                }
                temp = 1. / temp;
                if (alphai[jc] == 0.)
                {
                    i__2 = *n;
                    for (jr = 1;
                            jr <= i__2;
                            ++jr)
                    {
                        vl[jr + jc * vl_dim1] *= temp;
                        /* L30: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for (jr = 1;
                            jr <= i__2;
                            ++jr)
                    {
                        vl[jr + jc * vl_dim1] *= temp;
                        vl[jr + (jc + 1) * vl_dim1] *= temp;
                        /* L40: */
                    }
                }
L50:
                ;
            }
        }
        if (ilvr)
        {
            dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, & vr[vr_offset], ldvr, &ierr);
            i__1 = *n;
            for (jc = 1;
                    jc <= i__1;
                    ++jc)
            {
                if (alphai[jc] < 0.)
                {
                    goto L100;
                }
                temp = 0.;
                if (alphai[jc] == 0.)
                {
                    i__2 = *n;
                    for (jr = 1;
                            jr <= i__2;
                            ++jr)
                    {
                        /* Computing MAX */
                        d__2 = temp;
                        d__3 = (d__1 = vr[jr + jc * vr_dim1], f2c_abs(d__1)); // , expr subst
                        temp = max(d__2,d__3);
                        /* L60: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for (jr = 1;
                            jr <= i__2;
                            ++jr)
                    {
                        /* Computing MAX */
                        d__3 = temp;
                        d__4 = (d__1 = vr[jr + jc * vr_dim1], f2c_abs(d__1)) + (d__2 = vr[jr + (jc + 1) * vr_dim1], f2c_abs(d__2)); // , expr subst
                        temp = max(d__3,d__4);
                        /* L70: */
                    }
                }
                if (temp < smlnum)
                {
                    goto L100;
                }
                temp = 1. / temp;
                if (alphai[jc] == 0.)
                {
                    i__2 = *n;
                    for (jr = 1;
                            jr <= i__2;
                            ++jr)
                    {
                        vr[jr + jc * vr_dim1] *= temp;
                        /* L80: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for (jr = 1;
                            jr <= i__2;
                            ++jr)
                    {
                        vr[jr + jc * vr_dim1] *= temp;
                        vr[jr + (jc + 1) * vr_dim1] *= temp;
                        /* L90: */
                    }
                }
L100:
                ;
            }
        }
        /* End of eigenvector calculation */
    }
    /* Undo scaling if necessary */
L110:
    if (ilascl)
    {
        dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, & ierr);
        dlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, & ierr);
    }
    if (ilbscl)
    {
        dlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, & ierr);
    }
    work[1] = (doublereal) maxwrk;
    return 0;
    /* End of DGGEV */
}
/* dggev_ */
