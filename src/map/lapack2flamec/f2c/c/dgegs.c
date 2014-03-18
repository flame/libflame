/* ../netlib/dgegs.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b36 = 0.;
static doublereal c_b37 = 1.;
/* > \brief <b> DGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat rices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGEGS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgegs.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgegs.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgegs.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHAR, */
/* ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, */
/* LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVSL, JOBVSR */
/* INTEGER INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/* $ B( LDB, * ), BETA( * ), VSL( LDVSL, * ), */
/* $ VSR( LDVSR, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine DGGES. */
/* > */
/* > DGEGS computes the eigenvalues, real Schur form, and, optionally, */
/* > left and or/right Schur vectors of a real matrix pair (A,B). */
/* > Given two square matrices A and B, the generalized real Schur */
/* > factorization has the form */
/* > */
/* > A = Q*S*Z**T, B = Q*T*Z**T */
/* > */
/* > where Q and Z are orthogonal matrices, T is upper triangular, and S */
/* > is an upper quasi-triangular matrix with 1-by-1 and 2-by-2 diagonal */
/* > blocks, the 2-by-2 blocks corresponding to complex conjugate pairs */
/* > of eigenvalues of (A,B). The columns of Q are the left Schur vectors */
/* > and the columns of Z are the right Schur vectors. */
/* > */
/* > If only the eigenvalues of (A,B) are needed, the driver routine */
/* > DGEGV should be used instead. See DGEGV for a description of the */
/* > eigenvalues of the generalized nonsymmetric eigenvalue problem */
/* > (GNEP). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVSL */
/* > \verbatim */
/* > JOBVSL is CHARACTER*1 */
/* > = 'N': do not compute the left Schur vectors;
*/
/* > = 'V': compute the left Schur vectors (returned in VSL). */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVSR */
/* > \verbatim */
/* > JOBVSR is CHARACTER*1 */
/* > = 'N': do not compute the right Schur vectors;
*/
/* > = 'V': compute the right Schur vectors (returned in VSR). */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A, B, VSL, and VSR. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA, N) */
/* > On entry, the matrix A. */
/* > On exit, the upper quasi-triangular matrix S from the */
/* > generalized real Schur factorization. */
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
/* > On entry, the matrix B. */
/* > On exit, the upper triangular matrix T from the generalized */
/* > real Schur factorization. */
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
/* > The real parts of each scalar alpha defining an eigenvalue */
/* > of GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* > ALPHAI is DOUBLE PRECISION array, dimension (N) */
/* > The imaginary parts of each scalar alpha defining an */
/* > eigenvalue of GNEP. If ALPHAI(j) is zero, then the j-th */
/* > eigenvalue is real;
if positive, then the j-th and (j+1)-st */
/* > eigenvalues are a complex conjugate pair, with */
/* > ALPHAI(j+1) = -ALPHAI(j). */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is DOUBLE PRECISION array, dimension (N) */
/* > The scalars beta that define the eigenvalues of GNEP. */
/* > Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and */
/* > beta = BETA(j) represent the j-th eigenvalue of the matrix */
/* > pair (A,B), in one of the forms lambda = alpha/beta or */
/* > mu = beta/alpha. Since either lambda or mu may overflow, */
/* > they should not, in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[out] VSL */
/* > \verbatim */
/* > VSL is DOUBLE PRECISION array, dimension (LDVSL,N) */
/* > If JOBVSL = 'V', the matrix of left Schur vectors Q. */
/* > Not referenced if JOBVSL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSL */
/* > \verbatim */
/* > LDVSL is INTEGER */
/* > The leading dimension of the matrix VSL. LDVSL >=1, and */
/* > if JOBVSL = 'V', LDVSL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VSR */
/* > \verbatim */
/* > VSR is DOUBLE PRECISION array, dimension (LDVSR,N) */
/* > If JOBVSR = 'V', the matrix of right Schur vectors Z. */
/* > Not referenced if JOBVSR = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSR */
/* > \verbatim */
/* > LDVSR is INTEGER */
/* > The leading dimension of the matrix VSR. LDVSR >= 1, and */
/* > if JOBVSR = 'V', LDVSR >= N. */
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
/* > The dimension of the array WORK. LWORK >= max(1,4*N). */
/* > For good performance, LWORK must generally be larger. */
/* > To compute the optimal value of LWORK, call ILAENV to get */
/* > blocksizes (for DGEQRF, DORMQR, and DORGQR.) Then compute: */
/* > NB -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR */
/* > The optimal LWORK is 2*N + N*(NB+1). */
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
/* > The QZ iteration failed. (A,B) are not in Schur */
/* > form, but ALPHAR(j), ALPHAI(j), and BETA(j) should */
/* > be correct for j=INFO+1,...,N. */
/* > > N: errors that usually indicate LAPACK problems: */
/* > =N+1: error return from DGGBAL */
/* > =N+2: error return from DGEQRF */
/* > =N+3: error return from DORMQR */
/* > =N+4: error return from DORGQR */
/* > =N+5: error return from DGGHRD */
/* > =N+6: error return from DHGEQZ (other than failed */
/* > iteration) */
/* > =N+7: error return from DGGBAK (computing VSL) */
/* > =N+8: error return from DGGBAK (computing VSR) */
/* > =N+9: error return from DLASCL (various places) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup doubleGEeigen */
/* ===================================================================== */
/* Subroutine */
int dgegs_(char *jobvsl, char *jobvsr, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal * alphar, doublereal *alphai, doublereal *beta, doublereal *vsl, integer *ldvsl, doublereal *vsr, integer *ldvsr, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, vsr_dim1, vsr_offset, i__1, i__2;
    /* Local variables */
    integer nb, nb1, nb2, nb3, ihi, ilo;
    doublereal eps, anrm, bnrm;
    integer itau, lopt;
    extern logical lsame_(char *, char *);
    integer ileft, iinfo, icols;
    logical ilvsl;
    integer iwork;
    logical ilvsr;
    integer irows;
    extern /* Subroutine */
    int dggbak_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *), dggbal_(char *, integer *, doublereal *, integer *, doublereal *, integer *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *);
    extern doublereal dlamch_(char *), dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *);
    extern /* Subroutine */
    int dgghrd_(char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *), dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, integer *);
    logical ilascl, ilbscl;
    extern /* Subroutine */
    int dgeqrf_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, integer *), dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal safmin;
    extern /* Subroutine */
    int dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    doublereal bignum;
    extern /* Subroutine */
    int dhgeqz_(char *, char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, integer *, integer *);
    integer ijobvl, iright, ijobvr;
    extern /* Subroutine */
    int dorgqr_(integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, integer *);
    doublereal anrmto;
    integer lwkmin;
    doublereal bnrmto;
    extern /* Subroutine */
    int dormqr_(char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, integer *);
    doublereal smlnum;
    integer lwkopt;
    logical lquery;
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
    vsl_dim1 = *ldvsl;
    vsl_offset = 1 + vsl_dim1;
    vsl -= vsl_offset;
    vsr_dim1 = *ldvsr;
    vsr_offset = 1 + vsr_dim1;
    vsr -= vsr_offset;
    --work;
    /* Function Body */
    if (lsame_(jobvsl, "N"))
    {
        ijobvl = 1;
        ilvsl = FALSE_;
    }
    else if (lsame_(jobvsl, "V"))
    {
        ijobvl = 2;
        ilvsl = TRUE_;
    }
    else
    {
        ijobvl = -1;
        ilvsl = FALSE_;
    }
    if (lsame_(jobvsr, "N"))
    {
        ijobvr = 1;
        ilvsr = FALSE_;
    }
    else if (lsame_(jobvsr, "V"))
    {
        ijobvr = 2;
        ilvsr = TRUE_;
    }
    else
    {
        ijobvr = -1;
        ilvsr = FALSE_;
    }
    /* Test the input arguments */
    /* Computing MAX */
    i__1 = *n << 2;
    lwkmin = max(i__1,1);
    lwkopt = lwkmin;
    work[1] = (doublereal) lwkopt;
    lquery = *lwork == -1;
    *info = 0;
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
    else if (*ldvsl < 1 || ilvsl && *ldvsl < *n)
    {
        *info = -12;
    }
    else if (*ldvsr < 1 || ilvsr && *ldvsr < *n)
    {
        *info = -14;
    }
    else if (*lwork < lwkmin && ! lquery)
    {
        *info = -16;
    }
    if (*info == 0)
    {
        nb1 = ilaenv_(&c__1, "DGEQRF", " ", n, n, &c_n1, &c_n1);
        nb2 = ilaenv_(&c__1, "DORMQR", " ", n, n, n, &c_n1);
        nb3 = ilaenv_(&c__1, "DORGQR", " ", n, n, n, &c_n1);
        /* Computing MAX */
        i__1 = max(nb1,nb2);
        nb = max(i__1,nb3);
        lopt = (*n << 1) + *n * (nb + 1);
        work[1] = (doublereal) lopt;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGEGS ", &i__1);
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
    eps = dlamch_("E") * dlamch_("B");
    safmin = dlamch_("S");
    smlnum = *n * safmin / eps;
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
        dlascl_("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &a[a_offset], lda, & iinfo);
        if (iinfo != 0)
        {
            *info = *n + 9;
            return 0;
        }
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
        dlascl_("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, & iinfo);
        if (iinfo != 0)
        {
            *info = *n + 9;
            return 0;
        }
    }
    /* Permute the matrix to make it more nearly triangular */
    /* Workspace layout: (2*N words -- "work..." not actually used) */
    /* left_permutation, right_permutation, work... */
    ileft = 1;
    iright = *n + 1;
    iwork = iright + *n;
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[ ileft], &work[iright], &work[iwork], &iinfo);
    if (iinfo != 0)
    {
        *info = *n + 1;
        goto L10;
    }
    /* Reduce B to triangular form, and initialize VSL and/or VSR */
    /* Workspace layout: ("work..." must have at least N words) */
    /* left_permutation, right_permutation, tau, work... */
    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = iwork;
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;
    dgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[ iwork], &i__1, &iinfo);
    if (iinfo >= 0)
    {
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = (integer) work[iwork] + iwork - 1; // , expr subst
        lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0)
    {
        *info = *n + 2;
        goto L10;
    }
    i__1 = *lwork + 1 - iwork;
    dormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, & work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, & iinfo);
    if (iinfo >= 0)
    {
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = (integer) work[iwork] + iwork - 1; // , expr subst
        lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0)
    {
        *info = *n + 3;
        goto L10;
    }
    if (ilvsl)
    {
        dlaset_("Full", n, n, &c_b36, &c_b37, &vsl[vsl_offset], ldvsl);
        i__1 = irows - 1;
        i__2 = irows - 1;
        dlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vsl[ilo + 1 + ilo * vsl_dim1], ldvsl);
        i__1 = *lwork + 1 - iwork;
        dorgqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, & work[itau], &work[iwork], &i__1, &iinfo);
        if (iinfo >= 0)
        {
            /* Computing MAX */
            i__1 = lwkopt;
            i__2 = (integer) work[iwork] + iwork - 1; // , expr subst
            lwkopt = max(i__1,i__2);
        }
        if (iinfo != 0)
        {
            *info = *n + 4;
            goto L10;
        }
    }
    if (ilvsr)
    {
        dlaset_("Full", n, n, &c_b36, &c_b37, &vsr[vsr_offset], ldvsr);
    }
    /* Reduce to generalized Hessenberg form */
    dgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &iinfo);
    if (iinfo != 0)
    {
        *info = *n + 5;
        goto L10;
    }
    /* Perform QZ algorithm, computing Schur vectors if desired */
    /* Workspace layout: ("work..." must have at least 1 word) */
    /* left_permutation, right_permutation, work... */
    iwork = itau;
    i__1 = *lwork + 1 - iwork;
    dhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[ b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[vsl_offset] , ldvsl, &vsr[vsr_offset], ldvsr, &work[iwork], &i__1, &iinfo);
    if (iinfo >= 0)
    {
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = (integer) work[iwork] + iwork - 1; // , expr subst
        lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0)
    {
        if (iinfo > 0 && iinfo <= *n)
        {
            *info = iinfo;
        }
        else if (iinfo > *n && iinfo <= *n << 1)
        {
            *info = iinfo - *n;
        }
        else
        {
            *info = *n + 6;
        }
        goto L10;
    }
    /* Apply permutation to VSL and VSR */
    if (ilvsl)
    {
        dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsl[ vsl_offset], ldvsl, &iinfo);
        if (iinfo != 0)
        {
            *info = *n + 7;
            goto L10;
        }
    }
    if (ilvsr)
    {
        dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vsr[ vsr_offset], ldvsr, &iinfo);
        if (iinfo != 0)
        {
            *info = *n + 8;
            goto L10;
        }
    }
    /* Undo scaling */
    if (ilascl)
    {
        dlascl_("H", &c_n1, &c_n1, &anrmto, &anrm, n, n, &a[a_offset], lda, & iinfo);
        if (iinfo != 0)
        {
            *info = *n + 9;
            return 0;
        }
        dlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alphar[1], n, & iinfo);
        if (iinfo != 0)
        {
            *info = *n + 9;
            return 0;
        }
        dlascl_("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alphai[1], n, & iinfo);
        if (iinfo != 0)
        {
            *info = *n + 9;
            return 0;
        }
    }
    if (ilbscl)
    {
        dlascl_("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, & iinfo);
        if (iinfo != 0)
        {
            *info = *n + 9;
            return 0;
        }
        dlascl_("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &beta[1], n, & iinfo);
        if (iinfo != 0)
        {
            *info = *n + 9;
            return 0;
        }
    }
L10:
    work[1] = (doublereal) lwkopt;
    return 0;
    /* End of DGEGS */
}
/* dgegs_ */
