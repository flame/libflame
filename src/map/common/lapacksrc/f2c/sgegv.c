/* ../netlib/sgegv.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static real c_b27 = 1.f;
static real c_b38 = 0.f;
/* > \brief <b> SGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE mat rices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGEGV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgegv.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgegv.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgegv.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, */
/* BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVL, JOBVR */
/* INTEGER INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/* $ B( LDB, * ), BETA( * ), VL( LDVL, * ), */
/* $ VR( LDVR, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine SGGEV. */
/* > */
/* > SGEGV computes the eigenvalues and, optionally, the left and/or right */
/* > eigenvectors of a real matrix pair (A,B). */
/* > Given two square matrices A and B, */
/* > the generalized nonsymmetric eigenvalue problem (GNEP) is to find the */
/* > eigenvalues lambda and corresponding (non-zero) eigenvectors x such */
/* > that */
/* > */
/* > A*x = lambda*B*x. */
/* > */
/* > An alternate form is to find the eigenvalues mu and corresponding */
/* > eigenvectors y such that */
/* > */
/* > mu*A*y = B*y. */
/* > */
/* > These two forms are equivalent with mu = 1/lambda and x = y if */
/* > neither lambda nor mu is zero. In order to deal with the case that */
/* > lambda or mu is zero or small, two values alpha and beta are returned */
/* > for each eigenvalue, such that lambda = alpha/beta and */
/* > mu = beta/alpha. */
/* > */
/* > The vectors x and y in the above equations are right eigenvectors of */
/* > the matrix pair (A,B). Vectors u and v satisfying */
/* > */
/* > u**H*A = lambda*u**H*B or mu*v**H*A = v**H*B */
/* > */
/* > are left eigenvectors of (A,B). */
/* > */
/* > Note: this routine performs "full balancing" on A and B */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVL */
/* > \verbatim */
/* > JOBVL is CHARACTER*1 */
/* > = 'N': do not compute the left generalized eigenvectors;
*/
/* > = 'V': compute the left generalized eigenvectors (returned */
/* > in VL). */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVR */
/* > \verbatim */
/* > JOBVR is CHARACTER*1 */
/* > = 'N': do not compute the right generalized eigenvectors;
*/
/* > = 'V': compute the right generalized eigenvectors (returned */
/* > in VR). */
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
/* > A is REAL array, dimension (LDA, N) */
/* > On entry, the matrix A. */
/* > If JOBVL = 'V' or JOBVR = 'V', then on exit A */
/* > contains the real Schur form of A from the generalized Schur */
/* > factorization of the pair (A,B) after balancing. */
/* > If no eigenvectors were computed, then only the diagonal */
/* > blocks from the Schur form will be correct. See SGGHRD and */
/* > SHGEQZ for details. */
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
/* > B is REAL array, dimension (LDB, N) */
/* > On entry, the matrix B. */
/* > If JOBVL = 'V' or JOBVR = 'V', then on exit B contains the */
/* > upper triangular matrix obtained from B in the generalized */
/* > Schur factorization of the pair (A,B) after balancing. */
/* > If no eigenvectors were computed, then only those elements of */
/* > B corresponding to the diagonal blocks from the Schur form of */
/* > A will be correct. See SGGHRD and SHGEQZ for details. */
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
/* > ALPHAR is REAL array, dimension (N) */
/* > The real parts of each scalar alpha defining an eigenvalue of */
/* > GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* > ALPHAI is REAL array, dimension (N) */
/* > The imaginary parts of each scalar alpha defining an */
/* > eigenvalue of GNEP. If ALPHAI(j) is zero, then the j-th */
/* > eigenvalue is real;
if positive, then the j-th and */
/* > (j+1)-st eigenvalues are a complex conjugate pair, with */
/* > ALPHAI(j+1) = -ALPHAI(j). */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is REAL array, dimension (N) */
/* > The scalars beta that define the eigenvalues of GNEP. */
/* > */
/* > Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and */
/* > beta = BETA(j) represent the j-th eigenvalue of the matrix */
/* > pair (A,B), in one of the forms lambda = alpha/beta or */
/* > mu = beta/alpha. Since either lambda or mu may overflow, */
/* > they should not, in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* > VL is REAL array, dimension (LDVL,N) */
/* > If JOBVL = 'V', the left eigenvectors u(j) are stored */
/* > in the columns of VL, in the same order as their eigenvalues. */
/* > If the j-th eigenvalue is real, then u(j) = VL(:,j). */
/* > If the j-th and (j+1)-st eigenvalues form a complex conjugate */
/* > pair, then */
/* > u(j) = VL(:,j) + i*VL(:,j+1) */
/* > and */
/* > u(j+1) = VL(:,j) - i*VL(:,j+1). */
/* > */
/* > Each eigenvector is scaled so that its largest component has */
/* > f2c_abs(real part) + f2c_abs(imag. part) = 1, except for eigenvectors */
/* > corresponding to an eigenvalue with alpha = beta = 0, which */
/* > are set to zero. */
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
/* > VR is REAL array, dimension (LDVR,N) */
/* > If JOBVR = 'V', the right eigenvectors x(j) are stored */
/* > in the columns of VR, in the same order as their eigenvalues. */
/* > If the j-th eigenvalue is real, then x(j) = VR(:,j). */
/* > If the j-th and (j+1)-st eigenvalues form a complex conjugate */
/* > pair, then */
/* > x(j) = VR(:,j) + i*VR(:,j+1) */
/* > and */
/* > x(j+1) = VR(:,j) - i*VR(:,j+1). */
/* > */
/* > Each eigenvector is scaled so that its largest component has */
/* > f2c_abs(real part) + f2c_abs(imag. part) = 1, except for eigenvalues */
/* > corresponding to an eigenvalue with alpha = beta = 0, which */
/* > are set to zero. */
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
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= max(1,8*N). */
/* > For good performance, LWORK must generally be larger. */
/* > To compute the optimal value of LWORK, call ILAENV to get */
/* > blocksizes (for SGEQRF, SORMQR, and SORGQR.) Then compute: */
/* > NB -- MAX of the blocksizes for SGEQRF, SORMQR, and SORGQR;
*/
/* > The optimal LWORK is: */
/* > 2*N + MAX( 6*N, N*(NB+1) ). */
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
/* > > N: errors that usually indicate LAPACK problems: */
/* > =N+1: error return from SGGBAL */
/* > =N+2: error return from SGEQRF */
/* > =N+3: error return from SORMQR */
/* > =N+4: error return from SORGQR */
/* > =N+5: error return from SGGHRD */
/* > =N+6: error return from SHGEQZ (other than failed */
/* > iteration) */
/* > =N+7: error return from STGEVC */
/* > =N+8: error return from SGGBAK (computing VL) */
/* > =N+9: error return from SGGBAK (computing VR) */
/* > =N+10: error return from SLASCL (various calls) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realGEeigen */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Balancing */
/* > --------- */
/* > */
/* > This driver calls SGGBAL to both permute and scale rows and columns */
/* > of A and B. The permutations PL and PR are chosen so that PL*A*PR */
/* > and PL*B*R will be upper triangular except for the diagonal blocks */
/* > A(i:j,i:j) and B(i:j,i:j), with i and j as close together as */
/* > possible. The diagonal scaling matrices DL and DR are chosen so */
/* > that the pair DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to */
/* > one (except for the elements that start out zero.) */
/* > */
/* > After the eigenvalues and eigenvectors of the balanced matrices */
/* > have been computed, SGGBAK transforms the eigenvectors back to what */
/* > they would have been (in perfect arithmetic) if they had not been */
/* > balanced. */
/* > */
/* > Contents of A and B on Exit */
/* > -------- -- - --- - -- ---- */
/* > */
/* > If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or */
/* > both), then on exit the arrays A and B will contain the real Schur */
/* > form[*] of the "balanced" versions of A and B. If no eigenvectors */
/* > are computed, then only the diagonal blocks will be correct. */
/* > */
/* > [*] See SHGEQZ, SGEGS, or read the book "Matrix Computations", */
/* > by Golub & van Loan, pub. by Johns Hopkins U. Press. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int sgegv_(char *jobvl, char *jobvr, integer *n, real *a, integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real *beta, real *vl, integer *ldvl, real *vr, integer *ldvr, real *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2;
    real r__1, r__2, r__3, r__4;
    /* Local variables */
    integer jc, nb, in, jr, nb1, nb2, nb3, ihi, ilo;
    real eps;
    logical ilv;
    real absb, anrm, bnrm;
    integer itau;
    real temp;
    logical ilvl, ilvr;
    integer lopt;
    real anrm1, anrm2, bnrm1, bnrm2, absai, scale, absar, sbeta;
    extern logical lsame_(char *, char *);
    integer ileft, iinfo, icols, iwork, irows;
    real salfai;
    extern /* Subroutine */
    int sggbak_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, integer * ), sggbal_(char *, integer *, real *, integer *, real *, integer *, integer *, integer *, real *, real *, real *, integer *);
    real salfar;
    extern real slamch_(char *), slange_(char *, integer *, integer *, real *, integer *, real *);
    real safmin;
    extern /* Subroutine */
    int sgghrd_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, integer *, real *, integer * , real *, integer *, integer *);
    real safmax;
    char chtemp[1];
    logical ldumma[1];
    extern /* Subroutine */
    int slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *, integer *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer ijobvl, iright;
    logical ilimit;
    extern /* Subroutine */
    int sgeqrf_(integer *, integer *, real *, integer *, real *, real *, integer *, integer *);
    integer ijobvr;
    extern /* Subroutine */
    int slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *), slaset_(char *, integer *, integer *, real *, real *, real *, integer *), stgevc_( char *, char *, logical *, integer *, real *, integer *, real *, integer *, real *, integer *, real *, integer *, integer *, integer *, real *, integer *);
    real onepls;
    integer lwkmin;
    extern /* Subroutine */
    int shgeqz_(char *, char *, char *, integer *, integer *, integer *, real *, integer *, real *, integer *, real * , real *, real *, real *, integer *, real *, integer *, real *, integer *, integer *), sorgqr_(integer *, integer *, integer *, real *, integer *, real *, real *, integer * , integer *);
    integer lwkopt;
    logical lquery;
    extern /* Subroutine */
    int sormqr_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *);
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
    /* Computing MAX */
    i__1 = *n << 3;
    lwkmin = max(i__1,1);
    lwkopt = lwkmin;
    work[1] = (real) lwkopt;
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
    else if (*ldvl < 1 || ilvl && *ldvl < *n)
    {
        *info = -12;
    }
    else if (*ldvr < 1 || ilvr && *ldvr < *n)
    {
        *info = -14;
    }
    else if (*lwork < lwkmin && ! lquery)
    {
        *info = -16;
    }
    if (*info == 0)
    {
        nb1 = ilaenv_(&c__1, "SGEQRF", " ", n, n, &c_n1, &c_n1);
        nb2 = ilaenv_(&c__1, "SORMQR", " ", n, n, n, &c_n1);
        nb3 = ilaenv_(&c__1, "SORGQR", " ", n, n, n, &c_n1);
        /* Computing MAX */
        i__1 = max(nb1,nb2);
        nb = max(i__1,nb3);
        /* Computing MAX */
        i__1 = *n * 6;
        i__2 = *n * (nb + 1); // , expr subst
        lopt = (*n << 1) + max(i__1,i__2);
        work[1] = (real) lopt;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGEGV ", &i__1);
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
    eps = slamch_("E") * slamch_("B");
    safmin = slamch_("S");
    safmin += safmin;
    safmax = 1.f / safmin;
    onepls = eps * 4 + 1.f;
    /* Scale A */
    anrm = slange_("M", n, n, &a[a_offset], lda, &work[1]);
    anrm1 = anrm;
    anrm2 = 1.f;
    if (anrm < 1.f)
    {
        if (safmax * anrm < 1.f)
        {
            anrm1 = safmin;
            anrm2 = safmax * anrm;
        }
    }
    if (anrm > 0.f)
    {
        slascl_("G", &c_n1, &c_n1, &anrm, &c_b27, n, n, &a[a_offset], lda, & iinfo);
        if (iinfo != 0)
        {
            *info = *n + 10;
            return 0;
        }
    }
    /* Scale B */
    bnrm = slange_("M", n, n, &b[b_offset], ldb, &work[1]);
    bnrm1 = bnrm;
    bnrm2 = 1.f;
    if (bnrm < 1.f)
    {
        if (safmax * bnrm < 1.f)
        {
            bnrm1 = safmin;
            bnrm2 = safmax * bnrm;
        }
    }
    if (bnrm > 0.f)
    {
        slascl_("G", &c_n1, &c_n1, &bnrm, &c_b27, n, n, &b[b_offset], ldb, & iinfo);
        if (iinfo != 0)
        {
            *info = *n + 10;
            return 0;
        }
    }
    /* Permute the matrix to make it more nearly triangular */
    /* Workspace layout: (8*N words -- "work" requires 6*N words) */
    /* left_permutation, right_permutation, work... */
    ileft = 1;
    iright = *n + 1;
    iwork = iright + *n;
    sggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[ ileft], &work[iright], &work[iwork], &iinfo);
    if (iinfo != 0)
    {
        *info = *n + 1;
        goto L120;
    }
    /* Reduce B to triangular form, and initialize VL and/or VR */
    /* Workspace layout: ("work..." must have at least N words) */
    /* left_permutation, right_permutation, tau, work... */
    irows = ihi + 1 - ilo;
    if (ilv)
    {
        icols = *n + 1 - ilo;
    }
    else
    {
        icols = irows;
    }
    itau = iwork;
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;
    sgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[ iwork], &i__1, &iinfo);
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
        goto L120;
    }
    i__1 = *lwork + 1 - iwork;
    sormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, & work[itau], &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, & iinfo);
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
        goto L120;
    }
    if (ilvl)
    {
        slaset_("Full", n, n, &c_b38, &c_b27, &vl[vl_offset], ldvl) ;
        i__1 = irows - 1;
        i__2 = irows - 1;
        slacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb, &vl[ilo + 1 + ilo * vl_dim1], ldvl);
        i__1 = *lwork + 1 - iwork;
        sorgqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[ itau], &work[iwork], &i__1, &iinfo);
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
            goto L120;
        }
    }
    if (ilvr)
    {
        slaset_("Full", n, n, &c_b38, &c_b27, &vr[vr_offset], ldvr) ;
    }
    /* Reduce to generalized Hessenberg form */
    if (ilv)
    {
        /* Eigenvectors requested -- work on whole matrix. */
        sgghrd_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &iinfo);
    }
    else
    {
        sgghrd_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda, &b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[ vr_offset], ldvr, &iinfo);
    }
    if (iinfo != 0)
    {
        *info = *n + 5;
        goto L120;
    }
    /* Perform QZ algorithm */
    /* Workspace layout: ("work..." must have at least 1 word) */
    /* left_permutation, right_permutation, work... */
    iwork = itau;
    if (ilv)
    {
        *(unsigned char *)chtemp = 'S';
    }
    else
    {
        *(unsigned char *)chtemp = 'E';
    }
    i__1 = *lwork + 1 - iwork;
    shgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[ b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[iwork], &i__1, &iinfo);
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
        goto L120;
    }
    if (ilv)
    {
        /* Compute Eigenvectors (STGEVC requires 6*N words of workspace) */
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
        stgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[ iwork], &iinfo);
        if (iinfo != 0)
        {
            *info = *n + 7;
            goto L120;
        }
        /* Undo balancing on VL and VR, rescale */
        if (ilvl)
        {
            sggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, & vl[vl_offset], ldvl, &iinfo);
            if (iinfo != 0)
            {
                *info = *n + 8;
                goto L120;
            }
            i__1 = *n;
            for (jc = 1;
                    jc <= i__1;
                    ++jc)
            {
                if (alphai[jc] < 0.f)
                {
                    goto L50;
                }
                temp = 0.f;
                if (alphai[jc] == 0.f)
                {
                    i__2 = *n;
                    for (jr = 1;
                            jr <= i__2;
                            ++jr)
                    {
                        /* Computing MAX */
                        r__2 = temp;
                        r__3 = (r__1 = vl[jr + jc * vl_dim1], f2c_abs(r__1)); // , expr subst
                        temp = max(r__2,r__3);
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
                        r__3 = temp;
                        r__4 = (r__1 = vl[jr + jc * vl_dim1], f2c_abs(r__1)) + (r__2 = vl[jr + (jc + 1) * vl_dim1], f2c_abs(r__2)); // , expr subst
                        temp = max(r__3,r__4);
                        /* L20: */
                    }
                }
                if (temp < safmin)
                {
                    goto L50;
                }
                temp = 1.f / temp;
                if (alphai[jc] == 0.f)
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
            sggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, & vr[vr_offset], ldvr, &iinfo);
            if (iinfo != 0)
            {
                *info = *n + 9;
                goto L120;
            }
            i__1 = *n;
            for (jc = 1;
                    jc <= i__1;
                    ++jc)
            {
                if (alphai[jc] < 0.f)
                {
                    goto L100;
                }
                temp = 0.f;
                if (alphai[jc] == 0.f)
                {
                    i__2 = *n;
                    for (jr = 1;
                            jr <= i__2;
                            ++jr)
                    {
                        /* Computing MAX */
                        r__2 = temp;
                        r__3 = (r__1 = vr[jr + jc * vr_dim1], f2c_abs(r__1)); // , expr subst
                        temp = max(r__2,r__3);
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
                        r__3 = temp;
                        r__4 = (r__1 = vr[jr + jc * vr_dim1], f2c_abs(r__1)) + (r__2 = vr[jr + (jc + 1) * vr_dim1], f2c_abs(r__2)); // , expr subst
                        temp = max(r__3,r__4);
                        /* L70: */
                    }
                }
                if (temp < safmin)
                {
                    goto L100;
                }
                temp = 1.f / temp;
                if (alphai[jc] == 0.f)
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
    /* Undo scaling in alpha, beta */
    /* Note: this does not give the alpha and beta for the unscaled */
    /* problem. */
    /* Un-scaling is limited to avoid underflow in alpha and beta */
    /* if they are significant. */
    i__1 = *n;
    for (jc = 1;
            jc <= i__1;
            ++jc)
    {
        absar = (r__1 = alphar[jc], f2c_abs(r__1));
        absai = (r__1 = alphai[jc], f2c_abs(r__1));
        absb = (r__1 = beta[jc], f2c_abs(r__1));
        salfar = anrm * alphar[jc];
        salfai = anrm * alphai[jc];
        sbeta = bnrm * beta[jc];
        ilimit = FALSE_;
        scale = 1.f;
        /* Check for significant underflow in ALPHAI */
        /* Computing MAX */
        r__1 = safmin, r__2 = eps * absar;
        r__1 = max(r__1,r__2);
        r__2 = eps * absb; // ; expr subst
        if (f2c_abs(salfai) < safmin && absai >= max(r__1,r__2))
        {
            ilimit = TRUE_;
            /* Computing MAX */
            r__1 = onepls * safmin;
            r__2 = anrm2 * absai; // , expr subst
            scale = onepls * safmin / anrm1 / max(r__1,r__2);
        }
        else if (salfai == 0.f)
        {
            /* If insignificant underflow in ALPHAI, then make the */
            /* conjugate eigenvalue real. */
            if (alphai[jc] < 0.f && jc > 1)
            {
                alphai[jc - 1] = 0.f;
            }
            else if (alphai[jc] > 0.f && jc < *n)
            {
                alphai[jc + 1] = 0.f;
            }
        }
        /* Check for significant underflow in ALPHAR */
        /* Computing MAX */
        r__1 = safmin, r__2 = eps * absai;
        r__1 = max(r__1,r__2);
        r__2 = eps * absb; // ; expr subst
        if (f2c_abs(salfar) < safmin && absar >= max(r__1,r__2))
        {
            ilimit = TRUE_;
            /* Computing MAX */
            /* Computing MAX */
            r__3 = onepls * safmin;
            r__4 = anrm2 * absar; // , expr subst
            r__1 = scale;
            r__2 = onepls * safmin / anrm1 / max(r__3,r__4); // , expr subst
            scale = max(r__1,r__2);
        }
        /* Check for significant underflow in BETA */
        /* Computing MAX */
        r__1 = safmin, r__2 = eps * absar;
        r__1 = max(r__1,r__2);
        r__2 = eps * absai; // ; expr subst
        if (f2c_abs(sbeta) < safmin && absb >= max(r__1,r__2))
        {
            ilimit = TRUE_;
            /* Computing MAX */
            /* Computing MAX */
            r__3 = onepls * safmin;
            r__4 = bnrm2 * absb; // , expr subst
            r__1 = scale;
            r__2 = onepls * safmin / bnrm1 / max(r__3,r__4); // , expr subst
            scale = max(r__1,r__2);
        }
        /* Check for possible overflow when limiting scaling */
        if (ilimit)
        {
            /* Computing MAX */
            r__1 = f2c_abs(salfar), r__2 = f2c_abs(salfai);
            r__1 = max(r__1,r__2);
            r__2 = f2c_abs(sbeta); // ; expr subst
            temp = scale * safmin * max(r__1,r__2);
            if (temp > 1.f)
            {
                scale /= temp;
            }
            if (scale < 1.f)
            {
                ilimit = FALSE_;
            }
        }
        /* Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary. */
        if (ilimit)
        {
            salfar = scale * alphar[jc] * anrm;
            salfai = scale * alphai[jc] * anrm;
            sbeta = scale * beta[jc] * bnrm;
        }
        alphar[jc] = salfar;
        alphai[jc] = salfai;
        beta[jc] = sbeta;
        /* L110: */
    }
L120:
    work[1] = (real) lwkopt;
    return 0;
    /* End of SGEGV */
}
/* sgegv_ */
