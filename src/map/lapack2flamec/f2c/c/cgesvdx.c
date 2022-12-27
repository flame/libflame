/* ../netlib/v3.9.0/cgesvdx.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    0.f,0.f
}
;
static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief <b> CGESVDX computes the singular value decomposition (SVD) for GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGESVDX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesvdx .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesvdx .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesvdx .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, */
/* $ IL, IU, NS, S, U, LDU, VT, LDVT, WORK, */
/* $ LWORK, RWORK, IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBU, JOBVT, RANGE */
/* INTEGER IL, INFO, IU, LDA, LDU, LDVT, LWORK, M, N, NS */
/* REAL VL, VU */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL S( * ), RWORK( * ) */
/* COMPLEX A( LDA, * ), U( LDU, * ), VT( LDVT, * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGESVDX computes the singular value decomposition (SVD) of a complex */
/* > M-by-N matrix A, optionally computing the left and/or right singular */
/* > vectors. The SVD is written */
/* > */
/* > A = U * SIGMA * transpose(V) */
/* > */
/* > where SIGMA is an M-by-N matrix which is zero except for its */
/* > fla_min(m,n) diagonal elements, U is an M-by-M unitary matrix, and */
/* > V is an N-by-N unitary matrix. The diagonal elements of SIGMA */
/* > are the singular values of A;
they are real and non-negative, and */
/* > are returned in descending order. The first fla_min(m,n) columns of */
/* > U and V are the left and right singular vectors of A. */
/* > */
/* > CGESVDX uses an eigenvalue problem for obtaining the SVD, which */
/* > allows for the computation of a subset of singular values and */
/* > vectors. See SBDSVDX for details. */
/* > */
/* > Note that the routine returns V**T, not V. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBU */
/* > \verbatim */
/* > JOBU is CHARACTER*1 */
/* > Specifies options for computing all or part of the matrix U: */
/* > = 'V': the first fla_min(m,n) columns of U (the left singular */
/* > vectors) or as specified by RANGE are returned in */
/* > the array U;
*/
/* > = 'N': no columns of U (no left singular vectors) are */
/* > computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVT */
/* > \verbatim */
/* > JOBVT is CHARACTER*1 */
/* > Specifies options for computing all or part of the matrix */
/* > V**T: */
/* > = 'V': the first fla_min(m,n) rows of V**T (the right singular */
/* > vectors) or as specified by RANGE are returned in */
/* > the array VT;
*/
/* > = 'N': no rows of V**T (no right singular vectors) are */
/* > computed. */
/* > \endverbatim */
/* > */
/* > \param[in] RANGE */
/* > \verbatim */
/* > RANGE is CHARACTER*1 */
/* > = 'A': all singular values will be found. */
/* > = 'V': all singular values in the half-open interval (VL,VU] */
/* > will be found. */
/* > = 'I': the IL-th through IU-th singular values will be found. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the input matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the input matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the contents of A are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* > VL is REAL */
/* > If RANGE='V', the lower bound of the interval to */
/* > be searched for singular values. VU > VL. */
/* > Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* > VU is REAL */
/* > If RANGE='V', the upper bound of the interval to */
/* > be searched for singular values. VU > VL. */
/* > Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* > IL is INTEGER */
/* > If RANGE='I', the index of the */
/* > smallest singular value to be returned. */
/* > 1 <= IL <= IU <= fla_min(M,N), if fla_min(M,N) > 0. */
/* > Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* > IU is INTEGER */
/* > If RANGE='I', the index of the */
/* > largest singular value to be returned. */
/* > 1 <= IL <= IU <= fla_min(M,N), if fla_min(M,N) > 0. */
/* > Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[out] NS */
/* > \verbatim */
/* > NS is INTEGER */
/* > The total number of singular values found, */
/* > 0 <= NS <= fla_min(M,N). */
/* > If RANGE = 'A', NS = fla_min(M,N);
if RANGE = 'I', NS = IU-IL+1. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is REAL array, dimension (fla_min(M,N)) */
/* > The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* > U is COMPLEX array, dimension (LDU,UCOL) */
/* > If JOBU = 'V', U contains columns of U (the left singular */
/* > vectors, stored columnwise) as specified by RANGE;
if */
/* > JOBU = 'N', U is not referenced. */
/* > Note: The user must ensure that UCOL >= NS;
if RANGE = 'V', */
/* > the exact value of NS is not known in advance and an upper */
/* > bound must be used. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is INTEGER */
/* > The leading dimension of the array U. LDU >= 1;
if */
/* > JOBU = 'V', LDU >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* > VT is COMPLEX array, dimension (LDVT,N) */
/* > If JOBVT = 'V', VT contains the rows of V**T (the right singular */
/* > vectors, stored rowwise) as specified by RANGE;
if JOBVT = 'N', */
/* > VT is not referenced. */
/* > Note: The user must ensure that LDVT >= NS;
if RANGE = 'V', */
/* > the exact value of NS is not known in advance and an upper */
/* > bound must be used. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* > LDVT is INTEGER */
/* > The leading dimension of the array VT. LDVT >= 1;
if */
/* > JOBVT = 'V', LDVT >= NS (see above). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
*/
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > LWORK >= MAX(1,MIN(M,N)*(MIN(M,N)+4)) for the paths (see */
/* > comments inside the code): */
/* > - PATH 1 (M much larger than N) */
/* > - PATH 1t (N much larger than M) */
/* > LWORK >= MAX(1,MIN(M,N)*2+MAX(M,N)) for the other paths. */
/* > For good performance, LWORK should generally be larger. */
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
/* > RWORK is REAL array, dimension (MAX(1,LRWORK)) */
/* > LRWORK >= MIN(M,N)*(MIN(M,N)*2+15*MIN(M,N)). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (12*MIN(M,N)) */
/* > If INFO = 0, the first NS elements of IWORK are zero. If INFO > 0, */
/* > then IWORK contains the indices of the eigenvectors that failed */
/* > to converge in SBDSVDX/SSTEVX. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, then i eigenvectors failed to converge */
/* > in SBDSVDX/SSTEVX. */
/* > if INFO = N*2 + 1, an internal error occurred in */
/* > SBDSVDX */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2016 */
/* > \ingroup complexGEsing */
/* ===================================================================== */
/* Subroutine */
int cgesvdx_(char *jobu, char *jobvt, char *range, integer * m, integer *n, complex *a, integer *lda, real *vl, real *vu, integer * il, integer *iu, integer *ns, real *s, complex *u, integer *ldu, complex *vt, integer *ldvt, complex *work, integer *lwork, real * rwork, integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"cgesvdx inputs: jobu %c, jobvt %c, range %c, m %lld, n %lld, lda %lld, il %lld, iu %lld, ldu %lld, ldvt %lld, lwork %lld",*jobu, *jobvt, *range, *m, *n, *lda, *il, *iu, *ldu, *ldvt, *lwork);
#else
    snprintf(buffer, 256,"cgesvdx inputs: jobu %c, jobvt %c, range %c, m %d, n %d, lda %d, il %d, iu %d, ldu %d, ldvt %d, lwork %d",*jobu, *jobvt, *range, *m, *n, *lda, *il, *iu, *ldu, *ldvt, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1[2], i__2, i__3, i__4, i__5;
    real r__1;
    complex q__1;
    char ch__1[2];
    /* Builtin functions */
    /* Subroutine */

    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k, id, ie;
    real dum[1], eps;
    integer iscl;
    logical alls, inds;
    integer ilqf;
    real anrm;
    integer ierr, iqrf, itau;
    char jobz[1];
    logical vals;
    extern logical lsame_(char *, char *);
    integer iltgk, itemp, minmn, itaup, itauq, iutgk, itgkz, mnthr;
    logical wantu;
    extern /* Subroutine */
    int cgebrd_(integer *, integer *, complex *, integer *, real *, real *, complex *, complex *, complex *, integer *, integer *);
    extern real clange_(char *, integer *, integer *, complex *, integer *, real *);
    extern /* Subroutine */
    int cgelqf_(integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), clascl_( char *, integer *, integer *, real *, real *, integer *, integer *, complex *, integer *, integer *), cgeqrf_(integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *), clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    real bignum;
    extern /* Subroutine */
    int slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *, integer *, integer *);
    real abstol;
    extern /* Subroutine */
    int cunmbr_(char *, char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *);
    char rngtgk[1];
    extern /* Subroutine */
    int cunmlq_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *);
    integer itempr;
    extern /* Subroutine */
    int cunmqr_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *);
    integer minwrk, maxwrk;
    real smlnum;
    logical lquery, wantvt;
    extern /* Subroutine */
    int sbdsvdx_(char *, char *, char *, integer *, real *, real *, real *, real *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, integer *);
    /* -- LAPACK driver routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2016 */
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
    /* Test the input arguments. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --work;
    --rwork;
    --iwork;
    /* Function Body */
    *ns = 0;
    *info = 0;
    abstol = slamch_("S") * 2;
    lquery = *lwork == -1;
    minmn = fla_min(*m,*n);
    wantu = lsame_(jobu, "V");
    wantvt = lsame_(jobvt, "V");
    if (wantu || wantvt)
    {
        *(unsigned char *)jobz = 'V';
    }
    else
    {
        *(unsigned char *)jobz = 'N';
    }
    alls = lsame_(range, "A");
    vals = lsame_(range, "V");
    inds = lsame_(range, "I");
    *info = 0;
    if (! lsame_(jobu, "V") && ! lsame_(jobu, "N"))
    {
        *info = -1;
    }
    else if (! lsame_(jobvt, "V") && ! lsame_(jobvt, "N"))
    {
        *info = -2;
    }
    else if (! (alls || vals || inds))
    {
        *info = -3;
    }
    else if (*m < 0)
    {
        *info = -4;
    }
    else if (*n < 0)
    {
        *info = -5;
    }
    else if (*m > *lda)
    {
        *info = -7;
    }
    else if (minmn > 0)
    {
        if (vals)
        {
            if (*vl < 0.f)
            {
                *info = -8;
            }
            else if (*vu <= *vl)
            {
                *info = -9;
            }
        }
        else if (inds)
        {
            if (*il < 1 || *il > fla_max(1,minmn))
            {
                *info = -10;
            }
            else if (*iu < fla_min(minmn,*il) || *iu > minmn)
            {
                *info = -11;
            }
        }
        if (*info == 0)
        {
            if (wantu && *ldu < *m)
            {
                *info = -15;
            }
            else if (wantvt)
            {
                if (inds)
                {
                    if (*ldvt < *iu - *il + 1)
                    {
                        *info = -17;
                    }
                }
                else if (*ldvt < minmn)
                {
                    *info = -17;
                }
            }
        }
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* NB refers to the optimal block size for the immediately */
    /* following subroutine, as returned by ILAENV.) */
    if (*info == 0)
    {
        minwrk = 1;
        maxwrk = 1;
        if (minmn > 0)
        {
            if (*m >= *n)
            {
                mnthr = ilaenv_(&c__6, "CGESVD", ch__1, m, n, &c__0, &c__0);
                if (*m >= mnthr)
                {
                    /* Path 1 (M much larger than N) */
                    minwrk = *n * (*n + 5);
                    maxwrk = *n + *n * ilaenv_(&c__1, "CGEQRF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = *n * *n + (*n << 1) + (*n << 1) * ilaenv_(&c__1, "CGEBRD", " ", n, n, &c_n1, &c_n1); // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    if (wantu || wantvt)
                    {
                        /* Computing MAX */
                        i__2 = maxwrk;
                        i__3 = *n * *n + (*n << 1) + *n * ilaenv_(&c__1, "CUNMQR", "LN", n, n, n, &c_n1); // , expr subst
                        maxwrk = fla_max(i__2,i__3);
                    }
                }
                else
                {
                    /* Path 2 (M at least N, but not much larger) */
                    minwrk = *n * 3 + *m;
                    maxwrk = (*n << 1) + (*m + *n) * ilaenv_(&c__1, "CGEBRD", " ", m, n, &c_n1, &c_n1);
                    if (wantu || wantvt)
                    {
                        /* Computing MAX */
                        i__2 = maxwrk;
                        i__3 = (*n << 1) + *n * ilaenv_(&c__1, "CUNMQR", "LN", n, n, n, &c_n1); // , expr subst
                        maxwrk = fla_max(i__2,i__3);
                    }
                }
            }
            else
            {
                mnthr = ilaenv_(&c__6, "CGESVD", ch__1, m, n, &c__0, &c__0);
                if (*n >= mnthr)
                {
                    /* Path 1t (N much larger than M) */
                    minwrk = *m * (*m + 5);
                    maxwrk = *m + *m * ilaenv_(&c__1, "CGELQF", " ", m, n, & c_n1, &c_n1);
                    /* Computing MAX */
                    i__2 = maxwrk;
                    i__3 = *m * *m + (*m << 1) + (*m << 1) * ilaenv_(&c__1, "CGEBRD", " ", m, m, &c_n1, &c_n1); // , expr subst
                    maxwrk = fla_max(i__2,i__3);
                    if (wantu || wantvt)
                    {
                        /* Computing MAX */
                        i__2 = maxwrk;
                        i__3 = *m * *m + (*m << 1) + *m * ilaenv_(&c__1, "CUNMQR", "LN", m, m, m, &c_n1); // , expr subst
                        maxwrk = fla_max(i__2,i__3);
                    }
                }
                else
                {
                    /* Path 2t (N greater than M, but not much larger) */
                    minwrk = *m * 3 + *n;
                    maxwrk = (*m << 1) + (*m + *n) * ilaenv_(&c__1, "CGEBRD", " ", m, n, &c_n1, &c_n1);
                    if (wantu || wantvt)
                    {
                        /* Computing MAX */
                        i__2 = maxwrk;
                        i__3 = (*m << 1) + *m * ilaenv_(&c__1, "CUNMQR", "LN", m, m, m, &c_n1); // , expr subst
                        maxwrk = fla_max(i__2,i__3);
                    }
                }
            }
        }
        maxwrk = fla_max(maxwrk,minwrk);
        r__1 = (real) maxwrk;
        q__1.r = r__1;
        q__1.i = 0.f; // , expr subst
        work[1].r = q__1.r;
        work[1].i = q__1.i; // , expr subst
        if (*lwork < minwrk && ! lquery)
        {
            *info = -19;
        }
    }
    if (*info != 0)
    {
        i__2 = -(*info);
        xerbla_("CGESVDX", &i__2);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Set singular values indices accord to RANGE='A'. */
    if (alls)
    {
        *(unsigned char *)rngtgk = 'I';
        iltgk = 1;
        iutgk = fla_min(*m,*n);
    }
    else if (inds)
    {
        *(unsigned char *)rngtgk = 'I';
        iltgk = *il;
        iutgk = *iu;
    }
    else
    {
        *(unsigned char *)rngtgk = 'V';
        iltgk = 0;
        iutgk = 0;
    }
    /* Get machine constants */
    eps = slamch_("P");
    smlnum = sqrt(slamch_("S")) / eps;
    bignum = 1.f / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = clange_("M", m, n, &a[a_offset], lda, dum);
    iscl = 0;
    if (anrm > 0.f && anrm < smlnum)
    {
        iscl = 1;
        clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, info);
    }
    else if (anrm > bignum)
    {
        iscl = 1;
        clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, info);
    }
    if (*m >= *n)
    {
        /* A has at least as many rows as columns. If A has sufficiently */
        /* more rows than columns, first reduce A using the QR */
        /* decomposition. */
        if (*m >= mnthr)
        {
            /* Path 1 (M much larger than N): */
            /* A = Q * R = Q * ( QB * B * PB**T ) */
            /* = Q * ( QB * ( UB * S * VB**T ) * PB**T ) */
            /* U = Q * QB * UB;
            V**T = VB**T * PB**T */
            /* Compute A=Q*R */
            /* (Workspace: need 2*N, prefer N+N*NB) */
            itau = 1;
            itemp = itau + *n;
            i__2 = *lwork - itemp + 1;
            cgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2, info);
            /* Copy R into WORK and bidiagonalize it: */
            /* (Workspace: need N*N+3*N, prefer N*N+N+2*N*NB) */
            iqrf = itemp;
            itauq = itemp + *n * *n;
            itaup = itauq + *n;
            itemp = itaup + *n;
            id = 1;
            ie = id + *n;
            itgkz = ie + *n;
            clacpy_("U", n, n, &a[a_offset], lda, &work[iqrf], n);
            i__2 = *n - 1;
            i__3 = *n - 1;
            claset_("L", &i__2, &i__3, &c_b1, &c_b1, &work[iqrf + 1], n);
            i__2 = *lwork - itemp + 1;
            cgebrd_(n, n, &work[iqrf], n, &rwork[id], &rwork[ie], &work[itauq], &work[itaup], &work[itemp], &i__2, info);
            itempr = itgkz + *n * ((*n << 1) + 1);
            /* Solve eigenvalue problem TGK*Z=Z*S. */
            /* (Workspace: need 2*N*N+14*N) */
            i__2 = *n << 1;
            sbdsvdx_("U", jobz, rngtgk, n, &rwork[id], &rwork[ie], vl, vu, & iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[ itempr], &iwork[1], info) ;
            /* If needed, compute left singular vectors. */
            if (wantu)
            {
                k = itgkz;
                i__2 = *ns;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = *n;
                    for (j = 1;
                            j <= i__3;
                            ++j)
                    {
                        i__4 = j + i__ * u_dim1;
                        i__5 = k;
                        q__1.r = rwork[i__5];
                        q__1.i = 0.f; // , expr subst
                        u[i__4].r = q__1.r;
                        u[i__4].i = q__1.i; // , expr subst
                        ++k;
                    }
                    k += *n;
                }
                i__2 = *m - *n;
                claset_("A", &i__2, ns, &c_b1, &c_b1, &u[*n + 1 + u_dim1], ldu);
                /* Call CUNMBR to compute QB*UB. */
                /* (Workspace in WORK( ITEMP ): need N, prefer N*NB) */
                i__2 = *lwork - itemp + 1;
                cunmbr_("Q", "L", "N", n, ns, n, &work[iqrf], n, &work[itauq], &u[u_offset], ldu, &work[itemp], &i__2, info);
                /* Call CUNMQR to compute Q*(QB*UB). */
                /* (Workspace in WORK( ITEMP ): need N, prefer N*NB) */
                i__2 = *lwork - itemp + 1;
                cunmqr_("L", "N", m, ns, n, &a[a_offset], lda, &work[itau], & u[u_offset], ldu, &work[itemp], &i__2, info);
            }
            /* If needed, compute right singular vectors. */
            if (wantvt)
            {
                k = itgkz + *n;
                i__2 = *ns;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = *n;
                    for (j = 1;
                            j <= i__3;
                            ++j)
                    {
                        i__4 = i__ + j * vt_dim1;
                        i__5 = k;
                        q__1.r = rwork[i__5];
                        q__1.i = 0.f; // , expr subst
                        vt[i__4].r = q__1.r;
                        vt[i__4].i = q__1.i; // , expr subst
                        ++k;
                    }
                    k += *n;
                }
                /* Call CUNMBR to compute VB**T * PB**T */
                /* (Workspace in WORK( ITEMP ): need N, prefer N*NB) */
                i__2 = *lwork - itemp + 1;
                cunmbr_("P", "R", "C", ns, n, n, &work[iqrf], n, &work[itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, info);
            }
        }
        else
        {
            /* Path 2 (M at least N, but not much larger) */
            /* Reduce A to bidiagonal form without QR decomposition */
            /* A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
            /* U = QB * UB;
            V**T = VB**T * PB**T */
            /* Bidiagonalize A */
            /* (Workspace: need 2*N+M, prefer 2*N+(M+N)*NB) */
            itauq = 1;
            itaup = itauq + *n;
            itemp = itaup + *n;
            id = 1;
            ie = id + *n;
            itgkz = ie + *n;
            i__2 = *lwork - itemp + 1;
            cgebrd_(m, n, &a[a_offset], lda, &rwork[id], &rwork[ie], &work[ itauq], &work[itaup], &work[itemp], &i__2, info);
            itempr = itgkz + *n * ((*n << 1) + 1);
            /* Solve eigenvalue problem TGK*Z=Z*S. */
            /* (Workspace: need 2*N*N+14*N) */
            i__2 = *n << 1;
            sbdsvdx_("U", jobz, rngtgk, n, &rwork[id], &rwork[ie], vl, vu, & iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[ itempr], &iwork[1], info) ;
            /* If needed, compute left singular vectors. */
            if (wantu)
            {
                k = itgkz;
                i__2 = *ns;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = *n;
                    for (j = 1;
                            j <= i__3;
                            ++j)
                    {
                        i__4 = j + i__ * u_dim1;
                        i__5 = k;
                        q__1.r = rwork[i__5];
                        q__1.i = 0.f; // , expr subst
                        u[i__4].r = q__1.r;
                        u[i__4].i = q__1.i; // , expr subst
                        ++k;
                    }
                    k += *n;
                }
                i__2 = *m - *n;
                claset_("A", &i__2, ns, &c_b1, &c_b1, &u[*n + 1 + u_dim1], ldu);
                /* Call CUNMBR to compute QB*UB. */
                /* (Workspace in WORK( ITEMP ): need N, prefer N*NB) */
                i__2 = *lwork - itemp + 1;
                cunmbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[ itauq], &u[u_offset], ldu, &work[itemp], &i__2, &ierr);
            }
            /* If needed, compute right singular vectors. */
            if (wantvt)
            {
                k = itgkz + *n;
                i__2 = *ns;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = *n;
                    for (j = 1;
                            j <= i__3;
                            ++j)
                    {
                        i__4 = i__ + j * vt_dim1;
                        i__5 = k;
                        q__1.r = rwork[i__5];
                        q__1.i = 0.f; // , expr subst
                        vt[i__4].r = q__1.r;
                        vt[i__4].i = q__1.i; // , expr subst
                        ++k;
                    }
                    k += *n;
                }
                /* Call CUNMBR to compute VB**T * PB**T */
                /* (Workspace in WORK( ITEMP ): need N, prefer N*NB) */
                i__2 = *lwork - itemp + 1;
                cunmbr_("P", "R", "C", ns, n, n, &a[a_offset], lda, &work[ itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, & ierr);
            }
        }
    }
    else
    {
        /* A has more columns than rows. If A has sufficiently more */
        /* columns than rows, first reduce A using the LQ decomposition. */
        if (*n >= mnthr)
        {
            /* Path 1t (N much larger than M): */
            /* A = L * Q = ( QB * B * PB**T ) * Q */
            /* = ( QB * ( UB * S * VB**T ) * PB**T ) * Q */
            /* U = QB * UB ;
            V**T = VB**T * PB**T * Q */
            /* Compute A=L*Q */
            /* (Workspace: need 2*M, prefer M+M*NB) */
            itau = 1;
            itemp = itau + *m;
            i__2 = *lwork - itemp + 1;
            cgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[itemp], &i__2, info);
            /* Copy L into WORK and bidiagonalize it: */
            /* (Workspace in WORK( ITEMP ): need M*M+3*M, prefer M*M+M+2*M*NB) */
            ilqf = itemp;
            itauq = ilqf + *m * *m;
            itaup = itauq + *m;
            itemp = itaup + *m;
            id = 1;
            ie = id + *m;
            itgkz = ie + *m;
            clacpy_("L", m, m, &a[a_offset], lda, &work[ilqf], m);
            i__2 = *m - 1;
            i__3 = *m - 1;
            claset_("U", &i__2, &i__3, &c_b1, &c_b1, &work[ilqf + *m], m);
            i__2 = *lwork - itemp + 1;
            cgebrd_(m, m, &work[ilqf], m, &rwork[id], &rwork[ie], &work[itauq], &work[itaup], &work[itemp], &i__2, info);
            itempr = itgkz + *m * ((*m << 1) + 1);
            /* Solve eigenvalue problem TGK*Z=Z*S. */
            /* (Workspace: need 2*M*M+14*M) */
            i__2 = *m << 1;
            sbdsvdx_("U", jobz, rngtgk, m, &rwork[id], &rwork[ie], vl, vu, & iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[ itempr], &iwork[1], info) ;
            /* If needed, compute left singular vectors. */
            if (wantu)
            {
                k = itgkz;
                i__2 = *ns;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = *m;
                    for (j = 1;
                            j <= i__3;
                            ++j)
                    {
                        i__4 = j + i__ * u_dim1;
                        i__5 = k;
                        q__1.r = rwork[i__5];
                        q__1.i = 0.f; // , expr subst
                        u[i__4].r = q__1.r;
                        u[i__4].i = q__1.i; // , expr subst
                        ++k;
                    }
                    k += *m;
                }
                /* Call CUNMBR to compute QB*UB. */
                /* (Workspace in WORK( ITEMP ): need M, prefer M*NB) */
                i__2 = *lwork - itemp + 1;
                cunmbr_("Q", "L", "N", m, ns, m, &work[ilqf], m, &work[itauq], &u[u_offset], ldu, &work[itemp], &i__2, info);
            }
            /* If needed, compute right singular vectors. */
            if (wantvt)
            {
                k = itgkz + *m;
                i__2 = *ns;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = *m;
                    for (j = 1;
                            j <= i__3;
                            ++j)
                    {
                        i__4 = i__ + j * vt_dim1;
                        i__5 = k;
                        q__1.r = rwork[i__5];
                        q__1.i = 0.f; // , expr subst
                        vt[i__4].r = q__1.r;
                        vt[i__4].i = q__1.i; // , expr subst
                        ++k;
                    }
                    k += *m;
                }
                i__2 = *n - *m;
                claset_("A", ns, &i__2, &c_b1, &c_b1, &vt[(*m + 1) * vt_dim1 + 1], ldvt);
                /* Call CUNMBR to compute (VB**T)*(PB**T) */
                /* (Workspace in WORK( ITEMP ): need M, prefer M*NB) */
                i__2 = *lwork - itemp + 1;
                cunmbr_("P", "R", "C", ns, m, m, &work[ilqf], m, &work[itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, info);
                /* Call CUNMLQ to compute ((VB**T)*(PB**T))*Q. */
                /* (Workspace in WORK( ITEMP ): need M, prefer M*NB) */
                i__2 = *lwork - itemp + 1;
                cunmlq_("R", "N", ns, n, m, &a[a_offset], lda, &work[itau], & vt[vt_offset], ldvt, &work[itemp], &i__2, info);
            }
        }
        else
        {
            /* Path 2t (N greater than M, but not much larger) */
            /* Reduce to bidiagonal form without LQ decomposition */
            /* A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T */
            /* U = QB * UB;
            V**T = VB**T * PB**T */
            /* Bidiagonalize A */
            /* (Workspace: need 2*M+N, prefer 2*M+(M+N)*NB) */
            itauq = 1;
            itaup = itauq + *m;
            itemp = itaup + *m;
            id = 1;
            ie = id + *m;
            itgkz = ie + *m;
            i__2 = *lwork - itemp + 1;
            cgebrd_(m, n, &a[a_offset], lda, &rwork[id], &rwork[ie], &work[ itauq], &work[itaup], &work[itemp], &i__2, info);
            itempr = itgkz + *m * ((*m << 1) + 1);
            /* Solve eigenvalue problem TGK*Z=Z*S. */
            /* (Workspace: need 2*M*M+14*M) */
            i__2 = *m << 1;
            sbdsvdx_("L", jobz, rngtgk, m, &rwork[id], &rwork[ie], vl, vu, & iltgk, &iutgk, ns, &s[1], &rwork[itgkz], &i__2, &rwork[ itempr], &iwork[1], info) ;
            /* If needed, compute left singular vectors. */
            if (wantu)
            {
                k = itgkz;
                i__2 = *ns;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = *m;
                    for (j = 1;
                            j <= i__3;
                            ++j)
                    {
                        i__4 = j + i__ * u_dim1;
                        i__5 = k;
                        q__1.r = rwork[i__5];
                        q__1.i = 0.f; // , expr subst
                        u[i__4].r = q__1.r;
                        u[i__4].i = q__1.i; // , expr subst
                        ++k;
                    }
                    k += *m;
                }
                /* Call CUNMBR to compute QB*UB. */
                /* (Workspace in WORK( ITEMP ): need M, prefer M*NB) */
                i__2 = *lwork - itemp + 1;
                cunmbr_("Q", "L", "N", m, ns, n, &a[a_offset], lda, &work[ itauq], &u[u_offset], ldu, &work[itemp], &i__2, info);
            }
            /* If needed, compute right singular vectors. */
            if (wantvt)
            {
                k = itgkz + *m;
                i__2 = *ns;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = *m;
                    for (j = 1;
                            j <= i__3;
                            ++j)
                    {
                        i__4 = i__ + j * vt_dim1;
                        i__5 = k;
                        q__1.r = rwork[i__5];
                        q__1.i = 0.f; // , expr subst
                        vt[i__4].r = q__1.r;
                        vt[i__4].i = q__1.i; // , expr subst
                        ++k;
                    }
                    k += *m;
                }
                i__2 = *n - *m;
                claset_("A", ns, &i__2, &c_b1, &c_b1, &vt[(*m + 1) * vt_dim1 + 1], ldvt);
                /* Call CUNMBR to compute VB**T * PB**T */
                /* (Workspace in WORK( ITEMP ): need M, prefer M*NB) */
                i__2 = *lwork - itemp + 1;
                cunmbr_("P", "R", "C", ns, n, m, &a[a_offset], lda, &work[ itaup], &vt[vt_offset], ldvt, &work[itemp], &i__2, info);
            }
        }
    }
    /* Undo scaling if necessary */
    if (iscl == 1)
    {
        if (anrm > bignum)
        {
            slascl_("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], & minmn, info);
        }
        if (anrm < smlnum)
        {
            slascl_("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], & minmn, info);
        }
    }
    /* Return optimal workspace in WORK(1) */
    r__1 = (real) maxwrk;
    q__1.r = r__1;
    q__1.i = 0.f; // , expr subst
    work[1].r = q__1.r;
    work[1].i = q__1.i; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CGESVDX */
}
/* cgesvdx_ */
