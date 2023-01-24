/* slaqz4.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b4 = 0.f;
static real c_b5 = 1.f;
static integer c__1 = 1;
static logical c_true = TRUE_;
/* > \brief \b SLAQZ4 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAQZ4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqz4. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqz4. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqz4. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAQZ4( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NSHIFTS, */
/* $ NBLOCK_DESIRED, SR, SI, SS, A, LDA, B, LDB, Q, LDQ, Z, LDZ, */
/* $ QC, LDQC, ZC, LDZC, WORK, LWORK, INFO ) */
/* IMPLICIT NONE */
/* Function arguments */
/* LOGICAL, INTENT( IN ) :: ILSCHUR, ILQ, ILZ */
/* INTEGER, INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, */
/* $ NSHIFTS, NBLOCK_DESIRED, LDQC, LDZC */
/* REAL, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ Z( LDZ, * ), QC( LDQC, * ), ZC( LDZC, * ), WORK( * ), SR( * ), */
/* $ SI( * ), SS( * ) */
/* INTEGER, INTENT( OUT ) :: INFO */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAQZ4 Executes a single multishift QZ sweep */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ILSCHUR */
/* > \verbatim */
/* > ILSCHUR is LOGICAL */
/* > Determines whether or not to update the full Schur form */
/* > \endverbatim */
/* > */
/* > \param[in] ILQ */
/* > \verbatim */
/* > ILQ is LOGICAL */
/* > Determines whether or not to update the matrix Q */
/* > \endverbatim */
/* > */
/* > \param[in] ILZ */
/* > \verbatim */
/* > ILZ is LOGICAL */
/* > Determines whether or not to update the matrix Z */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A, B, Q, and Z. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] NSHIFTS */
/* > \verbatim */
/* > NSHIFTS is INTEGER */
/* > The desired number of shifts to use */
/* > \endverbatim */
/* > */
/* > \param[in] NBLOCK_DESIRED */
/* > \verbatim */
/* > NBLOCK_DESIRED is INTEGER */
/* > The desired size of the computational windows */
/* > \endverbatim */
/* > */
/* > \param[in] SR */
/* > \verbatim */
/* > SR is REAL array. SR contains */
/* > the real parts of the shifts to use. */
/* > \endverbatim */
/* > */
/* > \param[in] SI */
/* > \verbatim */
/* > SI is REAL array. SI contains */
/* > the imaginary parts of the shifts to use. */
/* > \endverbatim */
/* > */
/* > \param[in] SS */
/* > \verbatim */
/* > SS is REAL array. SS contains */
/* > the scale of the shifts to use. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (LDQ, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] QC */
/* > \verbatim */
/* > QC is REAL array, dimension (LDQC, NBLOCK_DESIRED) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQC */
/* > \verbatim */
/* > LDQC is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] ZC */
/* > \verbatim */
/* > ZC is REAL array, dimension (LDZC, NBLOCK_DESIRED) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZC */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO >= 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,N). */
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
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Thijs Steel, KU Leuven */
/* > \date May 2020 */
/* > \ingroup doubleGEcomputational */
/* > */
/* ===================================================================== */
/* Subroutine */
int slaqz4_(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi, integer *nshifts, integer * nblock_desired__, real *sr, real *si, real *ss, real *a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real *z__, integer * ldz, real *qc, integer *ldqc, real *zc, integer *ldzc, real *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slaqz4 inputs: n %" FLA_IS ", ilo %" FLA_IS ", ihi %" FLA_IS ", nshifts %" FLA_IS ", nblock_desired__ %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS ", ldq %" FLA_IS ", ldz %" FLA_IS ", ldqc %" FLA_IS ", ldzc %" FLA_IS "",*n, *ilo, *ihi, *nshifts, *nblock_desired__, *lda, *ldb, *ldq, *ldz, *ldqc, *ldzc);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, qc_dim1, qc_offset, zc_dim1, zc_offset, i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
    integer i__, j, k;
    real v[3], c1, c2, s1, s2;
    integer np, ns;
    real temp, swap;
    integer npos;
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *), sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *), slaqz1_(real *, integer *, real *, integer *, real *, real *, real *, real *, real *, real *), slaqz2_(logical *, logical *, integer *, integer *, integer *, integer *, real *, integer *, real *, integer *, integer *, integer *, real *, integer *, integer *, integer *, real *, integer *);
    integer nblock;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    integer ishift;
    extern /* Subroutine */
    int slaset_(char *, integer *, integer *, real *, real *, real *, integer *), slartg_(real *, real *, real *, real *, real *), slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *);
    integer istopb, swidth, istopm, sheight, istartb, istartm;
    /* Function arguments */
    /* Parameters */
    /* Local scalars */
    /* External functions */
    /* Parameter adjustments */
    --sr;
    --si;
    --ss;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    qc_dim1 = *ldqc;
    qc_offset = 1 + qc_dim1;
    qc -= qc_offset;
    zc_dim1 = *ldzc;
    zc_offset = 1 + zc_dim1;
    zc -= zc_offset;
    --work;
    /* Function Body */
    *info = 0;
    if (*nblock_desired__ < *nshifts + 1)
    {
        *info = -8;
    }
    if (*lwork == -1)
    {
        /* workspace query, quick return */
        work[1] = (real) (*n * *nblock_desired__);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (*lwork < *n * *nblock_desired__)
    {
        *info = -25;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SLAQZ4", &i__1);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Executable statements */
    if (*nshifts < 2)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (*ilo >= *ihi)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (*ilschur)
    {
        istartm = 1;
        istopm = *n;
    }
    else
    {
        istartm = *ilo;
        istopm = *ihi;
    }
    /* Shuffle shifts into pairs of real shifts and pairs */
    /* of complex conjugate shifts assuming complex */
    /* conjugate shifts are already adjacent to one */
    /* another */
    i__1 = *nshifts - 2;
    for (i__ = 1;
            i__ <= i__1;
            i__ += 2)
    {
        if (si[i__] != -si[i__ + 1])
        {
            swap = sr[i__];
            sr[i__] = sr[i__ + 1];
            sr[i__ + 1] = sr[i__ + 2];
            sr[i__ + 2] = swap;
            swap = si[i__];
            si[i__] = si[i__ + 1];
            si[i__ + 1] = si[i__ + 2];
            si[i__ + 2] = swap;
            swap = ss[i__];
            ss[i__] = ss[i__ + 1];
            ss[i__ + 1] = ss[i__ + 2];
            ss[i__ + 2] = swap;
        }
    }
    /* NSHFTS is supposed to be even, but if it is odd, */
    /* then simply reduce it by one. The shuffle above */
    /* ensures that the dropped shift is real and that */
    /* the remaining shifts are paired. */
    ns = *nshifts - *nshifts % 2;
    /* Computing MAX */
    i__1 = *nblock_desired__ - ns;
    npos = fla_max(i__1,1);
    /* The following block introduces the shifts and chases */
    /* them down one by one just enough to make space for */
    /* the other shifts. The near-the-diagonal block is */
    /* of size (ns+1) x ns. */
    i__1 = ns + 1;
    i__2 = ns + 1;
    slaset_("FULL", &i__1, &i__2, &c_b4, &c_b5, &qc[qc_offset], ldqc);
    slaset_("FULL", &ns, &ns, &c_b4, &c_b5, &zc[zc_offset], ldzc);
    i__1 = ns;
    for (i__ = 1;
            i__ <= i__1;
            i__ += 2)
    {
        /* Introduce the shift */
        slaqz1_(&a[*ilo + *ilo * a_dim1], lda, &b[*ilo + *ilo * b_dim1], ldb, &sr[i__], &sr[i__ + 1], &si[i__], &ss[i__], &ss[i__ + 1], v);
        temp = v[1];
        slartg_(&temp, &v[2], &c1, &s1, &v[1]);
        slartg_(v, &v[1], &c2, &s2, &temp);
        srot_(&ns, &a[*ilo + 1 + *ilo * a_dim1], lda, &a[*ilo + 2 + *ilo * a_dim1], lda, &c1, &s1);
        srot_(&ns, &a[*ilo + *ilo * a_dim1], lda, &a[*ilo + 1 + *ilo * a_dim1], lda, &c2, &s2);
        srot_(&ns, &b[*ilo + 1 + *ilo * b_dim1], ldb, &b[*ilo + 2 + *ilo * b_dim1], ldb, &c1, &s1);
        srot_(&ns, &b[*ilo + *ilo * b_dim1], ldb, &b[*ilo + 1 + *ilo * b_dim1], ldb, &c2, &s2);
        i__2 = ns + 1;
        srot_(&i__2, &qc[(qc_dim1 << 1) + 1], &c__1, &qc[qc_dim1 * 3 + 1], & c__1, &c1, &s1);
        i__2 = ns + 1;
        srot_(&i__2, &qc[qc_dim1 + 1], &c__1, &qc[(qc_dim1 << 1) + 1], &c__1, &c2, &s2);
        /* Chase the shift down */
        i__2 = ns - 1 - i__;
        for (j = 1;
                j <= i__2;
                ++j)
        {
            i__3 = *ihi - *ilo + 1;
            i__4 = ns + 1;
            slaqz2_(&c_true, &c_true, &j, &c__1, &ns, &i__3, &a[*ilo + *ilo * a_dim1], lda, &b[*ilo + *ilo * b_dim1], ldb, &i__4, &c__1, &qc[qc_offset], ldqc, &ns, &c__1, &zc[zc_offset], ldzc);
        }
    }
    /* Update the rest of the pencil */
    /* Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm) */
    /* from the left with Qc(1:ns+1,1:ns+1)' */
    sheight = ns + 1;
    swidth = istopm - (*ilo + ns) + 1;
    if (swidth > 0)
    {
        sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[qc_offset], ldqc, &a[*ilo + (*ilo + ns) * a_dim1], lda, &c_b4, &work[1], & sheight);
        slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[*ilo + (*ilo + ns) * a_dim1], lda);
        sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[qc_offset], ldqc, &b[*ilo + (*ilo + ns) * b_dim1], ldb, &c_b4, &work[1], & sheight);
        slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[*ilo + (*ilo + ns) * b_dim1], ldb);
    }
    if (*ilq)
    {
        sgemm_("N", "N", n, &sheight, &sheight, &c_b5, &q[*ilo * q_dim1 + 1], ldq, &qc[qc_offset], ldqc, &c_b4, &work[1], n);
        slacpy_("ALL", n, &sheight, &work[1], n, &q[*ilo * q_dim1 + 1], ldq);
    }
    /* Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1) */
    /* from the right with Zc(1:ns,1:ns) */
    sheight = *ilo - 1 - istartm + 1;
    swidth = ns;
    if (sheight > 0)
    {
        sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &a[istartm + *ilo * a_dim1], lda, &zc[zc_offset], ldzc, &c_b4, &work[1], & sheight);
        slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm + * ilo * a_dim1], lda);
        sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &b[istartm + *ilo * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b4, &work[1], & sheight);
        slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm + * ilo * b_dim1], ldb);
    }
    if (*ilz)
    {
        sgemm_("N", "N", n, &swidth, &swidth, &c_b5, &z__[*ilo * z_dim1 + 1], ldz, &zc[zc_offset], ldzc, &c_b4, &work[1], n);
        slacpy_("ALL", n, &swidth, &work[1], n, &z__[*ilo * z_dim1 + 1], ldz);
    }
    /* The following block chases the shifts down to the bottom */
    /* right block. If possible, a shift is moved down npos */
    /* positions at a time */
    k = *ilo;
    while(k < *ihi - ns)
    {
        /* Computing MIN */
        i__1 = *ihi - ns - k;
        np = fla_min(i__1,npos);
        /* Size of the near-the-diagonal block */
        nblock = ns + np;
        /* istartb points to the first row we will be updating */
        istartb = k + 1;
        /* istopb points to the last column we will be updating */
        istopb = k + nblock - 1;
        i__1 = ns + np;
        i__2 = ns + np;
        slaset_("FULL", &i__1, &i__2, &c_b4, &c_b5, &qc[qc_offset], ldqc);
        i__1 = ns + np;
        i__2 = ns + np;
        slaset_("FULL", &i__1, &i__2, &c_b4, &c_b5, &zc[zc_offset], ldzc);
        /* Near the diagonal shift chase */
        for (i__ = ns - 1;
                i__ >= 0;
                i__ += -2)
        {
            i__1 = np - 1;
            for (j = 0;
                    j <= i__1;
                    ++j)
            {
                /* Move down the block with index k+i+j-1, updating */
                /* the (ns+np x ns+np) block: */
                /* (k:k+ns+np,k:k+ns+np-1) */
                i__2 = k + i__ + j - 1;
                i__3 = k + 1;
                slaqz2_(&c_true, &c_true, &i__2, &istartb, &istopb, ihi, &a[ a_offset], lda, &b[b_offset], ldb, &nblock, &i__3, & qc[qc_offset], ldqc, &nblock, &k, &zc[zc_offset], ldzc);
            }
        }
        /* Update rest of the pencil */
        /* Update A(k+1:k+ns+np, k+ns+np:istopm) and */
        /* B(k+1:k+ns+np, k+ns+np:istopm) */
        /* from the left with Qc(1:ns+np,1:ns+np)' */
        sheight = ns + np;
        swidth = istopm - (k + ns + np) + 1;
        if (swidth > 0)
        {
            sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[ qc_offset], ldqc, &a[k + 1 + (k + ns + np) * a_dim1], lda, &c_b4, &work[1], &sheight);
            slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[k + 1 + ( k + ns + np) * a_dim1], lda);
            sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[ qc_offset], ldqc, &b[k + 1 + (k + ns + np) * b_dim1], ldb, &c_b4, &work[1], &sheight);
            slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[k + 1 + ( k + ns + np) * b_dim1], ldb);
        }
        if (*ilq)
        {
            sgemm_("N", "N", n, &nblock, &nblock, &c_b5, &q[(k + 1) * q_dim1 + 1], ldq, &qc[qc_offset], ldqc, &c_b4, &work[1], n);
            slacpy_("ALL", n, &nblock, &work[1], n, &q[(k + 1) * q_dim1 + 1], ldq);
        }
        /* Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1) */
        /* from the right with Zc(1:ns+np,1:ns+np) */
        sheight = k - istartm + 1;
        swidth = nblock;
        if (sheight > 0)
        {
            sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &a[istartm + k * a_dim1], lda, &zc[zc_offset], ldzc, &c_b4, &work[1], & sheight);
            slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm + k * a_dim1], lda);
            sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &b[istartm + k * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b4, &work[1], & sheight);
            slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm + k * b_dim1], ldb);
        }
        if (*ilz)
        {
            sgemm_("N", "N", n, &nblock, &nblock, &c_b5, &z__[k * z_dim1 + 1], ldz, &zc[zc_offset], ldzc, &c_b4, &work[1], n);
            slacpy_("ALL", n, &nblock, &work[1], n, &z__[k * z_dim1 + 1], ldz);
        }
        k += np;
    }
    /* The following block removes the shifts from the bottom right corner */
    /* one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi). */
    slaset_("FULL", &ns, &ns, &c_b4, &c_b5, &qc[qc_offset], ldqc);
    i__1 = ns + 1;
    i__2 = ns + 1;
    slaset_("FULL", &i__1, &i__2, &c_b4, &c_b5, &zc[zc_offset], ldzc);
    /* istartb points to the first row we will be updating */
    istartb = *ihi - ns + 1;
    /* istopb points to the last column we will be updating */
    istopb = *ihi;
    i__1 = ns;
    for (i__ = 1;
            i__ <= i__1;
            i__ += 2)
    {
        /* Chase the shift down to the bottom right corner */
        i__2 = *ihi - 2;
        for (ishift = *ihi - i__ - 1;
                ishift <= i__2;
                ++ishift)
        {
            i__3 = *ihi - ns + 1;
            i__4 = ns + 1;
            i__5 = *ihi - ns;
            slaqz2_(&c_true, &c_true, &ishift, &istartb, &istopb, ihi, &a[ a_offset], lda, &b[b_offset], ldb, &ns, &i__3, &qc[ qc_offset], ldqc, &i__4, &i__5, &zc[zc_offset], ldzc);
        }
    }
    /* Update rest of the pencil */
    /* Update A(ihi-ns+1:ihi, ihi+1:istopm) */
    /* from the left with Qc(1:ns,1:ns)' */
    sheight = ns;
    swidth = istopm - (*ihi + 1) + 1;
    if (swidth > 0)
    {
        sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[qc_offset], ldqc, &a[*ihi - ns + 1 + (*ihi + 1) * a_dim1], lda, &c_b4, & work[1], &sheight);
        slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[*ihi - ns + 1 + (*ihi + 1) * a_dim1], lda);
        sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[qc_offset], ldqc, &b[*ihi - ns + 1 + (*ihi + 1) * b_dim1], ldb, &c_b4, & work[1], &sheight);
        slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[*ihi - ns + 1 + (*ihi + 1) * b_dim1], ldb);
    }
    if (*ilq)
    {
        sgemm_("N", "N", n, &ns, &ns, &c_b5, &q[(*ihi - ns + 1) * q_dim1 + 1], ldq, &qc[qc_offset], ldqc, &c_b4, &work[1], n);
        slacpy_("ALL", n, &ns, &work[1], n, &q[(*ihi - ns + 1) * q_dim1 + 1], ldq);
    }
    /* Update A(istartm:ihi-ns,ihi-ns:ihi) */
    /* from the right with Zc(1:ns+1,1:ns+1) */
    sheight = *ihi - ns - istartm + 1;
    swidth = ns + 1;
    if (sheight > 0)
    {
        sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &a[istartm + (* ihi - ns) * a_dim1], lda, &zc[zc_offset], ldzc, &c_b4, &work[ 1], &sheight);
        slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm + (* ihi - ns) * a_dim1], lda);
        sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &b[istartm + (* ihi - ns) * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b4, &work[ 1], &sheight);
        slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm + (* ihi - ns) * b_dim1], ldb);
    }
    if (*ilz)
    {
        i__1 = ns + 1;
        i__2 = ns + 1;
        sgemm_("N", "N", n, &i__1, &i__2, &c_b5, &z__[(*ihi - ns) * z_dim1 + 1], ldz, &zc[zc_offset], ldzc, &c_b4, &work[1], n);
        i__1 = ns + 1;
        slacpy_("ALL", n, &i__1, &work[1], n, &z__[(*ihi - ns) * z_dim1 + 1], ldz);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
}
/* slaqz4_ */
