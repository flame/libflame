/* ../netlib/cgels.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    0.f,0.f
}
;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
/* > \brief <b> CGELS solves overdetermined or underdetermined systems for GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGELS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgels.f "> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgels.f "> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgels.f "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, LDA, LDB, LWORK, M, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGELS solves overdetermined or underdetermined complex linear systems */
/* > involving an M-by-N matrix A, or its conjugate-transpose, using a QR */
/* > or LQ factorization of A. It is assumed that A has full rank. */
/* > */
/* > The following options are provided: */
/* > */
/* > 1. If TRANS = 'N' and m >= n: find the least squares solution of */
/* > an overdetermined system, i.e., solve the least squares problem */
/* > minimize || B - A*X ||. */
/* > */
/* > 2. If TRANS = 'N' and m < n: find the minimum norm solution of */
/* > an underdetermined system A * X = B. */
/* > */
/* > 3. If TRANS = 'C' and m >= n: find the minimum norm solution of */
/* > an undetermined system A**H * X = B. */
/* > */
/* > 4. If TRANS = 'C' and m < n: find the least squares solution of */
/* > an overdetermined system, i.e., solve the least squares problem */
/* > minimize || B - A**H * X ||. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call;
they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution */
/* > matrix X. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': the linear system involves A;
*/
/* > = 'C': the linear system involves A**H. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of */
/* > columns of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > if M >= N, A is overwritten by details of its QR */
/* > factorization as returned by CGEQRF;
*/
/* > if M < N, A is overwritten by details of its LQ */
/* > factorization as returned by CGELQF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > On entry, the matrix B of right hand side vectors, stored */
/* > columnwise;
B is M-by-NRHS if TRANS = 'N', or N-by-NRHS */
/* > if TRANS = 'C'. */
/* > On exit, if INFO = 0, B is overwritten by the solution */
/* > vectors, stored columnwise: */
/* > if TRANS = 'N' and m >= n, rows 1 to n of B contain the least */
/* > squares solution vectors;
the residual sum of squares for the */
/* > solution in each column is given by the sum of squares of the */
/* > modulus of elements N+1 to M in that column;
*/
/* > if TRANS = 'N' and m < n, rows 1 to N of B contain the */
/* > minimum norm solution vectors;
*/
/* > if TRANS = 'C' and m >= n, rows 1 to M of B contain the */
/* > minimum norm solution vectors;
*/
/* > if TRANS = 'C' and m < n, rows 1 to M of B contain the */
/* > least squares solution vectors;
the residual sum of squares */
/* > for the solution in each column is given by the sum of */
/* > squares of the modulus of elements M+1 to N in that column. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= MAX(1,M,N). */
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
/* > The dimension of the array WORK. */
/* > LWORK >= max( 1, MN + max( MN, NRHS ) ). */
/* > For optimal performance, */
/* > LWORK >= max( 1, MN + max( MN, NRHS )*NB ). */
/* > where MN = min(M,N) and NB is the optimum block size. */
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
/* > > 0: if INFO = i, the i-th diagonal element of the */
/* > triangular factor of A is zero, so that A does not have */
/* > full rank;
the least squares solution could not be */
/* > computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexGEsolve */
/* ===================================================================== */
/* Subroutine */
int cgels_(char *trans, integer *m, integer *n, integer * nrhs, complex *a, integer *lda, complex *b, integer *ldb, complex * work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    real r__1;
    /* Local variables */
    integer i__, j, nb, mn;
    real anrm, bnrm;
    integer brow;
    logical tpsd;
    integer iascl, ibscl;
    extern logical lsame_(char *, char *);
    integer wsize;
    real rwork[1];
    extern /* Subroutine */
    int slabad_(real *, real *);
    extern real clange_(char *, integer *, integer *, complex *, integer *, real *);
    extern /* Subroutine */
    int cgelqf_(integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), clascl_( char *, integer *, integer *, real *, real *, integer *, integer * , complex *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int cgeqrf_(integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), claset_( char *, integer *, integer *, complex *, complex *, complex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer scllen;
    real bignum;
    extern /* Subroutine */
    int cunmlq_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *), cunmqr_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *);
    real smlnum;
    logical lquery;
    extern /* Subroutine */
    int ctrtrs_(char *, char *, char *, integer *, integer *, complex *, integer *, complex *, integer *, integer *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    /* Function Body */
    *info = 0;
    mn = min(*m,*n);
    lquery = *lwork == -1;
    if (! (lsame_(trans, "N") || lsame_(trans, "C")))
    {
        *info = -1;
    }
    else if (*m < 0)
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*nrhs < 0)
    {
        *info = -4;
    }
    else if (*lda < max(1,*m))
    {
        *info = -6;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = max(1,*m);
        if (*ldb < max(i__1,*n))
        {
            *info = -8;
        }
        else /* if(complicated condition) */
        {
            /* Computing MAX */
            i__1 = 1;
            i__2 = mn + max(mn,*nrhs); // , expr subst
            if (*lwork < max(i__1,i__2) && ! lquery)
            {
                *info = -10;
            }
        }
    }
    /* Figure out optimal block size */
    if (*info == 0 || *info == -10)
    {
        tpsd = TRUE_;
        if (lsame_(trans, "N"))
        {
            tpsd = FALSE_;
        }
        if (*m >= *n)
        {
            nb = ilaenv_(&c__1, "CGEQRF", " ", m, n, &c_n1, &c_n1);
            if (tpsd)
            {
                /* Computing MAX */
                i__1 = nb;
                i__2 = ilaenv_(&c__1, "CUNMQR", "LN", m, nrhs, n, & c_n1); // , expr subst
                nb = max(i__1,i__2);
            }
            else
            {
                /* Computing MAX */
                i__1 = nb;
                i__2 = ilaenv_(&c__1, "CUNMQR", "LC", m, nrhs, n, & c_n1); // , expr subst
                nb = max(i__1,i__2);
            }
        }
        else
        {
            nb = ilaenv_(&c__1, "CGELQF", " ", m, n, &c_n1, &c_n1);
            if (tpsd)
            {
                /* Computing MAX */
                i__1 = nb;
                i__2 = ilaenv_(&c__1, "CUNMLQ", "LC", n, nrhs, m, & c_n1); // , expr subst
                nb = max(i__1,i__2);
            }
            else
            {
                /* Computing MAX */
                i__1 = nb;
                i__2 = ilaenv_(&c__1, "CUNMLQ", "LN", n, nrhs, m, & c_n1); // , expr subst
                nb = max(i__1,i__2);
            }
        }
        /* Computing MAX */
        i__1 = 1;
        i__2 = mn + max(mn,*nrhs) * nb; // , expr subst
        wsize = max(i__1,i__2);
        r__1 = (real) wsize;
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGELS ", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    /* Computing MIN */
    i__1 = min(*m,*n);
    if (min(i__1,*nrhs) == 0)
    {
        i__1 = max(*m,*n);
        claset_("Full", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
        return 0;
    }
    /* Get machine parameters */
    smlnum = slamch_("S") / slamch_("P");
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);
    /* Scale A, B if max element outside range [SMLNUM,BIGNUM] */
    anrm = clange_("M", m, n, &a[a_offset], lda, rwork);
    iascl = 0;
    if (anrm > 0.f && anrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        clascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, info);
        iascl = 1;
    }
    else if (anrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        clascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, info);
        iascl = 2;
    }
    else if (anrm == 0.f)
    {
        /* Matrix all zero. Return zero solution. */
        i__1 = max(*m,*n);
        claset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
        goto L50;
    }
    brow = *m;
    if (tpsd)
    {
        brow = *n;
    }
    bnrm = clange_("M", &brow, nrhs, &b[b_offset], ldb, rwork);
    ibscl = 0;
    if (bnrm > 0.f && bnrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        clascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], ldb, info);
        ibscl = 1;
    }
    else if (bnrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        clascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], ldb, info);
        ibscl = 2;
    }
    if (*m >= *n)
    {
        /* compute QR factorization of A */
        i__1 = *lwork - mn;
        cgeqrf_(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info) ;
        /* workspace at least N, optimally N*NB */
        if (! tpsd)
        {
            /* Least-Squares Problem min || A * X - B || */
            /* B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS) */
            i__1 = *lwork - mn;
            cunmqr_("Left", "Conjugate transpose", m, nrhs, n, &a[a_offset], lda, &work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info);
            /* workspace at least NRHS, optimally NRHS*NB */
            /* B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */
            ctrtrs_("Upper", "No transpose", "Non-unit", n, nrhs, &a[a_offset] , lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
                return 0;
            }
            scllen = *n;
        }
        else
        {
            /* Overdetermined system of equations A**H * X = B */
            /* B(1:N,1:NRHS) := inv(R**H) * B(1:N,1:NRHS) */
            ctrtrs_("Upper", "Conjugate transpose", "Non-unit", n, nrhs, &a[ a_offset], lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
                return 0;
            }
            /* B(N+1:M,1:NRHS) = ZERO */
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = *m;
                for (i__ = *n + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    b[i__3].r = 0.f;
                    b[i__3].i = 0.f; // , expr subst
                    /* L10: */
                }
                /* L20: */
            }
            /* B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */
            i__1 = *lwork - mn;
            cunmqr_("Left", "No transpose", m, nrhs, n, &a[a_offset], lda, & work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info);
            /* workspace at least NRHS, optimally NRHS*NB */
            scllen = *m;
        }
    }
    else
    {
        /* Compute LQ factorization of A */
        i__1 = *lwork - mn;
        cgelqf_(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info) ;
        /* workspace at least M, optimally M*NB. */
        if (! tpsd)
        {
            /* underdetermined system of equations A * X = B */
            /* B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */
            ctrtrs_("Lower", "No transpose", "Non-unit", m, nrhs, &a[a_offset] , lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
                return 0;
            }
            /* B(M+1:N,1:NRHS) = 0 */
            i__1 = *nrhs;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = *n;
                for (i__ = *m + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    b[i__3].r = 0.f;
                    b[i__3].i = 0.f; // , expr subst
                    /* L30: */
                }
                /* L40: */
            }
            /* B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS) */
            i__1 = *lwork - mn;
            cunmlq_("Left", "Conjugate transpose", n, nrhs, m, &a[a_offset], lda, &work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info);
            /* workspace at least NRHS, optimally NRHS*NB */
            scllen = *n;
        }
        else
        {
            /* overdetermined system min || A**H * X - B || */
            /* B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */
            i__1 = *lwork - mn;
            cunmlq_("Left", "No transpose", n, nrhs, m, &a[a_offset], lda, & work[1], &b[b_offset], ldb, &work[mn + 1], &i__1, info);
            /* workspace at least NRHS, optimally NRHS*NB */
            /* B(1:M,1:NRHS) := inv(L**H) * B(1:M,1:NRHS) */
            ctrtrs_("Lower", "Conjugate transpose", "Non-unit", m, nrhs, &a[ a_offset], lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
                return 0;
            }
            scllen = *m;
        }
    }
    /* Undo scaling */
    if (iascl == 1)
    {
        clascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset] , ldb, info);
    }
    else if (iascl == 2)
    {
        clascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset] , ldb, info);
    }
    if (ibscl == 1)
    {
        clascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset] , ldb, info);
    }
    else if (ibscl == 2)
    {
        clascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset] , ldb, info);
    }
L50:
    r__1 = (real) wsize;
    work[1].r = r__1;
    work[1].i = 0.f; // , expr subst
    return 0;
    /* End of CGELS */
}
/* cgels_ */
