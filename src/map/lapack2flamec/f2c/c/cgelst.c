/* cgelst.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    0.f,0.f
}
;
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__0 = 0;
/* > \brief <b> CGELST solves overdetermined or underdetermined systems for GE matrices using QR or LQ factori zation with compact WY representation of Q.</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGELST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgelst. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgelst. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgelst. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGELST( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, */
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
/* > CGELST solves overdetermined or underdetermined real linear systems */
/* > involving an M-by-N matrix A, or its conjugate-transpose, using a QR */
/* > or LQ factorization of A with compact WY representation of Q. */
/* > It is assumed that A has full rank. */
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
/* > an underdetermined system A**T * X = B. */
/* > */
/* > 4. If TRANS = 'C' and m < n: find the least squares solution of */
/* > an overdetermined system, i.e., solve the least squares problem */
/* > minimize || B - A**T * X ||. */
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
/* > columns of the matrices B and X. NRHS >=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, */
/* > if M >= N, A is overwritten by details of its QR */
/* > factorization as returned by CGEQRT;
*/
/* > if M < N, A is overwritten by details of its LQ */
/* > factorization as returned by CGELQT. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
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
/* > solution in each column is given by the sum of squares of */
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
/* > The leading dimension of the array B. LDB >= fla_max(1,M,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (fla_max(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > LWORK >= fla_max( 1, MN + fla_max( MN, NRHS ) ). */
/* > For optimal performance, */
/* > LWORK >= fla_max( 1, (MN + fla_max( MN, NRHS ))*NB ). */
/* > where MN = fla_min(M,N) and NB is the optimum block size. */
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
/* > \ingroup complexGEsolve */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > November 2022, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int cgelst_(char *trans, integer *m, integer *n, integer * nrhs, complex *a, integer *lda, complex *b, integer *ldb, complex * work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgelst inputs: trans %c, m %" FLA_IS ", n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "",*trans, *m, *n, *nrhs, *lda, *ldb);
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
    integer nbmin;
    real rwork[1];
    integer lwopt;
    extern /* Subroutine */
    int slabad_(real *, real *);
    extern real clange_(char *, integer *, integer *, complex *, integer *, real *);
    extern /* Subroutine */
    int clascl_(char *, integer *, integer *, real *, real *, integer *, integer *, complex *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int cgelqt_(integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *);
    integer scllen;
    real bignum;
    extern /* Subroutine */
    int cgeqrt_(integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *);
    integer mnnrhs;
    real smlnum;
    logical lquery;
    extern /* Subroutine */
    int ctrtrs_(char *, char *, char *, integer *, integer *, complex *, integer *, complex *, integer *, integer *), cgemlqt_(char *, char *, integer *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *), cgemqrt_(char *, char *, integer *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *);
    /* -- LAPACK driver routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
    mn = fla_min(*m,*n);
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
    else if (*lda < fla_max(1,*m))
    {
        *info = -6;
    }
    else /* if(complicated condition) */
    {
        /* Computing fla_max */
        i__1 = fla_max(1,*m);
        if (*ldb < fla_max(i__1,*n))
        {
            *info = -8;
        }
        else /* if(complicated condition) */
        {
            /* Computing fla_max */
            i__1 = 1;
            i__2 = mn + fla_max(mn,*nrhs); // , expr subst
            if (*lwork < fla_max(i__1,i__2) && ! lquery)
            {
                *info = -10;
            }
        }
    }
    /* Figure out optimal block size and optimal workspace size */
    if (*info == 0 || *info == -10)
    {
        tpsd = TRUE_;
        if (lsame_(trans, "N"))
        {
            tpsd = FALSE_;
        }
        nb = ilaenv_(&c__1, "CGELST", " ", m, n, &c_n1, &c_n1);
        mnnrhs = fla_max(mn,*nrhs);
        /* Computing fla_max */
        i__1 = 1;
        i__2 = (mn + mnnrhs) * nb; // , expr subst
        lwopt = fla_max(i__1,i__2);
        r__1 = (real) lwopt;
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGELST ", &i__1);
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    else if (lquery)
    {
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Quick return if possible */
    /* Computing fla_min */
    i__1 = fla_min(*m,*n);
    if (fla_min(i__1,*nrhs) == 0)
    {
        i__1 = fla_max(*m,*n);
        claset_("Full", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
        r__1 = (real) lwopt;
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* *GEQRT and *GELQT routines cannot accept NB larger than fla_min(M,N) */
    if (nb > mn)
    {
        nb = mn;
    }
    /* Determine the block size from the supplied LWORK */
    /* ( at this stage we know that LWORK >= (minimum required workspace, */
    /* but it may be less than optimal) */
    /* Computing fla_min */
    i__1 = nb;
    i__2 = *lwork / (mn + mnnrhs); // , expr subst
    nb = fla_min(i__1,i__2);
    /* The minimum value of NB, when blocked code is used */
    /* Computing fla_max */
    i__1 = 2;
    i__2 = ilaenv_(&c__2, "CGELST", " ", m, n, &c_n1, &c_n1); // , expr subst
    nbmin = fla_max(i__1,i__2);
    if (nb < nbmin)
    {
        nb = 1;
    }
    /* Get machine parameters */
    smlnum = slamch_("S") / slamch_("P");
    bignum = 1.f / smlnum;
    slabad_(&smlnum, &bignum);
    /* Scale A, B if fla_max element outside range [SMLNUM,BIGNUM] */
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
        i__1 = fla_max(*m,*n);
        claset_("Full", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
        r__1 = (real) lwopt;
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
        return 0;
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
        /* M > N: */
        /* Compute the blocked QR factorization of A, */
        /* using the compact WY representation of Q, */
        /* workspace at least N, optimally N*NB. */
        cgeqrt_(m, n, &nb, &a[a_offset], lda, &work[1], &nb, &work[mn * nb + 1], info);
        if (! tpsd)
        {
            /* M > N, A is not transposed: */
            /* Overdetermined system of equations, */
            /* least-squares problem, fla_min || A * X - B ||. */
            /* Compute B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS), */
            /* using the compact WY representation of Q, */
            /* workspace at least NRHS, optimally NRHS*NB. */
            cgemqrt_("Left", "Conjugate transpose", m, nrhs, n, &nb, &a[ a_offset], lda, &work[1], &nb, &b[b_offset], ldb, &work[ mn * nb + 1], info);
            /* Compute B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */
            ctrtrs_("Upper", "No transpose", "Non-unit", n, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
    AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
            scllen = *n;
        }
        else
        {
            /* M > N, A is transposed: */
            /* Underdetermined system of equations, */
            /* minimum norm solution of A**T * X = B. */
            /* Compute B := inv(R**T) * B in two row blocks of B. */
            /* Block 1: B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS) */
            ctrtrs_("Upper", "Conjugate transpose", "Non-unit", n, nrhs, &a[ a_offset], lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
    AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
            /* Block 2: Zero out all rows below the N-th row in B: */
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
                }
            }
            /* Compute B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS), */
            /* using the compact WY representation of Q, */
            /* workspace at least NRHS, optimally NRHS*NB. */
            cgemqrt_("Left", "No transpose", m, nrhs, n, &nb, &a[a_offset], lda, &work[1], &nb, &b[b_offset], ldb, &work[mn * nb + 1], info);
            scllen = *m;
        }
    }
    else
    {
        /* M < N: */
        /* Compute the blocked LQ factorization of A, */
        /* using the compact WY representation of Q, */
        /* workspace at least M, optimally M*NB. */
        cgelqt_(m, n, &nb, &a[a_offset], lda, &work[1], &nb, &work[mn * nb + 1], info);
        if (! tpsd)
        {
            /* M < N, A is not transposed: */
            /* Underdetermined system of equations, */
            /* minimum norm solution of A * X = B. */
            /* Compute B := inv(L) * B in two row blocks of B. */
            /* Block 1: B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */
            ctrtrs_("Lower", "No transpose", "Non-unit", m, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
    AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
            /* Block 2: Zero out all rows below the M-th row in B: */
            /* B(M+1:N,1:NRHS) = ZERO */
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
                }
            }
            /* Compute B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS), */
            /* using the compact WY representation of Q, */
            /* workspace at least NRHS, optimally NRHS*NB. */
            cgemlqt_("Left", "Conjugate transpose", n, nrhs, m, &nb, &a[ a_offset], lda, &work[1], &nb, &b[b_offset], ldb, &work[ mn * nb + 1], info);
            scllen = *n;
        }
        else
        {
            /* M < N, A is transposed: */
            /* Overdetermined system of equations, */
            /* least-squares problem, fla_min || A**T * X - B ||. */
            /* Compute B(1:N,1:NRHS) := Q * B(1:N,1:NRHS), */
            /* using the compact WY representation of Q, */
            /* workspace at least NRHS, optimally NRHS*NB. */
            cgemlqt_("Left", "No transpose", n, nrhs, m, &nb, &a[a_offset], lda, &work[1], &nb, &b[b_offset], ldb, &work[mn * nb + 1], info);
            /* Compute B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS) */
            ctrtrs_("Lower", "Conjugate transpose", "Non-unit", m, nrhs, &a[ a_offset], lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
    AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
            scllen = *m;
        }
    }
    /* Undo scaling */
    if (iascl == 1)
    {
        clascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset], ldb, info);
    }
    else if (iascl == 2)
    {
        clascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset], ldb, info);
    }
    if (ibscl == 1)
    {
        clascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset], ldb, info);
    }
    else if (ibscl == 2)
    {
        clascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset], ldb, info);
    }
    r__1 = (real) lwopt;
    work[1].r = r__1;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of CGELST */
}
/* cgelst_ */
