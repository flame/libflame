/* zgetsls.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    0.,0.
}
;
static integer c_n1 = -1;
static integer c_n2 = -2;
static integer c__0 = 0;
/* > \brief \b ZGETSLS */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, */
/* $ WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, LDA, LDB, LWORK, M, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGETSLS solves overdetermined or underdetermined complex linear systems */
/* > involving an M-by-N matrix A, using a tall skinny QR or short wide LQ */
/* > factorization of A. It is assumed that A has full rank. */
/* > */
/* > */
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
/* > an undetermined system A**T * X = B. */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, */
/* > A is overwritten by details of its QR or LQ */
/* > factorization as returned by ZGEQR or ZGELQ. */
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
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > On entry, the matrix B of right hand side vectors, stored */
/* > columnwise;
B is M-by-NRHS if TRANS = 'N', or N-by-NRHS */
/* > if TRANS = 'C'. */
/* > On exit, if INFO = 0, B is overwritten by the solution */
/* > vectors, stored columnwise: */
/* > if TRANS = 'N' and m >= n, rows 1 to n of B contain the least */
/* > squares solution vectors. */
/* > if TRANS = 'N' and m < n, rows 1 to N of B contain the */
/* > minimum norm solution vectors;
*/
/* > if TRANS = 'C' and m >= n, rows 1 to M of B contain the */
/* > minimum norm solution vectors;
*/
/* > if TRANS = 'C' and m < n, rows 1 to M of B contain the */
/* > least squares solution vectors. */
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
/* > (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) contains optimal (or either minimal */
/* > or optimal, if query was assumed) LWORK. */
/* > See LWORK for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If LWORK = -1 or -2, then a workspace query is assumed. */
/* > If LWORK = -1, the routine calculates optimal size of WORK for the */
/* > optimal performance and returns this value in WORK(1). */
/* > If LWORK = -2, the routine calculates minimal size of WORK and */
/* > returns this value in WORK(1). */
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
/* > \ingroup complex16GEsolve */
/* ===================================================================== */
/* Subroutine */
int zgetsls_(char *trans, integer *m, integer *n, integer * nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgetsls inputs: trans %c, m %" FLA_IS ", n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "",*trans, *m, *n, *nrhs, *lda, *ldb);

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1;
    /* Local variables */
    integer i__, j;
    doublecomplex tq[5];
    integer lw1, lw2;
    doublereal dum[1];
    integer lwm, lwo;
    doublereal anrm, bnrm;
    logical tran;
    integer brow, tszm, tszo, info2, iascl, ibscl;
    extern logical lsame_(char *, char *);
    integer maxmn;
    extern /* Subroutine */
    int zgelq_(integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *), zgeqr_(integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    doublecomplex workq[1];
    extern /* Subroutine */
    int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    integer scllen;
    doublereal bignum;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, integer *, doublereal *);
    extern /* Subroutine */
    int zlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, doublecomplex *, integer *, integer *), zgemlq_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *), zlaset_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, integer *), zgemqr_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    doublereal smlnum;
    integer wsizem, wsizeo;
    logical lquery;
    extern /* Subroutine */
    int ztrtrs_(char *, char *, char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
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
    maxmn = fla_max(*m,*n);
    tran = lsame_(trans, "C");
    lquery = *lwork == -1 || *lwork == -2;
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
        /* Computing MAX */
        i__1 = fla_max(1,*m);
        if (*ldb < fla_max(i__1,*n))
        {
            *info = -8;
        }
    }
    if (*info == 0)
    {
        /* Determine the optimum and minimum LWORK */
        if (*m >= *n)
        {
            zgeqr_(m, n, &a[a_offset], lda, tq, &c_n1, workq, &c_n1, &info2);
            tszo = (integer) tq[0].r;
            lwo = (integer) workq[0].r;
            zgemqr_("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszo, &b[ b_offset], ldb, workq, &c_n1, &info2);
            /* Computing MAX */
            i__1 = lwo;
            i__2 = (integer) workq[0].r; // , expr subst
            lwo = fla_max(i__1,i__2);
            zgeqr_(m, n, &a[a_offset], lda, tq, &c_n2, workq, &c_n2, &info2);
            tszm = (integer) tq[0].r;
            lwm = (integer) workq[0].r;
            zgemqr_("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszm, &b[ b_offset], ldb, workq, &c_n1, &info2);
            /* Computing MAX */
            i__1 = lwm;
            i__2 = (integer) workq[0].r; // , expr subst
            lwm = fla_max(i__1,i__2);
            wsizeo = tszo + lwo;
            wsizem = tszm + lwm;
        }
        else
        {
            zgelq_(m, n, &a[a_offset], lda, tq, &c_n1, workq, &c_n1, &info2);
            tszo = (integer) tq[0].r;
            lwo = (integer) workq[0].r;
            zgemlq_("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszo, &b[ b_offset], ldb, workq, &c_n1, &info2);
            /* Computing MAX */
            i__1 = lwo;
            i__2 = (integer) workq[0].r; // , expr subst
            lwo = fla_max(i__1,i__2);
            zgelq_(m, n, &a[a_offset], lda, tq, &c_n2, workq, &c_n2, &info2);
            tszm = (integer) tq[0].r;
            lwm = (integer) workq[0].r;
            zgemlq_("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszm, &b[ b_offset], ldb, workq, &c_n1, &info2);
            /* Computing MAX */
            i__1 = lwm;
            i__2 = (integer) workq[0].r; // , expr subst
            lwm = fla_max(i__1,i__2);
            wsizeo = tszo + lwo;
            wsizem = tszm + lwm;
        }
        if (*lwork < wsizem && ! lquery)
        {
            *info = -10;
        }
        d__1 = (doublereal) wsizeo;
        work[1].r = d__1;
        work[1].i = 0.; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGETSLS", &i__1);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (lquery)
    {
        if (*lwork == -2)
        {
            d__1 = (doublereal) wsizem;
            work[1].r = d__1;
            work[1].i = 0.; // , expr subst
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    if (*lwork < wsizeo)
    {
        lw1 = tszm;
        lw2 = lwm;
    }
    else
    {
        lw1 = tszo;
        lw2 = lwo;
    }
    /* Quick return if possible */
    /* Computing MIN */
    i__1 = fla_min(*m,*n);
    if (fla_min(i__1,*nrhs) == 0)
    {
        i__1 = fla_max(*m,*n);
        zlaset_("FULL", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
    /* Get machine parameters */
    smlnum = dlamch_("S") / dlamch_("P");
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);
    /* Scale A, B if max element outside range [SMLNUM,BIGNUM] */
    anrm = zlange_("M", m, n, &a[a_offset], lda, dum);
    iascl = 0;
    if (anrm > 0. && anrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, info);
        iascl = 1;
    }
    else if (anrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, info);
        iascl = 2;
    }
    else if (anrm == 0.)
    {
        /* Matrix all zero. Return zero solution. */
        zlaset_("F", &maxmn, nrhs, &c_b1, &c_b1, &b[b_offset], ldb) ;
        goto L50;
    }
    brow = *m;
    if (tran)
    {
        brow = *n;
    }
    bnrm = zlange_("M", &brow, nrhs, &b[b_offset], ldb, dum);
    ibscl = 0;
    if (bnrm > 0. && bnrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        zlascl_("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], ldb, info);
        ibscl = 1;
    }
    else if (bnrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        zlascl_("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], ldb, info);
        ibscl = 2;
    }
    if (*m >= *n)
    {
        /* compute QR factorization of A */
        zgeqr_(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, info);
        if (! tran)
        {
            /* Least-Squares Problem min || A * X - B || */
            /* B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */
            zgemqr_("L", "C", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], & lw1, &b[b_offset], ldb, &work[1], &lw2, info);
            /* B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */
            ztrtrs_("U", "N", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
            scllen = *n;
        }
        else
        {
            /* Overdetermined system of equations A**T * X = B */
            /* B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS) */
            ztrtrs_("U", "C", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return 0;
            }
            /* B(N+1:M,1:NRHS) = CZERO */
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
                    b[i__3].r = 0.;
                    b[i__3].i = 0.; // , expr subst
                    /* L10: */
                }
                /* L20: */
            }
            /* B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */
            zgemqr_("L", "N", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], & lw1, &b[b_offset], ldb, &work[1], &lw2, info);
            scllen = *m;
        }
    }
    else
    {
        /* Compute LQ factorization of A */
        zgelq_(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, info);
        /* workspace at least M, optimally M*NB. */
        if (! tran)
        {
            /* underdetermined system of equations A * X = B */
            /* B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */
            ztrtrs_("L", "N", "N", m, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info);
            if (*info > 0)
            {
                AOCL_DTL_TRACE_LOG_EXIT
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
                    b[i__3].r = 0.;
                    b[i__3].i = 0.; // , expr subst
                    /* L30: */
                }
                /* L40: */
            }
            /* B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS) */
            zgemlq_("L", "C", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], & lw1, &b[b_offset], ldb, &work[1], &lw2, info);
            /* workspace at least NRHS, optimally NRHS*NB */
            scllen = *n;
        }
        else
        {
            /* overdetermined system min || A**T * X - B || */
            /* B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */
            zgemlq_("L", "N", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], & lw1, &b[b_offset], ldb, &work[1], &lw2, info);
            /* workspace at least NRHS, optimally NRHS*NB */
            /* B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS) */
            ztrtrs_("L", "C", "N", m, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info);
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
        zlascl_("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset], ldb, info);
    }
    else if (iascl == 2)
    {
        zlascl_("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset], ldb, info);
    }
    if (ibscl == 1)
    {
        zlascl_("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset], ldb, info);
    }
    else if (ibscl == 2)
    {
        zlascl_("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset], ldb, info);
    }
L50:
    d__1 = (doublereal) (tszo + lwo);
    work[1].r = d__1;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of ZGETSLS */
}
/* zgetsls_ */
