/* ../netlib/cgeqp3.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
/* > \brief \b CGEQP3 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGEQP3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqp3. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqp3. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqp3. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPVT( * ) */
/* REAL RWORK( * ) */
/* COMPLEX A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEQP3 computes a QR factorization with column pivoting of a */
/* > matrix A: A*P = Q*R using Level 3 BLAS. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the upper triangle of the array contains the */
/* > min(M,N)-by-N upper trapezoidal matrix R;
the elements below */
/* > the diagonal, together with the array TAU, represent the */
/* > unitary matrix Q as a product of min(M,N) elementary */
/* > reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] JPVT */
/* > \verbatim */
/* > JPVT is INTEGER array, dimension (N) */
/* > On entry, if JPVT(J).ne.0, the J-th column of A is permuted */
/* > to the front of A*P (a leading column);
if JPVT(J)=0, */
/* > the J-th column of A is a free column. */
/* > On exit, if JPVT(J)=K, then the J-th column of A*P was the */
/* > the K-th column of A. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (min(M,N)) */
/* > The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO=0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= N+1. */
/* > For optimal performance LWORK >= ( N+1 )*NB, where NB */
/* > is the optimal blocksize. */
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
/* > RWORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(k), where k = min(m,n). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar, and v is a real/complex vector */
/* > with v(1:i-1) = 0 and v(i) = 1;
v(i+1:m) is stored on exit in */
/* > A(i+1:m,i), and tau in TAU(i). */
/* > \endverbatim */
/* > \par Contributors: */
/* ================== */
/* > */
/* > G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain */
/* > X. Sun, Computer Science Dept., Duke University, USA */
/* > */
/* ===================================================================== */
/* Subroutine */
int cgeqp3_(integer *m, integer *n, complex *a, integer *lda, integer *jpvt, complex *tau, complex *work, integer *lwork, real * rwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    integer j, jb, na, nb, sm, sn, nx, fjb, iws, nfxd, nbmin;
    extern /* Subroutine */
    int cswap_(integer *, complex *, integer *, complex *, integer *);
    integer minmn, minws;
    extern /* Subroutine */
    int claqp2_(integer *, integer *, integer *, complex *, integer *, integer *, complex *, real *, real *, complex *);
    extern real scnrm2_(integer *, complex *, integer *);
    extern /* Subroutine */
    int cgeqrf_(integer *, integer *, complex *, integer *, complex *, complex *, integer *, integer *), xerbla_( char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int claqps_(integer *, integer *, integer *, integer *, integer *, complex *, integer *, integer *, complex *, real *, real *, complex *, complex *, integer *);
    integer topbmn, sminmn;
    extern /* Subroutine */
    int cunmqr_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, complex *, integer *, complex *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* Test input arguments */
    /* ==================== */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --tau;
    --work;
    --rwork;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < max(1,*m))
    {
        *info = -4;
    }
    if (*info == 0)
    {
        minmn = min(*m,*n);
        if (minmn == 0)
        {
            iws = 1;
            lwkopt = 1;
        }
        else
        {
            iws = *n + 1;
            nb = ilaenv_(&c__1, "CGEQRF", " ", m, n, &c_n1, &c_n1);
            lwkopt = (*n + 1) * nb;
        }
        work[1].r = (real) lwkopt;
        work[1].i = 0.f; // , expr subst
        if (*lwork < iws && ! lquery)
        {
            *info = -8;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGEQP3", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible. */
    if (minmn == 0)
    {
        return 0;
    }
    /* Move initial columns up front. */
    nfxd = 1;
    i__1 = *n;
    for (j = 1;
            j <= i__1;
            ++j)
    {
        if (jpvt[j] != 0)
        {
            if (j != nfxd)
            {
                cswap_(m, &a[j * a_dim1 + 1], &c__1, &a[nfxd * a_dim1 + 1], & c__1);
                jpvt[j] = jpvt[nfxd];
                jpvt[nfxd] = j;
            }
            else
            {
                jpvt[j] = j;
            }
            ++nfxd;
        }
        else
        {
            jpvt[j] = j;
        }
        /* L10: */
    }
    --nfxd;
    /* Factorize fixed columns */
    /* ======================= */
    /* Compute the QR factorization of fixed columns and update */
    /* remaining columns. */
    if (nfxd > 0)
    {
        na = min(*m,nfxd);
        /* CC CALL CGEQR2( M, NA, A, LDA, TAU, WORK, INFO ) */
        cgeqrf_(m, &na, &a[a_offset], lda, &tau[1], &work[1], lwork, info);
        /* Computing MAX */
        i__1 = iws;
        i__2 = (integer) work[1].r; // , expr subst
        iws = max(i__1,i__2);
        if (na < *n)
        {
            /* CC CALL CUNM2R( 'Left', 'Conjugate Transpose', M, N-NA, */
            /* CC $ NA, A, LDA, TAU, A( 1, NA+1 ), LDA, WORK, */
            /* CC $ INFO ) */
            i__1 = *n - na;
            cunmqr_("Left", "Conjugate Transpose", m, &i__1, &na, &a[a_offset] , lda, &tau[1], &a[(na + 1) * a_dim1 + 1], lda, &work[1], lwork, info);
            /* Computing MAX */
            i__1 = iws;
            i__2 = (integer) work[1].r; // , expr subst
            iws = max(i__1,i__2);
        }
    }
    /* Factorize free columns */
    /* ====================== */
    if (nfxd < minmn)
    {
        sm = *m - nfxd;
        sn = *n - nfxd;
        sminmn = minmn - nfxd;
        /* Determine the block size. */
        nb = ilaenv_(&c__1, "CGEQRF", " ", &sm, &sn, &c_n1, &c_n1);
        nbmin = 2;
        nx = 0;
        if (nb > 1 && nb < sminmn)
        {
            /* Determine when to cross over from blocked to unblocked code. */
            /* Computing MAX */
            i__1 = 0;
            i__2 = ilaenv_(&c__3, "CGEQRF", " ", &sm, &sn, &c_n1, & c_n1); // , expr subst
            nx = max(i__1,i__2);
            if (nx < sminmn)
            {
                /* Determine if workspace is large enough for blocked code. */
                minws = (sn + 1) * nb;
                iws = max(iws,minws);
                if (*lwork < minws)
                {
                    /* Not enough workspace to use optimal NB: Reduce NB and */
                    /* determine the minimum value of NB. */
                    nb = *lwork / (sn + 1);
                    /* Computing MAX */
                    i__1 = 2;
                    i__2 = ilaenv_(&c__2, "CGEQRF", " ", &sm, &sn, & c_n1, &c_n1); // , expr subst
                    nbmin = max(i__1,i__2);
                }
            }
        }
        /* Initialize partial column norms. The first N elements of work */
        /* store the exact column norms. */
        i__1 = *n;
        for (j = nfxd + 1;
                j <= i__1;
                ++j)
        {
            rwork[j] = scnrm2_(&sm, &a[nfxd + 1 + j * a_dim1], &c__1);
            rwork[*n + j] = rwork[j];
            /* L20: */
        }
        if (nb >= nbmin && nb < sminmn && nx < sminmn)
        {
            /* Use blocked code initially. */
            j = nfxd + 1;
            /* Compute factorization: while loop. */
            topbmn = minmn - nx;
L30:
            if (j <= topbmn)
            {
                /* Computing MIN */
                i__1 = nb;
                i__2 = topbmn - j + 1; // , expr subst
                jb = min(i__1,i__2);
                /* Factorize JB columns among columns J:N. */
                i__1 = *n - j + 1;
                i__2 = j - 1;
                i__3 = *n - j + 1;
                claqps_(m, &i__1, &i__2, &jb, &fjb, &a[j * a_dim1 + 1], lda, & jpvt[j], &tau[j], &rwork[j], &rwork[*n + j], &work[1], &work[jb + 1], &i__3);
                j += fjb;
                goto L30;
            }
        }
        else
        {
            j = nfxd + 1;
        }
        /* Use unblocked code to factor the last or only block. */
        if (j <= minmn)
        {
            i__1 = *n - j + 1;
            i__2 = j - 1;
            claqp2_(m, &i__1, &i__2, &a[j * a_dim1 + 1], lda, &jpvt[j], &tau[ j], &rwork[j], &rwork[*n + j], &work[1]);
        }
    }
    work[1].r = (real) iws;
    work[1].i = 0.f; // , expr subst
    return 0;
    /* End of CGEQP3 */
}
/* cgeqp3_ */
