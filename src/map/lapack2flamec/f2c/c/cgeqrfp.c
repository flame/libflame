/* ../netlib/cgeqrfp.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
/* > \brief \b CGEQRFP */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGEQRFP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqrfp .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqrfp .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqrfp .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEQRFP computes a QR factorization of a complex M-by-N matrix A: */
/* > A = Q * R. */
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
/* > On exit, the elements on and above the diagonal of the array */
/* > contain the fla_min(M,N)-by-N upper trapezoidal matrix R (R is */
/* > upper triangular if m >= n);
the elements below the diagonal, */
/* > with the array TAU, represent the unitary matrix Q as a */
/* > product of fla_min(m,n) elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors (see Further */
/* > Details). */
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
/* > The dimension of the array WORK. LWORK >= fla_max(1,N). */
/* > For optimum performance LWORK >= N*NB, where NB is */
/* > the optimal blocksize. */
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
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(k), where k = fla_min(m,n). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar, and v is a complex vector with */
/* > v(1:i-1) = 0 and v(i) = 1;
v(i+1:m) is stored on exit in A(i+1:m,i), */
/* > and tau in TAU(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int cgeqrfp_(integer *m, integer *n, complex *a, integer * lda, complex *tau, complex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"cgeqrfp inputs: m %lld, n %lld, lda %lld, lwork %lld",*m, *n, *lda, *lwork);
#else
    snprintf(buffer, 256,"cgeqrfp inputs: m %d, n %d, lda %d, lwork %d",*m, *n, *lda, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern /* Subroutine */
    int clarfb_(char *, char *, char *, char *, integer *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, integer *, complex *, integer *), clarft_(char *, char *, integer *, integer *, complex *, integer *, complex *, complex *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer ldwork, lwkopt;
    logical lquery;
    extern /* Subroutine */
    int cgeqr2p_(integer *, integer *, complex *, integer *, complex *, complex *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "CGEQRF", " ", m, n, &c_n1, &c_n1);
    lwkopt = *n * nb;
    work[1].r = (real) lwkopt;
    work[1].i = 0.f; // , expr subst
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -4;
    }
    else if (*lwork < fla_max(1,*n) && ! lquery)
    {
        *info = -7;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGEQRFP", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    k = fla_min(*m,*n);
    if (k == 0)
    {
        work[1].r = 1.f;
        work[1].i = 0.f; // , expr subst
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k)
    {
        /* Determine when to cross over from blocked to unblocked code. */
        /* Computing MAX */
        i__1 = 0;
        i__2 = ilaenv_(&c__3, "CGEQRF", " ", m, n, &c_n1, &c_n1); // , expr subst
        nx = fla_max(i__1,i__2);
        if (nx < k)
        {
            /* Determine if workspace is large enough for blocked code. */
            ldwork = *n;
            iws = ldwork * nb;
            if (*lwork < iws)
            {
                /* Not enough workspace to use optimal NB: reduce NB and */
                /* determine the minimum value of NB. */
                nb = *lwork / ldwork;
                /* Computing MAX */
                i__1 = 2;
                i__2 = ilaenv_(&c__2, "CGEQRF", " ", m, n, &c_n1, & c_n1); // , expr subst
                nbmin = fla_max(i__1,i__2);
            }
        }
    }
    if (nb >= nbmin && nb < k && nx < k)
    {
        /* Use blocked code initially */
        i__1 = k - nx;
        i__2 = nb;
        for (i__ = 1;
                i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
                i__ += i__2)
        {
            /* Computing MIN */
            i__3 = k - i__ + 1;
            ib = fla_min(i__3,nb);
            /* Compute the QR factorization of the current block */
            /* A(i:m,i:i+ib-1) */
            i__3 = *m - i__ + 1;
            cgeqr2p_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], & work[1], &iinfo);
            if (i__ + ib <= *n)
            {
                /* Form the triangular factor of the block reflector */
                /* H = H(i) H(i+1) . . . H(i+ib-1) */
                i__3 = *m - i__ + 1;
                clarft_("Forward", "Columnwise", &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &ldwork);
                /* Apply H**H to A(i:m,i+ib:n) from the left */
                i__3 = *m - i__ + 1;
                i__4 = *n - i__ - ib + 1;
                clarfb_("Left", "Conjugate transpose", "Forward", "Columnwise", &i__3, &i__4, &ib, &a[i__ + i__ * a_dim1], lda, & work[1], &ldwork, &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib + 1], &ldwork);
            }
            /* L10: */
        }
    }
    else
    {
        i__ = 1;
    }
    /* Use unblocked code to factor the last or only block. */
    if (i__ <= k)
    {
        i__2 = *m - i__ + 1;
        i__1 = *n - i__ + 1;
        cgeqr2p_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[ 1], &iinfo);
    }
    work[1].r = (real) iws;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CGEQRFP */
}
/* cgeqrfp_ */
