/* ../netlib/zgerqf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
/* > \brief \b ZGERQF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGERQF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgerqf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgerqf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgerqf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGERQF computes an RQ factorization of a complex M-by-N matrix A: */
/* > A = R * Q. */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, */
/* > if m <= n, the upper triangle of the subarray */
/* > A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R;
*/
/* > if m >= n, the elements on and above the (m-n)-th subdiagonal */
/* > contain the M-by-N upper trapezoidal matrix R;
*/
/* > the remaining elements, with the array TAU, represent the */
/* > unitary matrix Q as a product of min(m,n) elementary */
/* > reflectors (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 array, dimension (min(M,N)) */
/* > The scalar factors of the elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= max(1,M). */
/* > For optimum performance LWORK >= M*NB, where NB is */
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
/* > \ingroup complex16GEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of elementary reflectors */
/* > */
/* > Q = H(1)**H H(2)**H . . . H(k)**H, where k = min(m,n). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar, and v is a complex vector with */
/* > v(n-k+i+1:n) = 0 and v(n-k+i) = 1;
conjg(v(1:n-k+i-1)) is stored on */
/* > exit in A(m-k+i,1:n-k+i-1), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int zgerqf_(integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    integer i__, k, ib, nb, ki, kk, mu, nu, nx, iws, nbmin, iinfo;
    extern /* Subroutine */
    int zgerq2_(integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *), xerbla_( char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int zlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    integer ldwork;
    extern /* Subroutine */
    int zlarft_(char *, char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
    integer lwkopt;
    logical lquery;
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
        k = min(*m,*n);
        if (k == 0)
        {
            lwkopt = 1;
        }
        else
        {
            nb = ilaenv_(&c__1, "ZGERQF", " ", m, n, &c_n1, &c_n1);
            lwkopt = *m * nb;
        }
        work[1].r = (doublereal) lwkopt;
        work[1].i = 0.; // , expr subst
        if (*lwork < max(1,*m) && ! lquery)
        {
            *info = -7;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGERQF", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (k == 0)
    {
        return 0;
    }
    nbmin = 2;
    nx = 1;
    iws = *m;
    if (nb > 1 && nb < k)
    {
        /* Determine when to cross over from blocked to unblocked code. */
        /* Computing MAX */
        i__1 = 0;
        i__2 = ilaenv_(&c__3, "ZGERQF", " ", m, n, &c_n1, &c_n1); // , expr subst
        nx = max(i__1,i__2);
        if (nx < k)
        {
            /* Determine if workspace is large enough for blocked code. */
            ldwork = *m;
            iws = ldwork * nb;
            if (*lwork < iws)
            {
                /* Not enough workspace to use optimal NB: reduce NB and */
                /* determine the minimum value of NB. */
                nb = *lwork / ldwork;
                /* Computing MAX */
                i__1 = 2;
                i__2 = ilaenv_(&c__2, "ZGERQF", " ", m, n, &c_n1, & c_n1); // , expr subst
                nbmin = max(i__1,i__2);
            }
        }
    }
    if (nb >= nbmin && nb < k && nx < k)
    {
        /* Use blocked code initially. */
        /* The last kk rows are handled by the block method. */
        ki = (k - nx - 1) / nb * nb;
        /* Computing MIN */
        i__1 = k;
        i__2 = ki + nb; // , expr subst
        kk = min(i__1,i__2);
        i__1 = k - kk + 1;
        i__2 = -nb;
        for (i__ = k - kk + ki + 1;
                i__2 < 0 ? i__ >= i__1 : i__ <= i__1;
                i__ += i__2)
        {
            /* Computing MIN */
            i__3 = k - i__ + 1;
            ib = min(i__3,nb);
            /* Compute the RQ factorization of the current block */
            /* A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1) */
            i__3 = *n - k + i__ + ib - 1;
            zgerq2_(&ib, &i__3, &a[*m - k + i__ + a_dim1], lda, &tau[i__], & work[1], &iinfo);
            if (*m - k + i__ > 1)
            {
                /* Form the triangular factor of the block reflector */
                /* H = H(i+ib-1) . . . H(i+1) H(i) */
                i__3 = *n - k + i__ + ib - 1;
                zlarft_("Backward", "Rowwise", &i__3, &ib, &a[*m - k + i__ + a_dim1], lda, &tau[i__], &work[1], &ldwork);
                /* Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right */
                i__3 = *m - k + i__ - 1;
                i__4 = *n - k + i__ + ib - 1;
                zlarfb_("Right", "No transpose", "Backward", "Rowwise", &i__3, &i__4, &ib, &a[*m - k + i__ + a_dim1], lda, &work[1], &ldwork, &a[a_offset], lda, &work[ib + 1], &ldwork);
            }
            /* L10: */
        }
        mu = *m - k + i__ + nb - 1;
        nu = *n - k + i__ + nb - 1;
    }
    else
    {
        mu = *m;
        nu = *n;
    }
    /* Use unblocked code to factor the last or only block */
    if (mu > 0 && nu > 0)
    {
        zgerq2_(&mu, &nu, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
    }
    work[1].r = (doublereal) iws;
    work[1].i = 0.; // , expr subst
    return 0;
    /* End of ZGERQF */
}
/* zgerqf_ */
