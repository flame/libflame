/* ../netlib/sggrqf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b SGGRQF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGGRQF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggrqf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggrqf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggrqf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGGRQF( M, P, N, A, LDA, TAUA, B, LDB, TAUB, WORK, */
/* LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDB, LWORK, M, N, P */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGRQF computes a generalized RQ factorization of an M-by-N matrix A */
/* > and a P-by-N matrix B: */
/* > */
/* > A = R*Q, B = Z*T*Q, */
/* > */
/* > where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal */
/* > matrix, and R and T assume one of the forms: */
/* > */
/* > if M <= N, R = ( 0 R12 ) M, or if M > N, R = ( R11 ) M-N, */
/* > N-M M ( R21 ) N */
/* > N */
/* > */
/* > where R12 or R21 is upper triangular, and */
/* > */
/* > if P >= N, T = ( T11 ) N , or if P < N, T = ( T11 T12 ) P, */
/* > ( 0 ) P-N P N-P */
/* > N */
/* > */
/* > where T11 is upper triangular. */
/* > */
/* > In particular, if B is square and nonsingular, the GRQ factorization */
/* > of A and B implicitly gives the RQ factorization of A*inv(B): */
/* > */
/* > A*inv(B) = (R*inv(T))*Z**T */
/* > */
/* > where inv(B) denotes the inverse of the matrix B, and Z**T denotes the */
/* > transpose of the matrix Z. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* > P is INTEGER */
/* > The number of rows of the matrix B. P >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, if M <= N, the upper triangle of the subarray */
/* > A(1:M,N-M+1:N) contains the M-by-M upper triangular matrix R;
*/
/* > if M > N, the elements on and above the (M-N)-th subdiagonal */
/* > contain the M-by-N upper trapezoidal matrix R;
the remaining */
/* > elements, with the array TAUA, represent the orthogonal */
/* > matrix Q as a product of elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAUA */
/* > \verbatim */
/* > TAUA is REAL array, dimension (min(M,N)) */
/* > The scalar factors of the elementary reflectors which */
/* > represent the orthogonal matrix Q (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,N) */
/* > On entry, the P-by-N matrix B. */
/* > On exit, the elements on and above the diagonal of the array */
/* > contain the min(P,N)-by-N upper trapezoidal matrix T (T is */
/* > upper triangular if P >= N);
the elements below the diagonal, */
/* > with the array TAUB, represent the orthogonal matrix Z as a */
/* > product of elementary reflectors (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,P). */
/* > \endverbatim */
/* > */
/* > \param[out] TAUB */
/* > \verbatim */
/* > TAUB is REAL array, dimension (min(P,N)) */
/* > The scalar factors of the elementary reflectors which */
/* > represent the orthogonal matrix Z (see Further Details). */
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
/* > The dimension of the array WORK. LWORK >= max(1,N,M,P). */
/* > For optimum performance LWORK >= max(N,M,P)*max(NB1,NB2,NB3), */
/* > where NB1 is the optimal blocksize for the RQ factorization */
/* > of an M-by-N matrix, NB2 is the optimal blocksize for the */
/* > QR factorization of a P-by-N matrix, and NB3 is the optimal */
/* > blocksize for a call of SORMRQ. */
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
/* > < 0: if INF0= -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realOTHERcomputational */
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
/* > H(i) = I - taua * v * v**T */
/* > */
/* > where taua is a real scalar, and v is a real vector with */
/* > v(n-k+i+1:n) = 0 and v(n-k+i) = 1;
v(1:n-k+i-1) is stored on exit in */
/* > A(m-k+i,1:n-k+i-1), and taua in TAUA(i). */
/* > To form Q explicitly, use LAPACK subroutine SORGRQ. */
/* > To use Q to update another matrix, use LAPACK subroutine SORMRQ. */
/* > */
/* > The matrix Z is represented as a product of elementary reflectors */
/* > */
/* > Z = H(1) H(2) . . . H(k), where k = min(p,n). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - taub * v * v**T */
/* > */
/* > where taub is a real scalar, and v is a real vector with */
/* > v(1:i-1) = 0 and v(i) = 1;
v(i+1:p) is stored on exit in B(i+1:p,i), */
/* > and taub in TAUB(i). */
/* > To form Z explicitly, use LAPACK subroutine SORGQR. */
/* > To use Z to update another matrix, use LAPACK subroutine SORMQR. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int sggrqf_(integer *m, integer *p, integer *n, real *a, integer *lda, real *taua, real *b, integer *ldb, real *taub, real * work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    /* Local variables */
    integer nb, nb1, nb2, nb3, lopt;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int sgeqrf_(integer *, integer *, real *, integer *, real *, real *, integer *, integer *), sgerqf_(integer *, integer *, real *, integer *, real *, real *, integer *, integer * );
    integer lwkopt;
    logical lquery;
    extern /* Subroutine */
    int sormrq_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --taua;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --taub;
    --work;
    /* Function Body */
    *info = 0;
    nb1 = ilaenv_(&c__1, "SGERQF", " ", m, n, &c_n1, &c_n1);
    nb2 = ilaenv_(&c__1, "SGEQRF", " ", p, n, &c_n1, &c_n1);
    nb3 = ilaenv_(&c__1, "SORMRQ", " ", m, n, p, &c_n1);
    /* Computing MAX */
    i__1 = max(nb1,nb2);
    nb = max(i__1,nb3);
    /* Computing MAX */
    i__1 = max(*n,*m);
    lwkopt = max(i__1,*p) * nb;
    work[1] = (real) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*p < 0)
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*lda < max(1,*m))
    {
        *info = -5;
    }
    else if (*ldb < max(1,*p))
    {
        *info = -8;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = max(1,*m);
        i__1 = max(i__1,*p); // , expr subst
        if (*lwork < max(i__1,*n) && ! lquery)
        {
            *info = -11;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGGRQF", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* RQ factorization of M-by-N matrix A: A = R*Q */
    sgerqf_(m, n, &a[a_offset], lda, &taua[1], &work[1], lwork, info);
    lopt = work[1];
    /* Update B := B*Q**T */
    i__1 = min(*m,*n);
    /* Computing MAX */
    i__2 = 1;
    i__3 = *m - *n + 1; // , expr subst
    sormrq_("Right", "Transpose", p, n, &i__1, &a[max(i__2,i__3) + a_dim1], lda, &taua[1], &b[b_offset], ldb, &work[1], lwork, info);
    /* Computing MAX */
    i__1 = lopt;
    i__2 = (integer) work[1]; // , expr subst
    lopt = max(i__1,i__2);
    /* QR factorization of P-by-N matrix B: B = Z*T */
    sgeqrf_(p, n, &b[b_offset], ldb, &taub[1], &work[1], lwork, info);
    /* Computing MAX */
    i__1 = lopt;
    i__2 = (integer) work[1]; // , expr subst
    work[1] = (real) max(i__1,i__2);
    return 0;
    /* End of SGGRQF */
}
/* sggrqf_ */
