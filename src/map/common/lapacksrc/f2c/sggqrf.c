/* ../netlib/sggqrf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b SGGQRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGGQRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggqrf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggqrf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggqrf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGGQRF( N, M, P, A, LDA, TAUA, B, LDB, TAUB, WORK, */
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
/* > SGGQRF computes a generalized QR factorization of an N-by-M matrix A */
/* > and an N-by-P matrix B: */
/* > */
/* > A = Q*R, B = Q*T*Z, */
/* > */
/* > where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal */
/* > matrix, and R and T assume one of the forms: */
/* > */
/* > if N >= M, R = ( R11 ) M , or if N < M, R = ( R11 R12 ) N, */
/* > ( 0 ) N-M N M-N */
/* > M */
/* > */
/* > where R11 is upper triangular, and */
/* > */
/* > if N <= P, T = ( 0 T12 ) N, or if N > P, T = ( T11 ) N-P, */
/* > P-N N ( T21 ) P */
/* > P */
/* > */
/* > where T12 or T21 is upper triangular. */
/* > */
/* > In particular, if B is square and nonsingular, the GQR factorization */
/* > of A and B implicitly gives the QR factorization of inv(B)*A: */
/* > */
/* > inv(B)*A = Z**T*(inv(T)*R) */
/* > */
/* > where inv(B) denotes the inverse of the matrix B, and Z**T denotes the */
/* > transpose of the matrix Z. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of rows of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of columns of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* > P is INTEGER */
/* > The number of columns of the matrix B. P >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,M) */
/* > On entry, the N-by-M matrix A. */
/* > On exit, the elements on and above the diagonal of the array */
/* > contain the min(N,M)-by-M upper trapezoidal matrix R (R is */
/* > upper triangular if N >= M);
the elements below the diagonal, */
/* > with the array TAUA, represent the orthogonal matrix Q as a */
/* > product of min(N,M) elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] TAUA */
/* > \verbatim */
/* > TAUA is REAL array, dimension (min(N,M)) */
/* > The scalar factors of the elementary reflectors which */
/* > represent the orthogonal matrix Q (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,P) */
/* > On entry, the N-by-P matrix B. */
/* > On exit, if N <= P, the upper triangle of the subarray */
/* > B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T;
*/
/* > if N > P, the elements on and above the (N-P)-th subdiagonal */
/* > contain the N-by-P upper trapezoidal matrix T;
the remaining */
/* > elements, with the array TAUB, represent the orthogonal */
/* > matrix Z as a product of elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] TAUB */
/* > \verbatim */
/* > TAUB is REAL array, dimension (min(N,P)) */
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
/* > where NB1 is the optimal blocksize for the QR factorization */
/* > of an N-by-M matrix, NB2 is the optimal blocksize for the */
/* > RQ factorization of an N-by-P matrix, and NB3 is the optimal */
/* > blocksize for a call of SORMQR. */
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
/* > Q = H(1) H(2) . . . H(k), where k = min(n,m). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - taua * v * v**T */
/* > */
/* > where taua is a real scalar, and v is a real vector with */
/* > v(1:i-1) = 0 and v(i) = 1;
v(i+1:n) is stored on exit in A(i+1:n,i), */
/* > and taua in TAUA(i). */
/* > To form Q explicitly, use LAPACK subroutine SORGQR. */
/* > To use Q to update another matrix, use LAPACK subroutine SORMQR. */
/* > */
/* > The matrix Z is represented as a product of elementary reflectors */
/* > */
/* > Z = H(1) H(2) . . . H(k), where k = min(n,p). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - taub * v * v**T */
/* > */
/* > where taub is a real scalar, and v is a real vector with */
/* > v(p-k+i+1:p) = 0 and v(p-k+i) = 1;
v(1:p-k+i-1) is stored on exit in */
/* > B(n-k+i,1:p-k+i-1), and taub in TAUB(i). */
/* > To form Z explicitly, use LAPACK subroutine SORGRQ. */
/* > To use Z to update another matrix, use LAPACK subroutine SORMRQ. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int sggqrf_(integer *n, integer *m, integer *p, real *a, integer *lda, real *taua, real *b, integer *ldb, real *taub, real * work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
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
    int sormqr_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *, integer *, real *, integer *, integer *);
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
    nb1 = ilaenv_(&c__1, "SGEQRF", " ", n, m, &c_n1, &c_n1);
    nb2 = ilaenv_(&c__1, "SGERQF", " ", n, p, &c_n1, &c_n1);
    nb3 = ilaenv_(&c__1, "SORMQR", " ", n, m, p, &c_n1);
    /* Computing MAX */
    i__1 = max(nb1,nb2);
    nb = max(i__1,nb3);
    /* Computing MAX */
    i__1 = max(*n,*m);
    lwkopt = max(i__1,*p) * nb;
    work[1] = (real) lwkopt;
    lquery = *lwork == -1;
    if (*n < 0)
    {
        *info = -1;
    }
    else if (*m < 0)
    {
        *info = -2;
    }
    else if (*p < 0)
    {
        *info = -3;
    }
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -8;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = max(1,*n);
        i__1 = max(i__1,*m); // , expr subst
        if (*lwork < max(i__1,*p) && ! lquery)
        {
            *info = -11;
        }
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGGQRF", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* QR factorization of N-by-M matrix A: A = Q*R */
    sgeqrf_(n, m, &a[a_offset], lda, &taua[1], &work[1], lwork, info);
    lopt = work[1];
    /* Update B := Q**T*B. */
    i__1 = min(*n,*m);
    sormqr_("Left", "Transpose", n, p, &i__1, &a[a_offset], lda, &taua[1], &b[ b_offset], ldb, &work[1], lwork, info);
    /* Computing MAX */
    i__1 = lopt;
    i__2 = (integer) work[1]; // , expr subst
    lopt = max(i__1,i__2);
    /* RQ factorization of N-by-P matrix B: B = T*Z. */
    sgerqf_(n, p, &b[b_offset], ldb, &taub[1], &work[1], lwork, info);
    /* Computing MAX */
    i__1 = lopt;
    i__2 = (integer) work[1]; // , expr subst
    work[1] = (real) max(i__1,i__2);
    return 0;
    /* End of SGGQRF */
}
/* sggqrf_ */
