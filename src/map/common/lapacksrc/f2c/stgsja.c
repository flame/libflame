/* ../netlib/stgsja.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b13 = 0.f;
static real c_b14 = 1.f;
static integer c__1 = 1;
static real c_b43 = -1.f;
/* > \brief \b STGSJA */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download STGSJA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsja. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsja. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsja. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE STGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, */
/* LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, */
/* Q, LDQ, WORK, NCYCLE, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBQ, JOBU, JOBV */
/* INTEGER INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, */
/* $ NCYCLE, P */
/* REAL TOLA, TOLB */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/* $ BETA( * ), Q( LDQ, * ), U( LDU, * ), */
/* $ V( LDV, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > STGSJA computes the generalized singular value decomposition (GSVD) */
/* > of two real upper triangular (or trapezoidal) matrices A and B. */
/* > */
/* > On entry, it is assumed that matrices A and B have the following */
/* > forms, which may be obtained by the preprocessing subroutine SGGSVP */
/* > from a general M-by-N matrix A and P-by-N matrix B: */
/* > */
/* > N-K-L K L */
/* > A = K ( 0 A12 A13 ) if M-K-L >= 0;
*/
/* > L ( 0 0 A23 ) */
/* > M-K-L ( 0 0 0 ) */
/* > */
/* > N-K-L K L */
/* > A = K ( 0 A12 A13 ) if M-K-L < 0;
*/
/* > M-K ( 0 0 A23 ) */
/* > */
/* > N-K-L K L */
/* > B = L ( 0 0 B13 ) */
/* > P-L ( 0 0 0 ) */
/* > */
/* > where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular */
/* > upper triangular;
A23 is L-by-L upper triangular if M-K-L >= 0, */
/* > otherwise A23 is (M-K)-by-L upper trapezoidal. */
/* > */
/* > On exit, */
/* > */
/* > U**T *A*Q = D1*( 0 R ), V**T *B*Q = D2*( 0 R ), */
/* > */
/* > where U, V and Q are orthogonal matrices. */
/* > R is a nonsingular upper triangular matrix, and D1 and D2 are */
/* > ``diagonal'' matrices, which are of the following structures: */
/* > */
/* > If M-K-L >= 0, */
/* > */
/* > K L */
/* > D1 = K ( I 0 ) */
/* > L ( 0 C ) */
/* > M-K-L ( 0 0 ) */
/* > */
/* > K L */
/* > D2 = L ( 0 S ) */
/* > P-L ( 0 0 ) */
/* > */
/* > N-K-L K L */
/* > ( 0 R ) = K ( 0 R11 R12 ) K */
/* > L ( 0 0 R22 ) L */
/* > */
/* > where */
/* > */
/* > C = diag( ALPHA(K+1), ... , ALPHA(K+L) ), */
/* > S = diag( BETA(K+1), ... , BETA(K+L) ), */
/* > C**2 + S**2 = I. */
/* > */
/* > R is stored in A(1:K+L,N-K-L+1:N) on exit. */
/* > */
/* > If M-K-L < 0, */
/* > */
/* > K M-K K+L-M */
/* > D1 = K ( I 0 0 ) */
/* > M-K ( 0 C 0 ) */
/* > */
/* > K M-K K+L-M */
/* > D2 = M-K ( 0 S 0 ) */
/* > K+L-M ( 0 0 I ) */
/* > P-L ( 0 0 0 ) */
/* > */
/* > N-K-L K M-K K+L-M */
/* > ( 0 R ) = K ( 0 R11 R12 R13 ) */
/* > M-K ( 0 0 R22 R23 ) */
/* > K+L-M ( 0 0 0 R33 ) */
/* > */
/* > where */
/* > C = diag( ALPHA(K+1), ... , ALPHA(M) ), */
/* > S = diag( BETA(K+1), ... , BETA(M) ), */
/* > C**2 + S**2 = I. */
/* > */
/* > R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored */
/* > ( 0 R22 R23 ) */
/* > in B(M-K+1:L,N+M-K-L+1:N) on exit. */
/* > */
/* > The computation of the orthogonal transformation matrices U, V or Q */
/* > is optional. These matrices may either be formed explicitly, or they */
/* > may be postmultiplied into input matrices U1, V1, or Q1. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBU */
/* > \verbatim */
/* > JOBU is CHARACTER*1 */
/* > = 'U': U must contain an orthogonal matrix U1 on entry, and */
/* > the product U1*U is returned;
*/
/* > = 'I': U is initialized to the unit matrix, and the */
/* > orthogonal matrix U is returned;
*/
/* > = 'N': U is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* > JOBV is CHARACTER*1 */
/* > = 'V': V must contain an orthogonal matrix V1 on entry, and */
/* > the product V1*V is returned;
*/
/* > = 'I': V is initialized to the unit matrix, and the */
/* > orthogonal matrix V is returned;
*/
/* > = 'N': V is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBQ */
/* > \verbatim */
/* > JOBQ is CHARACTER*1 */
/* > = 'Q': Q must contain an orthogonal matrix Q1 on entry, and */
/* > the product Q1*Q is returned;
*/
/* > = 'I': Q is initialized to the unit matrix, and the */
/* > orthogonal matrix Q is returned;
*/
/* > = 'N': Q is not computed. */
/* > \endverbatim */
/* > */
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
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is INTEGER */
/* > */
/* > K and L specify the subblocks in the input matrices A and B: */
/* > A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,N-L+1:N) */
/* > of A and B, whose GSVD is going to be computed by STGSJA. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular */
/* > matrix R or part of R. See Purpose for details. */
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
/* > B is REAL array, dimension (LDB,N) */
/* > On entry, the P-by-N matrix B. */
/* > On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains */
/* > a part of R. See Purpose for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in] TOLA */
/* > \verbatim */
/* > TOLA is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] TOLB */
/* > \verbatim */
/* > TOLB is REAL */
/* > */
/* > TOLA and TOLB are the convergence criteria for the Jacobi- */
/* > Kogbetliantz iteration procedure. Generally, they are the */
/* > same as used in the preprocessing step, say */
/* > TOLA = max(M,N)*norm(A)*MACHEPS, */
/* > TOLB = max(P,N)*norm(B)*MACHEPS. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* > ALPHA is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is REAL array, dimension (N) */
/* > */
/* > On exit, ALPHA and BETA contain the generalized singular */
/* > value pairs of A and B;
*/
/* > ALPHA(1:K) = 1, */
/* > BETA(1:K) = 0, */
/* > and if M-K-L >= 0, */
/* > ALPHA(K+1:K+L) = diag(C), */
/* > BETA(K+1:K+L) = diag(S), */
/* > or if M-K-L < 0, */
/* > ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0 */
/* > BETA(K+1:M) = S, BETA(M+1:K+L) = 1. */
/* > Furthermore, if K+L < N, */
/* > ALPHA(K+L+1:N) = 0 and */
/* > BETA(K+L+1:N) = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U */
/* > \verbatim */
/* > U is REAL array, dimension (LDU,M) */
/* > On entry, if JOBU = 'U', U must contain a matrix U1 (usually */
/* > the orthogonal matrix returned by SGGSVP). */
/* > On exit, */
/* > if JOBU = 'I', U contains the orthogonal matrix U;
*/
/* > if JOBU = 'U', U contains the product U1*U. */
/* > If JOBU = 'N', U is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is INTEGER */
/* > The leading dimension of the array U. LDU >= max(1,M) if */
/* > JOBU = 'U';
LDU >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* > V is REAL array, dimension (LDV,P) */
/* > On entry, if JOBV = 'V', V must contain a matrix V1 (usually */
/* > the orthogonal matrix returned by SGGSVP). */
/* > On exit, */
/* > if JOBV = 'I', V contains the orthogonal matrix V;
*/
/* > if JOBV = 'V', V contains the product V1*V. */
/* > If JOBV = 'N', V is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V. LDV >= max(1,P) if */
/* > JOBV = 'V';
LDV >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (LDQ,N) */
/* > On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually */
/* > the orthogonal matrix returned by SGGSVP). */
/* > On exit, */
/* > if JOBQ = 'I', Q contains the orthogonal matrix Q;
*/
/* > if JOBQ = 'Q', Q contains the product Q1*Q. */
/* > If JOBQ = 'N', Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= max(1,N) if */
/* > JOBQ = 'Q';
LDQ >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] NCYCLE */
/* > \verbatim */
/* > NCYCLE is INTEGER */
/* > The number of cycles required for convergence. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > = 1: the procedure does not converge after MAXIT cycles. */
/* > \endverbatim */
/* > */
/* > \verbatim */
/* > Internal Parameters */
/* > =================== */
/* > */
/* > MAXIT INTEGER */
/* > MAXIT specifies the total loops that the iterative procedure */
/* > may take. If after MAXIT cycles, the routine fails to */
/* > converge, we return INFO = 1. */
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
/* > STGSJA essentially uses a variant of Kogbetliantz algorithm to reduce */
/* > min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L */
/* > matrix B13 to the form: */
/* > */
/* > U1**T *A13*Q1 = C1*R1;
V1**T *B13*Q1 = S1*R1, */
/* > */
/* > where U1, V1 and Q1 are orthogonal matrix, and Z**T is the transpose */
/* > of Z. C1 and S1 are diagonal matrices satisfying */
/* > */
/* > C1**2 + S1**2 = I, */
/* > */
/* > and R1 is an L-by-L nonsingular upper triangular matrix. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int stgsja_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k, integer *l, real *a, integer *lda, real *b, integer *ldb, real *tola, real *tolb, real *alpha, real * beta, real *u, integer *ldu, real *v, integer *ldv, real *q, integer * ldq, real *work, integer *ncycle, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    real r__1;
    /* Local variables */
    integer i__, j;
    real a1, a2, a3, b1, b2, b3, csq, csu, csv, snq, rwk, snu, snv;
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *);
    real gamma;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    logical initq, initu, initv, wantq, upper;
    real error, ssmin;
    logical wantu, wantv;
    extern /* Subroutine */
    int scopy_(integer *, real *, integer *, real *, integer *), slags2_(logical *, real *, real *, real *, real *, real *, real *, real *, real *, real *, real *, real *, real *);
    integer kcycle;
    extern /* Subroutine */
    int xerbla_(char *, integer *), slapll_( integer *, real *, integer *, real *, integer *, real *), slartg_( real *, real *, real *, real *, real *), slaset_(char *, integer * , integer *, real *, real *, real *, integer *);
    /* -- LAPACK computational routine (version 3.4.0) -- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Decode and test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alpha;
    --beta;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --work;
    /* Function Body */
    initu = lsame_(jobu, "I");
    wantu = initu || lsame_(jobu, "U");
    initv = lsame_(jobv, "I");
    wantv = initv || lsame_(jobv, "V");
    initq = lsame_(jobq, "I");
    wantq = initq || lsame_(jobq, "Q");
    *info = 0;
    if (! (initu || wantu || lsame_(jobu, "N")))
    {
        *info = -1;
    }
    else if (! (initv || wantv || lsame_(jobv, "N")))
    {
        *info = -2;
    }
    else if (! (initq || wantq || lsame_(jobq, "N")))
    {
        *info = -3;
    }
    else if (*m < 0)
    {
        *info = -4;
    }
    else if (*p < 0)
    {
        *info = -5;
    }
    else if (*n < 0)
    {
        *info = -6;
    }
    else if (*lda < max(1,*m))
    {
        *info = -10;
    }
    else if (*ldb < max(1,*p))
    {
        *info = -12;
    }
    else if (*ldu < 1 || wantu && *ldu < *m)
    {
        *info = -18;
    }
    else if (*ldv < 1 || wantv && *ldv < *p)
    {
        *info = -20;
    }
    else if (*ldq < 1 || wantq && *ldq < *n)
    {
        *info = -22;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("STGSJA", &i__1);
        return 0;
    }
    /* Initialize U, V and Q, if necessary */
    if (initu)
    {
        slaset_("Full", m, m, &c_b13, &c_b14, &u[u_offset], ldu);
    }
    if (initv)
    {
        slaset_("Full", p, p, &c_b13, &c_b14, &v[v_offset], ldv);
    }
    if (initq)
    {
        slaset_("Full", n, n, &c_b13, &c_b14, &q[q_offset], ldq);
    }
    /* Loop until convergence */
    upper = FALSE_;
    for (kcycle = 1;
            kcycle <= 40;
            ++kcycle)
    {
        upper = ! upper;
        i__1 = *l - 1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = *l;
            for (j = i__ + 1;
                    j <= i__2;
                    ++j)
            {
                a1 = 0.f;
                a2 = 0.f;
                a3 = 0.f;
                if (*k + i__ <= *m)
                {
                    a1 = a[*k + i__ + (*n - *l + i__) * a_dim1];
                }
                if (*k + j <= *m)
                {
                    a3 = a[*k + j + (*n - *l + j) * a_dim1];
                }
                b1 = b[i__ + (*n - *l + i__) * b_dim1];
                b3 = b[j + (*n - *l + j) * b_dim1];
                if (upper)
                {
                    if (*k + i__ <= *m)
                    {
                        a2 = a[*k + i__ + (*n - *l + j) * a_dim1];
                    }
                    b2 = b[i__ + (*n - *l + j) * b_dim1];
                }
                else
                {
                    if (*k + j <= *m)
                    {
                        a2 = a[*k + j + (*n - *l + i__) * a_dim1];
                    }
                    b2 = b[j + (*n - *l + i__) * b_dim1];
                }
                slags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, & csv, &snv, &csq, &snq);
                /* Update (K+I)-th and (K+J)-th rows of matrix A: U**T *A */
                if (*k + j <= *m)
                {
                    srot_(l, &a[*k + j + (*n - *l + 1) * a_dim1], lda, &a[*k + i__ + (*n - *l + 1) * a_dim1], lda, &csu, &snu);
                }
                /* Update I-th and J-th rows of matrix B: V**T *B */
                srot_(l, &b[j + (*n - *l + 1) * b_dim1], ldb, &b[i__ + (*n - * l + 1) * b_dim1], ldb, &csv, &snv);
                /* Update (N-L+I)-th and (N-L+J)-th columns of matrices */
                /* A and B: A*Q and B*Q */
                /* Computing MIN */
                i__4 = *k + *l;
                i__3 = min(i__4,*m);
                srot_(&i__3, &a[(*n - *l + j) * a_dim1 + 1], &c__1, &a[(*n - * l + i__) * a_dim1 + 1], &c__1, &csq, &snq);
                srot_(l, &b[(*n - *l + j) * b_dim1 + 1], &c__1, &b[(*n - *l + i__) * b_dim1 + 1], &c__1, &csq, &snq);
                if (upper)
                {
                    if (*k + i__ <= *m)
                    {
                        a[*k + i__ + (*n - *l + j) * a_dim1] = 0.f;
                    }
                    b[i__ + (*n - *l + j) * b_dim1] = 0.f;
                }
                else
                {
                    if (*k + j <= *m)
                    {
                        a[*k + j + (*n - *l + i__) * a_dim1] = 0.f;
                    }
                    b[j + (*n - *l + i__) * b_dim1] = 0.f;
                }
                /* Update orthogonal matrices U, V, Q, if desired. */
                if (wantu && *k + j <= *m)
                {
                    srot_(m, &u[(*k + j) * u_dim1 + 1], &c__1, &u[(*k + i__) * u_dim1 + 1], &c__1, &csu, &snu);
                }
                if (wantv)
                {
                    srot_(p, &v[j * v_dim1 + 1], &c__1, &v[i__ * v_dim1 + 1], &c__1, &csv, &snv);
                }
                if (wantq)
                {
                    srot_(n, &q[(*n - *l + j) * q_dim1 + 1], &c__1, &q[(*n - * l + i__) * q_dim1 + 1], &c__1, &csq, &snq);
                }
                /* L10: */
            }
            /* L20: */
        }
        if (! upper)
        {
            /* The matrices A13 and B13 were lower triangular at the start */
            /* of the cycle, and are now upper triangular. */
            /* Convergence test: test the parallelism of the corresponding */
            /* rows of A and B. */
            error = 0.f;
            /* Computing MIN */
            i__2 = *l;
            i__3 = *m - *k; // , expr subst
            i__1 = min(i__2,i__3);
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                i__2 = *l - i__ + 1;
                scopy_(&i__2, &a[*k + i__ + (*n - *l + i__) * a_dim1], lda, & work[1], &c__1);
                i__2 = *l - i__ + 1;
                scopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &work[* l + 1], &c__1);
                i__2 = *l - i__ + 1;
                slapll_(&i__2, &work[1], &c__1, &work[*l + 1], &c__1, &ssmin);
                error = max(error,ssmin);
                /* L30: */
            }
            if (f2c_abs(error) <= min(*tola,*tolb))
            {
                goto L50;
            }
        }
        /* End of cycle loop */
        /* L40: */
    }
    /* The algorithm has not converged after MAXIT cycles. */
    *info = 1;
    goto L100;
L50: /* If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged. */
    /* Compute the generalized singular value pairs (ALPHA, BETA), and */
    /* set the triangular matrix R to array A. */
    i__1 = *k;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        alpha[i__] = 1.f;
        beta[i__] = 0.f;
        /* L60: */
    }
    /* Computing MIN */
    i__2 = *l;
    i__3 = *m - *k; // , expr subst
    i__1 = min(i__2,i__3);
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        a1 = a[*k + i__ + (*n - *l + i__) * a_dim1];
        b1 = b[i__ + (*n - *l + i__) * b_dim1];
        if (a1 != 0.f)
        {
            gamma = b1 / a1;
            /* change sign if necessary */
            if (gamma < 0.f)
            {
                i__2 = *l - i__ + 1;
                sscal_(&i__2, &c_b43, &b[i__ + (*n - *l + i__) * b_dim1], ldb) ;
                if (wantv)
                {
                    sscal_(p, &c_b43, &v[i__ * v_dim1 + 1], &c__1);
                }
            }
            r__1 = f2c_abs(gamma);
            slartg_(&r__1, &c_b14, &beta[*k + i__], &alpha[*k + i__], &rwk);
            if (alpha[*k + i__] >= beta[*k + i__])
            {
                i__2 = *l - i__ + 1;
                r__1 = 1.f / alpha[*k + i__];
                sscal_(&i__2, &r__1, &a[*k + i__ + (*n - *l + i__) * a_dim1], lda);
            }
            else
            {
                i__2 = *l - i__ + 1;
                r__1 = 1.f / beta[*k + i__];
                sscal_(&i__2, &r__1, &b[i__ + (*n - *l + i__) * b_dim1], ldb);
                i__2 = *l - i__ + 1;
                scopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &a[*k + i__ + (*n - *l + i__) * a_dim1], lda);
            }
        }
        else
        {
            alpha[*k + i__] = 0.f;
            beta[*k + i__] = 1.f;
            i__2 = *l - i__ + 1;
            scopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &a[*k + i__ + (*n - *l + i__) * a_dim1], lda);
        }
        /* L70: */
    }
    /* Post-assignment */
    i__1 = *k + *l;
    for (i__ = *m + 1;
            i__ <= i__1;
            ++i__)
    {
        alpha[i__] = 0.f;
        beta[i__] = 1.f;
        /* L80: */
    }
    if (*k + *l < *n)
    {
        i__1 = *n;
        for (i__ = *k + *l + 1;
                i__ <= i__1;
                ++i__)
        {
            alpha[i__] = 0.f;
            beta[i__] = 0.f;
            /* L90: */
        }
    }
L100:
    *ncycle = kcycle;
    return 0;
    /* End of STGSJA */
}
/* stgsja_ */
