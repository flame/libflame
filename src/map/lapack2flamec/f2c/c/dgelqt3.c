/* ../netlib/v3.9.0/dgelqt3.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static doublereal c_b7 = 1.;
 static doublereal c_b19 = -1.;
 /* > \brief \b DGELQT3 recursively computes a LQ factorization of a general real or complex matrix using the c ompact WY representation of Q. */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download DGEQRT3 + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelqt3 .f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelqt3 .f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelqt3 .f"> */
 /* > [TXT]</a> */
 /* > \endhtmlonly */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE DGELQT3( M, N, A, LDA, T, LDT, INFO ) */
 /* .. Scalar Arguments .. */
 /* INTEGER INFO, LDA, M, N, LDT */
 /* .. */
 /* .. Array Arguments .. */
 /* DOUBLE PRECISION A( LDA, * ), T( LDT, * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > DGELQT3 recursively computes a LQ factorization of a real M-by-N */
 /* > matrix A, using the compact WY representation of Q. */
 /* > */
 /* > Based on the algorithm of Elmroth and Gustavson, */
 /* > IBM J. Res. Develop. Vol 44 No. 4 July 2000. */
 /* > \endverbatim */
 /* Arguments: */
 /* ========== */
 /* > \param[in] M */
 /* > \verbatim */
 /* > M is INTEGER */
 /* > The number of rows of the matrix A. M =< N. */
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
 /* > A is DOUBLE PRECISION array, dimension (LDA,N) */
 /* > On entry, the real M-by-N matrix A. On exit, the elements on and */
 /* > below the diagonal contain the N-by-N lower triangular matrix L;
 the */
 /* > elements above the diagonal are the rows of V. See below for */
 /* > further details. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] T */
 /* > \verbatim */
 /* > T is DOUBLE PRECISION array, dimension (LDT,N) */
 /* > The N-by-N upper triangular factor of the block reflector. */
 /* > The elements on and above the diagonal contain the block */
 /* > reflector T;
 the elements below the diagonal are not used. */
 /* > See below for further details. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDT */
 /* > \verbatim */
 /* > LDT is INTEGER */
 /* > The leading dimension of the array T. LDT >= max(1,N). */
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
 /* > \date November 2017 */
 /* > \ingroup doubleGEcomputational */
 /* > \par Further Details: */
 /* ===================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > The matrix V stores the elementary reflectors H(i) in the i-th row */
 /* > above the diagonal. For example, if M=5 and N=3, the matrix V is */
 /* > */
 /* > V = ( 1 v1 v1 v1 v1 ) */
 /* > ( 1 v2 v2 v2 ) */
 /* > ( 1 v3 v3 v3 ) */
 /* > */
 /* > */
 /* > where the vi's represent the vectors which define H(i), which are returned */
 /* > in the matrix A. The 1's along the diagonal of V are not stored in A. The */
 /* > block reflector H is then given by */
 /* > */
 /* > H = I - V * T * V**T */
 /* > */
 /* > where V**T is the transpose of V. */
 /* > */
 /* > For details of the algorithm, see Elmroth and Gustavson (cited above). */
 /* > \endverbatim */
 /* > */
 /* ===================================================================== */
 /* Subroutine */
 int dgelqt3_(integer *m, integer *n, doublereal *a, integer * lda, doublereal *t, integer *ldt, integer *info) {
 AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE 
 char buffer[256]; 
 snprintf(buffer, 256,"dgelqt3 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldt %" FLA_IS "",*m, *n, *lda, *ldt);
 AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
 /* System generated locals */
 integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2;
 /* Local variables */
 integer i__, j, i1, j1, m1, m2;
 extern /* Subroutine */
 int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
 integer iinfo;
 extern /* Subroutine */
 int dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *), dlarfg_( integer *, doublereal *, doublereal *, integer *, doublereal *), xerbla_(char *, integer *);
 /* -- LAPACK computational routine (version 3.8.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* November 2017 */
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
 /* .. Executable Statements .. */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 t_dim1 = *ldt;
 t_offset = 1 + t_dim1;
 t -= t_offset;
 /* Function Body */
 *info = 0;
 if (*m < 0) {
 *info = -1;
 }
 else if (*n < *m) {
 *info = -2;
 }
 else if (*lda < max(1,*m)) {
 *info = -4;
 }
 else if (*ldt < max(1,*m)) {
 *info = -6;
 }
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("DGELQT3", &i__1);
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 }
 if (*m == 1) {
 /* Compute Householder transform when N=1 */
 dlarfg_(n, &a[a_offset], &a[min(2,*n) * a_dim1 + 1], lda, &t[t_offset] );
 }
 else {
 /* Otherwise, split A into blocks... */
 m1 = *m / 2;
 m2 = *m - m1;
 /* Computing MIN */
 i__1 = m1 + 1;
 i1 = min(i__1,*m);
 /* Computing MIN */
 i__1 = *m + 1;
 j1 = min(i__1,*n);
 /* Compute A(1:M1,1:N) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H */
 dgelqt3_(&m1, n, &a[a_offset], lda, &t[t_offset], ldt, &iinfo);
 /* Compute A(J1:M,1:N) = Q1^H A(J1:M,1:N) [workspace: T(1:N1,J1:N)] */
 i__1 = m2;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = m1;
 for (j = 1;
 j <= i__2;
 ++j) {
 t[i__ + m1 + j * t_dim1] = a[i__ + m1 + j * a_dim1];
 }
 }
 dtrmm_("R", "U", "T", "U", &m2, &m1, &c_b7, &a[a_offset], lda, &t[i1 + t_dim1], ldt);
 i__1 = *n - m1;
 dgemm_("N", "T", &m2, &m1, &i__1, &c_b7, &a[i1 + i1 * a_dim1], lda, & a[i1 * a_dim1 + 1], lda, &c_b7, &t[i1 + t_dim1], ldt);
 dtrmm_("R", "U", "N", "N", &m2, &m1, &c_b7, &t[t_offset], ldt, &t[i1 + t_dim1], ldt);
 i__1 = *n - m1;
 dgemm_("N", "N", &m2, &i__1, &m1, &c_b19, &t[i1 + t_dim1], ldt, &a[i1 * a_dim1 + 1], lda, &c_b7, &a[i1 + i1 * a_dim1], lda);
 dtrmm_("R", "U", "N", "U", &m2, &m1, &c_b7, &a[a_offset], lda, &t[i1 + t_dim1], ldt);
 i__1 = m2;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = m1;
 for (j = 1;
 j <= i__2;
 ++j) {
 a[i__ + m1 + j * a_dim1] -= t[i__ + m1 + j * t_dim1];
 t[i__ + m1 + j * t_dim1] = 0.;
 }
 }
 /* Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H */
 i__1 = *n - m1;
 dgelqt3_(&m2, &i__1, &a[i1 + i1 * a_dim1], lda, &t[i1 + i1 * t_dim1], ldt, &iinfo);
 /* Compute T3 = T(J1:N1,1:N) = -T1 Y1^H Y2 T2 */
 i__1 = m2;
 for (i__ = 1;
 i__ <= i__1;
 ++i__) {
 i__2 = m1;
 for (j = 1;
 j <= i__2;
 ++j) {
 t[j + (i__ + m1) * t_dim1] = a[j + (i__ + m1) * a_dim1];
 }
 }
 dtrmm_("R", "U", "T", "U", &m1, &m2, &c_b7, &a[i1 + i1 * a_dim1], lda, &t[i1 * t_dim1 + 1], ldt);
 i__1 = *n - *m;
 dgemm_("N", "T", &m1, &m2, &i__1, &c_b7, &a[j1 * a_dim1 + 1], lda, &a[ i1 + j1 * a_dim1], lda, &c_b7, &t[i1 * t_dim1 + 1], ldt);
 dtrmm_("L", "U", "N", "N", &m1, &m2, &c_b19, &t[t_offset], ldt, &t[i1 * t_dim1 + 1], ldt);
 dtrmm_("R", "U", "N", "N", &m1, &m2, &c_b7, &t[i1 + i1 * t_dim1], ldt, &t[i1 * t_dim1 + 1], ldt);
 /* Y = (Y1,Y2);
 L = [ L1 0 ];
 T = [T1 T3] */
 /* [ A(1:N1,J1:N) L2 ] [ 0 T2] */
 }
 AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
 return 0;
 /* End of DGELQT3 */
 }
 /* dgelqt3_ */
 