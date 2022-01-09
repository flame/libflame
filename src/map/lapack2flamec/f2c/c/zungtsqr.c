/* ../netlib/v3.9.0/zungtsqr.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
 #include "FLA_f2c.h" /* Table of constant values */
 static doublecomplex c_b1 = {
1.,0.}
;
 static doublecomplex c_b2 = {
0.,0.}
;
 static integer c__1 = 1;
 /* > \brief \b ZUNGTSQR */
 /* =========== DOCUMENTATION =========== */
 /* Online html documentation available at */
 /* http://www.netlib.org/lapack/explore-html/ */
 /* > \htmlonly */
 /* > Download ZUNGTSQR + dependencies */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zuntsqr .f"> */
 /* > [TGZ]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungtsq r.f"> */
 /* > [ZIP]</a> */
 /* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungtsq r.f"> */
 /* > [TXT]</a> */
 /* > */
 /* Definition: */
 /* =========== */
 /* SUBROUTINE ZUNGTSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, */
 /* $ INFO ) */
 /* .. Scalar Arguments .. */
 /* INTEGER INFO, LDA, LDT, LWORK, M, N, MB, NB */
 /* .. */
 /* .. Array Arguments .. */
 /* COMPLEX*16 A( LDA, * ), T( LDT, * ), WORK( * ) */
 /* .. */
 /* > \par Purpose: */
 /* ============= */
 /* > */
 /* > \verbatim */
 /* > */
 /* > ZUNGTSQR generates an M-by-N complex matrix Q_out with orthonormal */
 /* > columns, which are the first N columns of a product of comlpex unitary */
 /* > matrices of order M which are returned by ZLATSQR */
 /* > */
 /* > Q_out = first_N_columns_of( Q(1)_in * Q(2)_in * ... * Q(k)_in ). */
 /* > */
 /* > See the documentation for ZLATSQR. */
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
 /* > The number of columns of the matrix A. M >= N >= 0. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] MB */
 /* > \verbatim */
 /* > MB is INTEGER */
 /* > The row block size used by DLATSQR to return */
 /* > arrays A and T. MB > N. */
 /* > (Note that if MB > M, then M is used instead of MB */
 /* > as the row block size). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] NB */
 /* > \verbatim */
 /* > NB is INTEGER */
 /* > The column block size used by ZLATSQR to return */
 /* > arrays A and T. NB >= 1. */
 /* > (Note that if NB > N, then N is used instead of NB */
 /* > as the column block size). */
 /* > \endverbatim */
 /* > */
 /* > \param[in,out] A */
 /* > \verbatim */
 /* > A is COMPLEX*16 array, dimension (LDA,N) */
 /* > */
 /* > On entry: */
 /* > */
 /* > The elements on and above the diagonal are not accessed. */
 /* > The elements below the diagonal represent the unit */
 /* > lower-trapezoidal blocked matrix V computed by ZLATSQR */
 /* > that defines the input matrices Q_in(k) (ones on the */
 /* > diagonal are not stored) (same format as the output A */
 /* > below the diagonal in ZLATSQR). */
 /* > */
 /* > On exit: */
 /* > */
 /* > The array A contains an M-by-N orthonormal matrix Q_out, */
 /* > i.e the columns of A are orthogonal unit vectors. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDA */
 /* > \verbatim */
 /* > LDA is INTEGER */
 /* > The leading dimension of the array A. LDA >= max(1,M). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] T */
 /* > \verbatim */
 /* > T is COMPLEX*16 array, */
 /* > dimension (LDT, N * NIRB) */
 /* > where NIRB = Number_of_input_row_blocks */
 /* > = MAX( 1, CEIL((M-N)/(MB-N)) ) */
 /* > Let NICB = Number_of_input_col_blocks */
 /* > = CEIL(N/NB) */
 /* > */
 /* > The upper-triangular block reflectors used to define the */
 /* > input matrices Q_in(k), k=(1:NIRB*NICB). The block */
 /* > reflectors are stored in compact form in NIRB block */
 /* > reflector sequences. Each of NIRB block reflector sequences */
 /* > is stored in a larger NB-by-N column block of T and consists */
 /* > of NICB smaller NB-by-NB upper-triangular column blocks. */
 /* > (same format as the output T in ZLATSQR). */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LDT */
 /* > \verbatim */
 /* > LDT is INTEGER */
 /* > The leading dimension of the array T. */
 /* > LDT >= max(1,min(NB1,N)). */
 /* > \endverbatim */
 /* > */
 /* > \param[out] WORK */
 /* > \verbatim */
 /* > (workspace) COMPLEX*16 array, dimension (MAX(2,LWORK)) */
 /* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
 /* > \endverbatim */
 /* > */
 /* > \param[in] LWORK */
 /* > \verbatim */
 /* > The dimension of the array WORK. LWORK >= (M+NB)*N. */
 /* > If LWORK = -1, then a workspace query is assumed. */
 /* > The routine only calculates the optimal size of the WORK */
 /* > array, returns this value as the first entry of the WORK */
 /* > array, and no error message related to LWORK is issued */
 /* > by XERBLA. */
 /* > \endverbatim */
 /* > */
 /* > \param[out] INFO */
 /* > \verbatim */
 /* > INFO is INTEGER */
 /* > = 0: successful exit */
 /* > < 0: if INFO = -i, the i-th argument had an illegal value */
 /* > \endverbatim */
 /* > */
 /* Authors: */
 /* ======== */
 /* > \author Univ. of Tennessee */
 /* > \author Univ. of California Berkeley */
 /* > \author Univ. of Colorado Denver */
 /* > \author NAG Ltd. */
 /* > \date November 2019 */
 /* > \ingroup comlex16OTHERcomputational */
 /* > \par Contributors: */
 /* ================== */
 /* > */
 /* > \verbatim */
 /* > */
 /* > November 2019, Igor Kozachenko, */
 /* > Computer Science Division, */
 /* > University of California, Berkeley */
 /* > */
 /* > \endverbatim */
 /* ===================================================================== */
 /* Subroutine */
 int zungtsqr_(integer *m, integer *n, integer *mb, integer * nb, doublecomplex *a, integer *lda, doublecomplex *t, integer *ldt, doublecomplex *work, integer *lwork, integer *info) {
 /* System generated locals */
 integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2;
 doublecomplex z__1;
 /* Local variables */
 extern /* Subroutine */
 int zlamtsqr_(char *, char *, integer *, integer * , integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
 integer lworkopt, j, lc, lw, ldc, iinfo;
 extern /* Subroutine */
 int zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(char *, integer *), zlaset_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *, integer *);
 logical lquery;
 integer nblocal;
 /* -- LAPACK computational routine (version 3.9.0) -- */
 /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
 /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
 /* November 2019 */
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
 /* .. Intrinsic Functions .. */
 /* .. */
 /* .. Executable Statements .. */
 /* Test the input parameters */
 /* Parameter adjustments */
 a_dim1 = *lda;
 a_offset = 1 + a_dim1;
 a -= a_offset;
 t_dim1 = *ldt;
 t_offset = 1 + t_dim1;
 t -= t_offset;
 --work;
 /* Function Body */
 lquery = *lwork == -1;
 *info = 0;
 if (*m < 0) {
 *info = -1;
 }
 else if (*n < 0 || *m < *n) {
 *info = -2;
 }
 else if (*mb <= *n) {
 *info = -3;
 }
 else if (*nb < 1) {
 *info = -4;
 }
 else if (*lda < max(1,*m)) {
 *info = -6;
 }
 else /* if(complicated condition) */
 {
 /* Computing MAX */
 i__1 = 1; i__2 = min(*nb,*n); // , expr subst  
 if (*ldt < max(i__1,i__2)) {
 *info = -8;
 }
 else {
 /* Test the input LWORK for the dimension of the array WORK. */
 /* This workspace is used to store array C(LDC, N) and WORK(LWORK) */
 /* in the call to ZLAMTSQR. See the documentation for ZLAMTSQR. */
 if (*lwork < 2 && ! lquery) {
 *info = -10;
 }
 else {
 /* Set block size for column blocks */
 nblocal = min(*nb,*n);
 /* LWORK = -1, then set the size for the array C(LDC,N) */
 /* in ZLAMTSQR call and set the optimal size of the work array */
 /* WORK(LWORK) in ZLAMTSQR call. */
 ldc = *m;
 lc = ldc * *n;
 lw = *n * nblocal;
 lworkopt = lc + lw;
 if (*lwork < max(1,lworkopt) && ! lquery) {
 *info = -10;
 }
 }
 }
 }
 /* Handle error in the input parameters and return workspace query. */
 if (*info != 0) {
 i__1 = -(*info);
 xerbla_("ZUNGTSQR", &i__1);
 return 0;
 }
 else if (lquery) {
 z__1.r = (doublereal) lworkopt; z__1.i = 0.; // , expr subst  
 work[1].r = z__1.r; work[1].i = z__1.i; // , expr subst  
 return 0;
 }
 /* Quick return if possible */
 if (min(*m,*n) == 0) {
 z__1.r = (doublereal) lworkopt; z__1.i = 0.; // , expr subst  
 work[1].r = z__1.r; work[1].i = z__1.i; // , expr subst  
 return 0;
 }
 /* (1) Form explicitly the tall-skinny M-by-N left submatrix Q1_in */
 /* of M-by-M orthogonal matrix Q_in, which is implicitly stored in */
 /* the subdiagonal part of input array A and in the input array T. */
 /* Perform by the following operation using the routine ZLAMTSQR. */
 /* Q1_in = Q_in * ( I ), where I is a N-by-N identity matrix, */
 /* ( 0 ) 0 is a (M-N)-by-N zero matrix. */
 /* (1a) Form M-by-N matrix in the array WORK(1:LDC*N) with ones */
 /* on the diagonal and zeros elsewhere. */
 zlaset_("F", m, n, &c_b2, &c_b1, &work[1], &ldc);
 /* (1b) On input, WORK(1:LDC*N) stores ( I );
 */
 /* ( 0 ) */
 /* On output, WORK(1:LDC*N) stores Q1_in. */
 zlamtsqr_("L", "N", m, n, n, mb, &nblocal, &a[a_offset], lda, &t[t_offset] , ldt, &work[1], &ldc, &work[lc + 1], &lw, &iinfo);
 /* (2) Copy the result from the part of the work array (1:M,1:N) */
 /* with the leading dimension LDC that starts at WORK(1) into */
 /* the output array A(1:M,1:N) column-by-column. */
 i__1 = *n;
 for (j = 1;
 j <= i__1;
 ++j) {
 zcopy_(m, &work[(j - 1) * ldc + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
 }
 z__1.r = (doublereal) lworkopt; z__1.i = 0.; // , expr subst  
 work[1].r = z__1.r; work[1].i = z__1.i; // , expr subst  
 return 0;
 /* End of ZUNGTSQR */
 }
 /* zungtsqr_ */
 