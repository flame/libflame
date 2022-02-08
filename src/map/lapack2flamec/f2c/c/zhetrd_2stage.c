/* ../netlib/v3.9.0/zhetrd_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
/* > \brief \b ZHETRD_2STAGE */
/* @precisions fortran z -> s d c */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHETRD_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrd_ 2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrd_ 2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrd_ 2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHETRD_2STAGE( VECT, UPLO, N, A, LDA, D, E, TAU, */
/* HOUS2, LHOUS2, WORK, LWORK, INFO ) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* CHARACTER VECT, UPLO */
/* INTEGER N, LDA, LWORK, LHOUS2, INFO */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION D( * ), E( * ) */
/* COMPLEX*16 A( LDA, * ), TAU( * ), */
/* HOUS2( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRD_2STAGE reduces a complex Hermitian matrix A to real symmetric */
/* > tridiagonal form T by a unitary similarity transformation: */
/* > Q1**H Q2**H* A * Q2 * Q1 = T. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] VECT */
/* > \verbatim */
/* > VECT is CHARACTER*1 */
/* > = 'N': No need for the Housholder representation, */
/* > in particular for the second stage (Band to */
/* > tridiagonal) and thus LHOUS2 is of size max(1, 4*N);
*/
/* > = 'V': the Householder representation is needed to */
/* > either generate Q1 Q2 or to apply Q1 Q2, */
/* > then LHOUS2 is to be queried and computed. */
/* > (NOT AVAILABLE IN THIS RELEASE). */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
*/
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the leading */
/* > N-by-N upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > On exit, if UPLO = 'U', the band superdiagonal */
/* > of A are overwritten by the corresponding elements of the */
/* > internal band-diagonal matrix AB, and the elements above */
/* > the KD superdiagonal, with the array TAU, represent the unitary */
/* > matrix Q1 as a product of elementary reflectors;
if UPLO */
/* > = 'L', the diagonal and band subdiagonal of A are over- */
/* > written by the corresponding elements of the internal band-diagonal */
/* > matrix AB, and the elements below the KD subdiagonal, with */
/* > the array TAU, represent the unitary matrix Q1 as a product */
/* > of elementary reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is DOUBLE PRECISION array, dimension (N) */
/* > The diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is DOUBLE PRECISION array, dimension (N-1) */
/* > The off-diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX*16 array, dimension (N-KD) */
/* > The scalar factors of the elementary reflectors of */
/* > the first stage (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[out] HOUS2 */
/* > \verbatim */
/* > HOUS2 is COMPLEX*16 array, dimension (LHOUS2) */
/* > Stores the Householder representation of the stage2 */
/* > band to tridiagonal. */
/* > \endverbatim */
/* > */
/* > \param[in] LHOUS2 */
/* > \verbatim */
/* > LHOUS2 is INTEGER */
/* > The dimension of the array HOUS2. */
/* > If LWORK = -1, or LHOUS2 = -1, */
/* > then a query is assumed;
the routine */
/* > only calculates the optimal size of the HOUS2 array, returns */
/* > this value as the first entry of the HOUS2 array, and no error */
/* > message related to LHOUS2 is issued by XERBLA. */
/* > If VECT='N', LHOUS2 = max(1, 4*n);
*/
/* > if VECT='V', option not yet available. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK = MAX(1, dimension) */
/* > If LWORK = -1, or LHOUS2=-1, */
/* > then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > LWORK = MAX(1, dimension) where */
/* > dimension = max(stage1,stage2) + (KD+1)*N */
/* > = N*KD + N*max(KD+1,FACTOPTNB) */
/* > + max(2*KD*KD, KD*NTHREADS) */
/* > + (KD+1)*N */
/* > where KD is the blocking size of the reduction, */
/* > FACTOPTNB is the blocking used by the QR or LQ */
/* > algorithm, usually FACTOPTNB=128 is a good choice */
/* > NTHREADS is the number of threads used when */
/* > openMP compilation is enabled, otherwise =1. */
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
/* > \ingroup complex16HEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Implemented by Azzam Haidar. */
/* > */
/* > All details are available on technical report, SC11, SC13 papers. */
/* > */
/* > Azzam Haidar, Hatem Ltaief, and Jack Dongarra. */
/* > Parallel reduction to condensed forms for symmetric eigenvalue problems */
/* > using aggregated fine-grained and memory-aware kernels. In Proceedings */
/* > of 2011 International Conference for High Performance Computing, */
/* > Networking, Storage and Analysis (SC '11), New York, NY, USA, */
/* > Article 8 , 11 pages. */
/* > http://doi.acm.org/10.1145/2063384.2063394 */
/* > */
/* > A. Haidar, J. Kurzak, P. Luszczek, 2013. */
/* > An improved parallel singular value algorithm and its implementation */
/* > for multicore hardware, In Proceedings of 2013 International Conference */
/* > for High Performance Computing, Networking, Storage and Analysis (SC '13). */
/* > Denver, Colorado, USA, 2013. */
/* > Article 90, 12 pages. */
/* > http://doi.acm.org/10.1145/2503210.2503292 */
/* > */
/* > A. Haidar, R. Solca, S. Tomov, T. Schulthess and J. Dongarra. */
/* > A novel hybrid CPU-GPU generalized eigensolver for electronic structure */
/* > calculations based on fine-grained memory aware tasks. */
/* > International Journal of High Performance Computing Applications. */
/* > Volume 28 Issue 2, Pages 196-209, May 2014. */
/* > http://hpc.sagepub.com/content/28/2/196 */
/* > */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int zhetrd_2stage_(char *vect, char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau, doublecomplex *hous2, integer *lhous2, doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    integer ib, kd, ldab;
    extern /* Subroutine */
    int zhetrd_he2hb_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, integer *);
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int zhetrd_hb2st_(char *, char *, char *, integer *, integer *, doublecomplex *, integer *, doublereal *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    integer lwrk, wpos;
    extern logical lsame_(char *, char *);
    integer abpos, lhmin, lwmin;
    logical wantq, upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    logical lquery;
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    --hous2;
    --work;
    /* Function Body */
    *info = 0;
    wantq = lsame_(vect, "V");
    upper = lsame_(uplo, "U");
    lquery = *lwork == -1 || *lhous2 == -1;
    /* Determine the block size, the workspace size and the hous size. */
    kd = ilaenv2stage_(&c__1, "ZHETRD_2STAGE", vect, n, &c_n1, &c_n1, &c_n1);
    ib = ilaenv2stage_(&c__2, "ZHETRD_2STAGE", vect, n, &kd, &c_n1, &c_n1);
    lhmin = ilaenv2stage_(&c__3, "ZHETRD_2STAGE", vect, n, &kd, &ib, &c_n1);
    lwmin = ilaenv2stage_(&c__4, "ZHETRD_2STAGE", vect, n, &kd, &ib, &c_n1);
    /* WRITE(*,*),'ZHETRD_2STAGE N KD UPLO LHMIN LWMIN ',N, KD, UPLO, */
    /* $ LHMIN, LWMIN */
    if (! lsame_(vect, "N"))
    {
        *info = -1;
    }
    else if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    else if (*lhous2 < lhmin && ! lquery)
    {
        *info = -10;
    }
    else if (*lwork < lwmin && ! lquery)
    {
        *info = -12;
    }
    if (*info == 0)
    {
        hous2[1].r = (doublereal) lhmin;
        hous2[1].i = 0.; // , expr subst
        work[1].r = (doublereal) lwmin;
        work[1].i = 0.; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHETRD_2STAGE", &i__1);
        return 0;
    }
    else if (lquery)
    {
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        work[1].r = 1.;
        work[1].i = 0.; // , expr subst
        return 0;
    }
    /* Determine pointer position */
    ldab = kd + 1;
    lwrk = *lwork - ldab * *n;
    abpos = 1;
    wpos = abpos + ldab * *n;
    zhetrd_he2hb_(uplo, n, &kd, &a[a_offset], lda, &work[abpos], &ldab, &tau[ 1], &work[wpos], &lwrk, info);
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHETRD_HE2HB", &i__1);
        return 0;
    }
    zhetrd_hb2st_("Y", vect, uplo, n, &kd, &work[abpos], &ldab, &d__[1], &e[ 1], &hous2[1], lhous2, &work[wpos], &lwrk, info);
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHETRD_HB2ST", &i__1);
        return 0;
    }
    hous2[1].r = (doublereal) lhmin;
    hous2[1].i = 0.; // , expr subst
    work[1].r = (doublereal) lwmin;
    work[1].i = 0.; // , expr subst
    return 0;
    /* End of ZHETRD_2STAGE */
}
/* zhetrd_2stage__ */
