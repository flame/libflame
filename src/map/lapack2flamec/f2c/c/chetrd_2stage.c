/* ../netlib/v3.9.0/chetrd_2stage.f -- translated by f2c (version 20160102). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
/* > \brief \b CHETRD_2STAGE */
/* @generated from zhetrd_2stage.f, fortran z -> c, Sun Nov 6 19:34:06 2016 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHETRD_2STAGE + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetrd_ 2stage.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetrd_ 2stage.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetrd_ 2stage.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHETRD_2STAGE( VECT, UPLO, N, A, LDA, D, E, TAU, */
/* HOUS2, LHOUS2, WORK, LWORK, INFO ) */
/* IMPLICIT NONE */
/* .. Scalar Arguments .. */
/* CHARACTER VECT, UPLO */
/* INTEGER N, LDA, LWORK, LHOUS2, INFO */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), E( * ) */
/* COMPLEX A( LDA, * ), TAU( * ), */
/* HOUS2( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHETRD_2STAGE reduces a complex Hermitian matrix A to real symmetric */
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
/* > tridiagonal) and thus LHOUS2 is of size fla_max(1, 4*N);
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
/* > A is COMPLEX array, dimension (LDA,N) */
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
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > The diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > The off-diagonal elements of the tridiagonal matrix T. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (N-KD) */
/* > The scalar factors of the elementary reflectors of */
/* > the first stage (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[out] HOUS2 */
/* > \verbatim */
/* > HOUS2 is COMPLEX array, dimension (LHOUS2) */
/* > Stores the Householder representation of the stage2 */
/* > band to tridiagonal. */
/* > \endverbatim */
/* > */
/* > \param[in] LHOUS2 */
/* > \verbatim */
/* > LHOUS2 is INTEGER */
/* > The dimension of the array HOUS2. */
/* > If LWORK = -1, or LHOUS2=-1, */
/* > then a query is assumed;
the routine */
/* > only calculates the optimal size of the HOUS2 array, returns */
/* > this value as the first entry of the HOUS2 array, and no error */
/* > message related to LHOUS2 is issued by XERBLA. */
/* > If VECT='N', LHOUS2 = fla_max(1, 4*n);
*/
/* > if VECT='V', option not yet available. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK = MAX(1, dimension) */
/* > If LWORK = -1, or LHOUS2 = -1, */
/* > then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > LWORK = MAX(1, dimension) where */
/* > dimension = fla_max(stage1,stage2) + (KD+1)*N */
/* > = N*KD + N*max(KD+1,FACTOPTNB) */
/* > + fla_max(2*KD*KD, KD*NTHREADS) */
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
/* > \ingroup complexHEcomputational */
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
int chetrd_2stage_(char *vect, char *uplo, integer *n, complex *a, integer *lda, real *d__, real *e, complex *tau, complex * hous2, integer *lhous2, complex *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,"chetrd_2stage inputs: vect %c, uplo %c, n %lld, lda %lld, lhous2 %lld, lwork %lld",*vect, *uplo, *n, *lda, *lhous2, *lwork);
#else
    snprintf(buffer, 256,"chetrd_2stage inputs: vect %c, uplo %c, n %d, lda %d, lhous2 %d, lwork %d",*vect, *uplo, *n, *lda, *lhous2, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    integer ib, kd;
    extern /* Subroutine */
    int chetrd_he2hb_(char *, integer *, integer *, complex *, integer *, complex *, integer *, complex *, complex *, integer *, integer *);
    integer ldab;
    extern /* Subroutine */
    int chetrd_hb2st_(char *, char *, char *, integer *, integer *, complex *, integer *, real *, real *, complex *, integer *, complex *, integer *, integer *);
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer lwrk, wpos;
    extern logical lsame_(char *, char *);
    integer abpos, lhmin, lwmin;
    logical upper;
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
    upper = lsame_(uplo, "U");
    lquery = *lwork == -1 || *lhous2 == -1;
    /* Determine the block size, the workspace size and the hous size. */
    kd = ilaenv2stage_(&c__1, "CHETRD_2STAGE", vect, n, &c_n1, &c_n1, &c_n1);
    ib = ilaenv2stage_(&c__2, "CHETRD_2STAGE", vect, n, &kd, &c_n1, &c_n1);
    lhmin = ilaenv2stage_(&c__3, "CHETRD_2STAGE", vect, n, &kd, &ib, &c_n1);
    lwmin = ilaenv2stage_(&c__4, "CHETRD_2STAGE", vect, n, &kd, &ib, &c_n1);
    /* WRITE(*,*),'CHETRD_2STAGE N KD UPLO LHMIN LWMIN ',N, KD, UPLO, */
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
    else if (*lda < fla_max(1,*n))
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
        hous2[1].r = (real) lhmin;
        hous2[1].i = 0.f; // , expr subst
        work[1].r = (real) lwmin;
        work[1].i = 0.f; // , expr subst
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHETRD_2STAGE", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    else if (lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Quick return if possible */
    if (*n == 0)
    {
        work[1].r = 1.f;
        work[1].i = 0.f; // , expr subst
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    /* Determine pointer position */
    ldab = kd + 1;
    lwrk = *lwork - ldab * *n;
    abpos = 1;
    wpos = abpos + ldab * *n;
    chetrd_he2hb_(uplo, n, &kd, &a[a_offset], lda, &work[abpos], &ldab, &tau[ 1], &work[wpos], &lwrk, info);
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHETRD_HE2HB", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    chetrd_hb2st_("Y", vect, uplo, n, &kd, &work[abpos], &ldab, &d__[1], &e[ 1], &hous2[1], lhous2, &work[wpos], &lwrk, info);
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CHETRD_HB2ST", &i__1);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return 0;
    }
    hous2[1].r = (real) lhmin;
    hous2[1].i = 0.f; // , expr subst
    work[1].r = (real) lwmin;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CHETRD_2STAGE */
}
/* chetrd_2stage__ */
