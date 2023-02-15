/* ../netlib/zhetri2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static integer c_n1 = -1;
/* > \brief \b ZHETRI2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHETRI2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri2 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri2 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri2 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHETRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETRI2 computes the inverse of a COMPLEX*16 hermitian indefinite matrix */
/* > A using the factorization A = U*D*U**T or A = L*D*L**T computed by */
/* > ZHETRF. ZHETRI2 set the LEADING DIMENSION of the workspace */
/* > before calling ZHETRI2X that actually computes the inverse. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**T;
*/
/* > = 'L': Lower triangular, form is A = L*D*L**T. */
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
/* > On entry, the NB diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by ZHETRF. */
/* > */
/* > On exit, if INFO = 0, the (symmetric) inverse of the original */
/* > matrix. If UPLO = 'U', the upper triangular part of the */
/* > inverse is formed and the part of A below the diagonal is not */
/* > referenced;
if UPLO = 'L' the lower triangular part of the */
/* > inverse is formed and the part of A above the diagonal is */
/* > not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the NB structure of D */
/* > as determined by ZHETRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N+NB+1)*(NB+3) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > WORK is size >= (N+NB+1)*(NB+3) */
/* > If LDWORK = -1, then a workspace query is assumed;
the routine */
/* > calculates: */
/* > - the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, */
/* > - and no error message related to LDWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, D(i,i) = 0;
the matrix is singular and its */
/* > inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complex16HEcomputational */
/* ===================================================================== */
/* Subroutine */
int zhetri2_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    /* Local variables */
    extern /* Subroutine */
    int zhetri2x_(char *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *, integer *);
    extern logical lsame_(char *, char *);
    integer nbmax;
    logical upper;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern /* Subroutine */
    int zhetri_(char *, integer *, doublecomplex *, integer *, integer *, doublecomplex *, integer *);
    logical lquery;
    integer minsize;
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    lquery = *lwork == -1;
    /* Get blocksize */
    nbmax = ilaenv_(&c__1, "ZHETRF", uplo, n, &c_n1, &c_n1, &c_n1);
    if (nbmax >= *n)
    {
        minsize = *n;
    }
    else
    {
        minsize = (*n + nbmax + 1) * (nbmax + 3);
    }
    if (! upper && ! lsame_(uplo, "L"))
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < max(1,*n))
    {
        *info = -4;
    }
    else if (*lwork < minsize && ! lquery)
    {
        *info = -7;
    }
    /* Quick return if possible */
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHETRI2", &i__1);
        return 0;
    }
    else if (lquery)
    {
        work[1].r = (doublereal) minsize;
        work[1].i = 0.; // , expr subst
        return 0;
    }
    if (*n == 0)
    {
        return 0;
    }
    if (nbmax >= *n)
    {
        zhetri_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], info);
    }
    else
    {
        zhetri2x_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], &nbmax, info);
    }
    return 0;
    /* End of ZHETRI2 */
}
/* zhetri2_ */
