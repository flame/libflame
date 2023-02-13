/* ../netlib/sgghrd.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static real c_b10 = 0.f;
static real c_b11 = 1.f;
static integer c__1 = 1;
/* > \brief \b SGGHRD */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGGHRD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgghrd. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgghrd. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgghrd. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, */
/* LDQ, Z, LDZ, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER COMPQ, COMPZ */
/* INTEGER IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ Z( LDZ, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGHRD reduces a pair of real matrices (A,B) to generalized upper */
/* > Hessenberg form using orthogonal transformations, where A is a */
/* > general matrix and B is upper triangular. The form of the */
/* > generalized eigenvalue problem is */
/* > A*x = lambda*B*x, */
/* > and B is typically made upper triangular by computing its QR */
/* > factorization and moving the orthogonal matrix Q to the left side */
/* > of the equation. */
/* > */
/* > This subroutine simultaneously reduces A to a Hessenberg matrix H: */
/* > Q**T*A*Z = H */
/* > and transforms B to another upper triangular matrix T: */
/* > Q**T*B*Z = T */
/* > in order to reduce the problem to its standard form */
/* > H*y = lambda*T*y */
/* > where y = Z**T*x. */
/* > */
/* > The orthogonal matrices Q and Z are determined as products of Givens */
/* > rotations. They may either be formed explicitly, or they may be */
/* > postmultiplied into input matrices Q1 and Z1, so that */
/* > */
/* > Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T */
/* > */
/* > Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T */
/* > */
/* > If Q1 is the orthogonal matrix from the QR factorization of B in the */
/* > original equation A*x = lambda*B*x, then SGGHRD reduces the original */
/* > problem to generalized Hessenberg form. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] COMPQ */
/* > \verbatim */
/* > COMPQ is CHARACTER*1 */
/* > = 'N': do not compute Q;
*/
/* > = 'I': Q is initialized to the unit matrix, and the */
/* > orthogonal matrix Q is returned;
*/
/* > = 'V': Q must contain an orthogonal matrix Q1 on entry, */
/* > and the product Q1*Q is returned. */
/* > \endverbatim */
/* > */
/* > \param[in] COMPZ */
/* > \verbatim */
/* > COMPZ is CHARACTER*1 */
/* > = 'N': do not compute Z;
*/
/* > = 'I': Z is initialized to the unit matrix, and the */
/* > orthogonal matrix Z is returned;
*/
/* > = 'V': Z must contain an orthogonal matrix Z1 on entry, */
/* > and the product Z1*Z is returned. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > */
/* > ILO and IHI mark the rows and columns of A which are to be */
/* > reduced. It is assumed that A is already upper triangular */
/* > in rows and columns 1:ILO-1 and IHI+1:N. ILO and IHI are */
/* > normally set by a previous call to SGGBAL;
otherwise they */
/* > should be set to 1 and N respectively. */
/* > 1 <= ILO <= IHI <= N, if N > 0;
ILO=1 and IHI=0, if N=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, N) */
/* > On entry, the N-by-N general matrix to be reduced. */
/* > On exit, the upper triangle and the first subdiagonal of A */
/* > are overwritten with the upper Hessenberg matrix H, and the */
/* > rest is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB, N) */
/* > On entry, the N-by-N upper triangular matrix B. */
/* > On exit, the upper triangular matrix T = Q**T B Z. The */
/* > elements below the diagonal are set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is REAL array, dimension (LDQ, N) */
/* > On entry, if COMPQ = 'V', the orthogonal matrix Q1, */
/* > typically from the QR factorization of B. */
/* > On exit, if COMPQ='I', the orthogonal matrix Q, and if */
/* > COMPQ = 'V', the product Q1*Q. */
/* > Not referenced if COMPQ='N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. */
/* > LDQ >= N if COMPQ='V' or 'I';
LDQ >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is REAL array, dimension (LDZ, N) */
/* > On entry, if COMPZ = 'V', the orthogonal matrix Z1. */
/* > On exit, if COMPZ='I', the orthogonal matrix Z, and if */
/* > COMPZ = 'V', the product Z1*Z. */
/* > Not referenced if COMPZ='N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > The leading dimension of the array Z. */
/* > LDZ >= N if COMPZ='V' or 'I';
LDZ >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
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
/* > This routine reduces A to Hessenberg and B to triangular form by */
/* > an unblocked reduction, as described in _Matrix_Computations_, */
/* > by Golub and Van Loan (Johns Hopkins Press.) */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int sgghrd_(char *compq, char *compz, integer *n, integer * ilo, integer *ihi, real *a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real *z__, integer *ldz, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, i__1, i__2, i__3;
    /* Local variables */
    real c__, s;
    logical ilq, ilz;
    integer jcol;
    real temp;
    integer jrow;
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    integer icompq;
    extern /* Subroutine */
    int slaset_(char *, integer *, integer *, real *, real *, real *, integer *), slartg_(real *, real *, real * , real *, real *);
    integer icompz;
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
    /* Decode COMPQ */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    /* Function Body */
    if (lsame_(compq, "N"))
    {
        ilq = FALSE_;
        icompq = 1;
    }
    else if (lsame_(compq, "V"))
    {
        ilq = TRUE_;
        icompq = 2;
    }
    else if (lsame_(compq, "I"))
    {
        ilq = TRUE_;
        icompq = 3;
    }
    else
    {
        icompq = 0;
    }
    /* Decode COMPZ */
    if (lsame_(compz, "N"))
    {
        ilz = FALSE_;
        icompz = 1;
    }
    else if (lsame_(compz, "V"))
    {
        ilz = TRUE_;
        icompz = 2;
    }
    else if (lsame_(compz, "I"))
    {
        ilz = TRUE_;
        icompz = 3;
    }
    else
    {
        icompz = 0;
    }
    /* Test the input parameters. */
    *info = 0;
    if (icompq <= 0)
    {
        *info = -1;
    }
    else if (icompz <= 0)
    {
        *info = -2;
    }
    else if (*n < 0)
    {
        *info = -3;
    }
    else if (*ilo < 1)
    {
        *info = -4;
    }
    else if (*ihi > *n || *ihi < *ilo - 1)
    {
        *info = -5;
    }
    else if (*lda < max(1,*n))
    {
        *info = -7;
    }
    else if (*ldb < max(1,*n))
    {
        *info = -9;
    }
    else if (ilq && *ldq < *n || *ldq < 1)
    {
        *info = -11;
    }
    else if (ilz && *ldz < *n || *ldz < 1)
    {
        *info = -13;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGGHRD", &i__1);
        return 0;
    }
    /* Initialize Q and Z if desired. */
    if (icompq == 3)
    {
        slaset_("Full", n, n, &c_b10, &c_b11, &q[q_offset], ldq);
    }
    if (icompz == 3)
    {
        slaset_("Full", n, n, &c_b10, &c_b11, &z__[z_offset], ldz);
    }
    /* Quick return if possible */
    if (*n <= 1)
    {
        return 0;
    }
    /* Zero out lower triangle of B */
    i__1 = *n - 1;
    for (jcol = 1;
            jcol <= i__1;
            ++jcol)
    {
        i__2 = *n;
        for (jrow = jcol + 1;
                jrow <= i__2;
                ++jrow)
        {
            b[jrow + jcol * b_dim1] = 0.f;
            /* L10: */
        }
        /* L20: */
    }
    /* Reduce A and B */
    i__1 = *ihi - 2;
    for (jcol = *ilo;
            jcol <= i__1;
            ++jcol)
    {
        i__2 = jcol + 2;
        for (jrow = *ihi;
                jrow >= i__2;
                --jrow)
        {
            /* Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL) */
            temp = a[jrow - 1 + jcol * a_dim1];
            slartg_(&temp, &a[jrow + jcol * a_dim1], &c__, &s, &a[jrow - 1 + jcol * a_dim1]);
            a[jrow + jcol * a_dim1] = 0.f;
            i__3 = *n - jcol;
            srot_(&i__3, &a[jrow - 1 + (jcol + 1) * a_dim1], lda, &a[jrow + ( jcol + 1) * a_dim1], lda, &c__, &s);
            i__3 = *n + 2 - jrow;
            srot_(&i__3, &b[jrow - 1 + (jrow - 1) * b_dim1], ldb, &b[jrow + ( jrow - 1) * b_dim1], ldb, &c__, &s);
            if (ilq)
            {
                srot_(n, &q[(jrow - 1) * q_dim1 + 1], &c__1, &q[jrow * q_dim1 + 1], &c__1, &c__, &s);
            }
            /* Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1) */
            temp = b[jrow + jrow * b_dim1];
            slartg_(&temp, &b[jrow + (jrow - 1) * b_dim1], &c__, &s, &b[jrow + jrow * b_dim1]);
            b[jrow + (jrow - 1) * b_dim1] = 0.f;
            srot_(ihi, &a[jrow * a_dim1 + 1], &c__1, &a[(jrow - 1) * a_dim1 + 1], &c__1, &c__, &s);
            i__3 = jrow - 1;
            srot_(&i__3, &b[jrow * b_dim1 + 1], &c__1, &b[(jrow - 1) * b_dim1 + 1], &c__1, &c__, &s);
            if (ilz)
            {
                srot_(n, &z__[jrow * z_dim1 + 1], &c__1, &z__[(jrow - 1) * z_dim1 + 1], &c__1, &c__, &s);
            }
            /* L30: */
        }
        /* L40: */
    }
    return 0;
    /* End of SGGHRD */
}
/* sgghrd_ */
