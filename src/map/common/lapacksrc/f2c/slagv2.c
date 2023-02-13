/* ../netlib/slagv2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__2 = 2;
static integer c__1 = 1;
/* > \brief \b SLAGV2 computes the Generalized Schur factorization of a real 2-by-2 matrix pencil (A,B) where B is upper triangular. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAGV2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagv2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagv2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagv2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAGV2( A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, CSL, SNL, */
/* CSR, SNR ) */
/* .. Scalar Arguments .. */
/* INTEGER LDA, LDB */
/* REAL CSL, CSR, SNL, SNR */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), ALPHAI( 2 ), ALPHAR( 2 ), */
/* $ B( LDB, * ), BETA( 2 ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAGV2 computes the Generalized Schur factorization of a real 2-by-2 */
/* > matrix pencil (A,B) where B is upper triangular. This routine */
/* > computes orthogonal (rotation) matrices given by CSL, SNL and CSR, */
/* > SNR such that */
/* > */
/* > 1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0 */
/* > types), then */
/* > */
/* > [ a11 a12 ] := [ CSL SNL ] [ a11 a12 ] [ CSR -SNR ] */
/* > [ 0 a22 ] [ -SNL CSL ] [ a21 a22 ] [ SNR CSR ] */
/* > */
/* > [ b11 b12 ] := [ CSL SNL ] [ b11 b12 ] [ CSR -SNR ] */
/* > [ 0 b22 ] [ -SNL CSL ] [ 0 b22 ] [ SNR CSR ], */
/* > */
/* > 2) if the pencil (A,B) has a pair of complex conjugate eigenvalues, */
/* > then */
/* > */
/* > [ a11 a12 ] := [ CSL SNL ] [ a11 a12 ] [ CSR -SNR ] */
/* > [ a21 a22 ] [ -SNL CSL ] [ a21 a22 ] [ SNR CSR ] */
/* > */
/* > [ b11 0 ] := [ CSL SNL ] [ b11 b12 ] [ CSR -SNR ] */
/* > [ 0 b22 ] [ -SNL CSL ] [ 0 b22 ] [ SNR CSR ] */
/* > */
/* > where b11 >= b22 > 0. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, 2) */
/* > On entry, the 2 x 2 matrix A. */
/* > On exit, A is overwritten by the ``A-part'' of the */
/* > generalized Schur form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > THe leading dimension of the array A. LDA >= 2. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB, 2) */
/* > On entry, the upper triangular 2 x 2 matrix B. */
/* > On exit, B is overwritten by the ``B-part'' of the */
/* > generalized Schur form. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > THe leading dimension of the array B. LDB >= 2. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAR */
/* > \verbatim */
/* > ALPHAR is REAL array, dimension (2) */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* > ALPHAI is REAL array, dimension (2) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is REAL array, dimension (2) */
/* > (ALPHAR(k)+i*ALPHAI(k))/BETA(k) are the eigenvalues of the */
/* > pencil (A,B), k=1,2, i = sqrt(-1). Note that BETA(k) may */
/* > be zero. */
/* > \endverbatim */
/* > */
/* > \param[out] CSL */
/* > \verbatim */
/* > CSL is REAL */
/* > The cosine of the left rotation matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] SNL */
/* > \verbatim */
/* > SNL is REAL */
/* > The sine of the left rotation matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] CSR */
/* > \verbatim */
/* > CSR is REAL */
/* > The cosine of the right rotation matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] SNR */
/* > \verbatim */
/* > SNR is REAL */
/* > The sine of the right rotation matrix. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup realOTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA */
/* ===================================================================== */
/* Subroutine */
int slagv2_(real *a, integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real *beta, real *csl, real *snl, real * csr, real *snr)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset;
    real r__1, r__2, r__3, r__4, r__5, r__6;
    /* Local variables */
    real r__, t, h1, h2, h3, wi, qq, rr, wr1, wr2, ulp;
    extern /* Subroutine */
    int srot_(integer *, real *, integer *, real *, integer *, real *, real *), slag2_(real *, integer *, real *, integer *, real *, real *, real *, real *, real *, real *);
    real anorm, bnorm, scale1, scale2;
    extern /* Subroutine */
    int slasv2_(real *, real *, real *, real *, real * , real *, real *, real *, real *);
    extern real slapy2_(real *, real *);
    real ascale, bscale;
    extern real slamch_(char *);
    real safmin;
    extern /* Subroutine */
    int slartg_(real *, real *, real *, real *, real * );
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alphar;
    --alphai;
    --beta;
    /* Function Body */
    safmin = slamch_("S");
    ulp = slamch_("P");
    /* Scale A */
    /* Computing MAX */
    r__5 = (r__1 = a[a_dim1 + 1], f2c_abs(r__1)) + (r__2 = a[a_dim1 + 2], f2c_abs( r__2));
    r__6 = (r__3 = a[(a_dim1 << 1) + 1], f2c_abs(r__3)) + (r__4 = a[(a_dim1 << 1) + 2], f2c_abs(r__4));
    r__5 = max(r__5,r__6); // ; expr subst
    anorm = max(r__5,safmin);
    ascale = 1.f / anorm;
    a[a_dim1 + 1] = ascale * a[a_dim1 + 1];
    a[(a_dim1 << 1) + 1] = ascale * a[(a_dim1 << 1) + 1];
    a[a_dim1 + 2] = ascale * a[a_dim1 + 2];
    a[(a_dim1 << 1) + 2] = ascale * a[(a_dim1 << 1) + 2];
    /* Scale B */
    /* Computing MAX */
    r__4 = (r__3 = b[b_dim1 + 1], f2c_abs(r__3));
    r__5 = (r__1 = b[(b_dim1 << 1) + 1], f2c_abs(r__1)) + (r__2 = b[(b_dim1 << 1) + 2], f2c_abs(r__2));
    r__4 = max(r__4,r__5); // ; expr subst
    bnorm = max(r__4,safmin);
    bscale = 1.f / bnorm;
    b[b_dim1 + 1] = bscale * b[b_dim1 + 1];
    b[(b_dim1 << 1) + 1] = bscale * b[(b_dim1 << 1) + 1];
    b[(b_dim1 << 1) + 2] = bscale * b[(b_dim1 << 1) + 2];
    /* Check if A can be deflated */
    if ((r__1 = a[a_dim1 + 2], f2c_abs(r__1)) <= ulp)
    {
        *csl = 1.f;
        *snl = 0.f;
        *csr = 1.f;
        *snr = 0.f;
        a[a_dim1 + 2] = 0.f;
        b[b_dim1 + 2] = 0.f;
        wi = 0.f;
        /* Check if B is singular */
    }
    else if ((r__1 = b[b_dim1 + 1], f2c_abs(r__1)) <= ulp)
    {
        slartg_(&a[a_dim1 + 1], &a[a_dim1 + 2], csl, snl, &r__);
        *csr = 1.f;
        *snr = 0.f;
        srot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
        srot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csl, snl);
        a[a_dim1 + 2] = 0.f;
        b[b_dim1 + 1] = 0.f;
        b[b_dim1 + 2] = 0.f;
        wi = 0.f;
    }
    else if ((r__1 = b[(b_dim1 << 1) + 2], f2c_abs(r__1)) <= ulp)
    {
        slartg_(&a[(a_dim1 << 1) + 2], &a[a_dim1 + 2], csr, snr, &t);
        *snr = -(*snr);
        srot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, csr, snr);
        srot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, csr, snr);
        *csl = 1.f;
        *snl = 0.f;
        a[a_dim1 + 2] = 0.f;
        b[b_dim1 + 2] = 0.f;
        b[(b_dim1 << 1) + 2] = 0.f;
        wi = 0.f;
    }
    else
    {
        /* B is nonsingular, first compute the eigenvalues of (A,B) */
        slag2_(&a[a_offset], lda, &b[b_offset], ldb, &safmin, &scale1, & scale2, &wr1, &wr2, &wi);
        if (wi == 0.f)
        {
            /* two real eigenvalues, compute s*A-w*B */
            h1 = scale1 * a[a_dim1 + 1] - wr1 * b[b_dim1 + 1];
            h2 = scale1 * a[(a_dim1 << 1) + 1] - wr1 * b[(b_dim1 << 1) + 1];
            h3 = scale1 * a[(a_dim1 << 1) + 2] - wr1 * b[(b_dim1 << 1) + 2];
            rr = slapy2_(&h1, &h2);
            r__1 = scale1 * a[a_dim1 + 2];
            qq = slapy2_(&r__1, &h3);
            if (rr > qq)
            {
                /* find right rotation matrix to zero 1,1 element of */
                /* (sA - wB) */
                slartg_(&h2, &h1, csr, snr, &t);
            }
            else
            {
                /* find right rotation matrix to zero 2,1 element of */
                /* (sA - wB) */
                r__1 = scale1 * a[a_dim1 + 2];
                slartg_(&h3, &r__1, csr, snr, &t);
            }
            *snr = -(*snr);
            srot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, csr, snr);
            srot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, csr, snr);
            /* compute inf norms of A and B */
            /* Computing MAX */
            r__5 = (r__1 = a[a_dim1 + 1], f2c_abs(r__1)) + (r__2 = a[(a_dim1 << 1) + 1], f2c_abs(r__2));
            r__6 = (r__3 = a[a_dim1 + 2], f2c_abs(r__3) ) + (r__4 = a[(a_dim1 << 1) + 2], f2c_abs(r__4)); // , expr subst
            h1 = max(r__5,r__6);
            /* Computing MAX */
            r__5 = (r__1 = b[b_dim1 + 1], f2c_abs(r__1)) + (r__2 = b[(b_dim1 << 1) + 1], f2c_abs(r__2));
            r__6 = (r__3 = b[b_dim1 + 2], f2c_abs(r__3) ) + (r__4 = b[(b_dim1 << 1) + 2], f2c_abs(r__4)); // , expr subst
            h2 = max(r__5,r__6);
            if (scale1 * h1 >= f2c_abs(wr1) * h2)
            {
                /* find left rotation matrix Q to zero out B(2,1) */
                slartg_(&b[b_dim1 + 1], &b[b_dim1 + 2], csl, snl, &r__);
            }
            else
            {
                /* find left rotation matrix Q to zero out A(2,1) */
                slartg_(&a[a_dim1 + 1], &a[a_dim1 + 2], csl, snl, &r__);
            }
            srot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
            srot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csl, snl);
            a[a_dim1 + 2] = 0.f;
            b[b_dim1 + 2] = 0.f;
        }
        else
        {
            /* a pair of complex conjugate eigenvalues */
            /* first compute the SVD of the matrix B */
            slasv2_(&b[b_dim1 + 1], &b[(b_dim1 << 1) + 1], &b[(b_dim1 << 1) + 2], &r__, &t, snr, csr, snl, csl);
            /* Form (A,B) := Q(A,B)Z**T where Q is left rotation matrix and */
            /* Z is right rotation matrix computed from SLASV2 */
            srot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
            srot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csl, snl);
            srot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, csr, snr);
            srot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, csr, snr);
            b[b_dim1 + 2] = 0.f;
            b[(b_dim1 << 1) + 1] = 0.f;
        }
    }
    /* Unscaling */
    a[a_dim1 + 1] = anorm * a[a_dim1 + 1];
    a[a_dim1 + 2] = anorm * a[a_dim1 + 2];
    a[(a_dim1 << 1) + 1] = anorm * a[(a_dim1 << 1) + 1];
    a[(a_dim1 << 1) + 2] = anorm * a[(a_dim1 << 1) + 2];
    b[b_dim1 + 1] = bnorm * b[b_dim1 + 1];
    b[b_dim1 + 2] = bnorm * b[b_dim1 + 2];
    b[(b_dim1 << 1) + 1] = bnorm * b[(b_dim1 << 1) + 1];
    b[(b_dim1 << 1) + 2] = bnorm * b[(b_dim1 << 1) + 2];
    if (wi == 0.f)
    {
        alphar[1] = a[a_dim1 + 1];
        alphar[2] = a[(a_dim1 << 1) + 2];
        alphai[1] = 0.f;
        alphai[2] = 0.f;
        beta[1] = b[b_dim1 + 1];
        beta[2] = b[(b_dim1 << 1) + 2];
    }
    else
    {
        alphar[1] = anorm * wr1 / scale1 / bnorm;
        alphai[1] = anorm * wi / scale1 / bnorm;
        alphar[2] = alphar[1];
        alphai[2] = -alphai[1];
        beta[1] = 1.f;
        beta[2] = 1.f;
    }
    return 0;
    /* End of SLAGV2 */
}
/* slagv2_ */
