/* ../netlib/dlaqp2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DLAQP2 computes a QR factorization with column pivoting of the matrix block. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAQP2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqp2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqp2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqp2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2, */
/* WORK ) */
/* .. Scalar Arguments .. */
/* INTEGER LDA, M, N, OFFSET */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPVT( * ) */
/* DOUBLE PRECISION A( LDA, * ), TAU( * ), VN1( * ), VN2( * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQP2 computes a QR factorization with column pivoting of */
/* > the block A(OFFSET+1:M,1:N). */
/* > The block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized. */
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
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] OFFSET */
/* > \verbatim */
/* > OFFSET is INTEGER */
/* > The number of rows of the matrix A that must be pivoted */
/* > but no factorized. OFFSET >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the upper triangle of block A(OFFSET+1:M,1:N) is */
/* > the triangular factor obtained;
the elements in block */
/* > A(OFFSET+1:M,1:N) below the diagonal, together with the */
/* > array TAU, represent the orthogonal matrix Q as a product of */
/* > elementary reflectors. Block A(1:OFFSET,1:N) has been */
/* > accordingly pivoted, but no factorized. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] JPVT */
/* > \verbatim */
/* > JPVT is INTEGER array, dimension (N) */
/* > On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted */
/* > to the front of A*P (a leading column);
if JPVT(i) = 0, */
/* > the i-th column of A is a free column. */
/* > On exit, if JPVT(i) = k, then the i-th column of A*P */
/* > was the k-th column of A. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (min(M,N)) */
/* > The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN1 */
/* > \verbatim */
/* > VN1 is DOUBLE PRECISION array, dimension (N) */
/* > The vector with the partial column norms. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN2 */
/* > \verbatim */
/* > VN2 is DOUBLE PRECISION array, dimension (N) */
/* > The vector with the exact column norms. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup doubleOTHERauxiliary */
/* > \par Contributors: */
/* ================== */
/* > */
/* > G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain */
/* > X. Sun, Computer Science Dept., Duke University, USA */
/* > \n */
/* > Partial column norm updating strategy modified on April 2011 */
/* > Z. Drmac and Z. Bujanovic, Dept. of Mathematics, */
/* > University of Zagreb, Croatia. */
/* > \par References: */
/* ================ */
/* > */
/* > LAPACK Working Note 176 */
/* > \htmlonly */
/* > <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">[PDF]</a> */
/* > \endhtmlonly */
/* ===================================================================== */
/* Subroutine */
int dlaqp2_(integer *m, integer *n, integer *offset, doublereal *a, integer *lda, integer *jpvt, doublereal *tau, doublereal *vn1, doublereal *vn2, doublereal *work)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, mn;
    doublereal aii;
    integer pvt;
    doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    doublereal temp2, tol3z;
    extern /* Subroutine */
    int dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *);
    integer offpi, itemp;
    extern /* Subroutine */
    int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
    int dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    /* -- LAPACK auxiliary routine (version 3.5.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2013 */
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
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --tau;
    --vn1;
    --vn2;
    --work;
    /* Function Body */
    /* Computing MIN */
    i__1 = *m - *offset;
    mn = min(i__1,*n);
    tol3z = sqrt(dlamch_("Epsilon"));
    /* Compute factorization. */
    i__1 = mn;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        offpi = *offset + i__;
        /* Determine ith pivot column and swap if necessary. */
        i__2 = *n - i__ + 1;
        pvt = i__ - 1 + idamax_(&i__2, &vn1[i__], &c__1);
        if (pvt != i__)
        {
            dswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], & c__1);
            itemp = jpvt[pvt];
            jpvt[pvt] = jpvt[i__];
            jpvt[i__] = itemp;
            vn1[pvt] = vn1[i__];
            vn2[pvt] = vn2[i__];
        }
        /* Generate elementary reflector H(i). */
        if (offpi < *m)
        {
            i__2 = *m - offpi + 1;
            dlarfg_(&i__2, &a[offpi + i__ * a_dim1], &a[offpi + 1 + i__ * a_dim1], &c__1, &tau[i__]);
        }
        else
        {
            dlarfg_(&c__1, &a[*m + i__ * a_dim1], &a[*m + i__ * a_dim1], & c__1, &tau[i__]);
        }
        if (i__ < *n)
        {
            /* Apply H(i)**T to A(offset+i:m,i+1:n) from the left. */
            aii = a[offpi + i__ * a_dim1];
            a[offpi + i__ * a_dim1] = 1.;
            i__2 = *m - offpi + 1;
            i__3 = *n - i__;
            dlarf_("Left", &i__2, &i__3, &a[offpi + i__ * a_dim1], &c__1, & tau[i__], &a[offpi + (i__ + 1) * a_dim1], lda, &work[1]);
            a[offpi + i__ * a_dim1] = aii;
        }
        /* Update partial column norms. */
        i__2 = *n;
        for (j = i__ + 1;
                j <= i__2;
                ++j)
        {
            if (vn1[j] != 0.)
            {
                /* NOTE: The following 4 lines follow from the analysis in */
                /* Lapack Working Note 176. */
                /* Computing 2nd power */
                d__2 = (d__1 = a[offpi + j * a_dim1], f2c_abs(d__1)) / vn1[j];
                temp = 1. - d__2 * d__2;
                temp = max(temp,0.);
                /* Computing 2nd power */
                d__1 = vn1[j] / vn2[j];
                temp2 = temp * (d__1 * d__1);
                if (temp2 <= tol3z)
                {
                    if (offpi < *m)
                    {
                        i__3 = *m - offpi;
                        vn1[j] = dnrm2_(&i__3, &a[offpi + 1 + j * a_dim1], & c__1);
                        vn2[j] = vn1[j];
                    }
                    else
                    {
                        vn1[j] = 0.;
                        vn2[j] = 0.;
                    }
                }
                else
                {
                    vn1[j] *= sqrt(temp);
                }
            }
            /* L10: */
        }
        /* L20: */
    }
    return 0;
    /* End of DLAQP2 */
}
/* dlaqp2_ */
