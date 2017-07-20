/* ../netlib/slaqp2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SLAQP2 computes a QR factorization with column pivoting of the matrix block. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLAQP2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqp2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqp2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqp2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2, */
/* WORK ) */
/* .. Scalar Arguments .. */
/* INTEGER LDA, M, N, OFFSET */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPVT( * ) */
/* REAL A( LDA, * ), TAU( * ), VN1( * ), VN2( * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAQP2 computes a QR factorization with column pivoting of */
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
/* > A is REAL array, dimension (LDA,N) */
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
/* > TAU is REAL array, dimension (min(M,N)) */
/* > The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN1 */
/* > \verbatim */
/* > VN1 is REAL array, dimension (N) */
/* > The vector with the partial column norms. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VN2 */
/* > \verbatim */
/* > VN2 is REAL array, dimension (N) */
/* > The vector with the exact column norms. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (N) */
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
int slaqp2_(integer *m, integer *n, integer *offset, real *a, integer *lda, integer *jpvt, real *tau, real *vn1, real *vn2, real * work)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, mn;
    real aii;
    integer pvt;
    real temp, temp2;
    extern real snrm2_(integer *, real *, integer *);
    real tol3z;
    integer offpi;
    extern /* Subroutine */
    int slarf_(char *, integer *, integer *, real *, integer *, real *, real *, integer *, real *);
    integer itemp;
    extern /* Subroutine */
    int sswap_(integer *, real *, integer *, real *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
    int slarfg_(integer *, real *, real *, integer *, real *);
    extern integer isamax_(integer *, real *, integer *);
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
    tol3z = sqrt(slamch_("Epsilon"));
    /* Compute factorization. */
    i__1 = mn;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        offpi = *offset + i__;
        /* Determine ith pivot column and swap if necessary. */
        i__2 = *n - i__ + 1;
        pvt = i__ - 1 + isamax_(&i__2, &vn1[i__], &c__1);
        if (pvt != i__)
        {
            sswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], & c__1);
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
            slarfg_(&i__2, &a[offpi + i__ * a_dim1], &a[offpi + 1 + i__ * a_dim1], &c__1, &tau[i__]);
        }
        else
        {
            slarfg_(&c__1, &a[*m + i__ * a_dim1], &a[*m + i__ * a_dim1], & c__1, &tau[i__]);
        }
        if (i__ < *n)
        {
            /* Apply H(i)**T to A(offset+i:m,i+1:n) from the left. */
            aii = a[offpi + i__ * a_dim1];
            a[offpi + i__ * a_dim1] = 1.f;
            i__2 = *m - offpi + 1;
            i__3 = *n - i__;
            slarf_("Left", &i__2, &i__3, &a[offpi + i__ * a_dim1], &c__1, & tau[i__], &a[offpi + (i__ + 1) * a_dim1], lda, &work[1]);
            a[offpi + i__ * a_dim1] = aii;
        }
        /* Update partial column norms. */
        i__2 = *n;
        for (j = i__ + 1;
                j <= i__2;
                ++j)
        {
            if (vn1[j] != 0.f)
            {
                /* NOTE: The following 4 lines follow from the analysis in */
                /* Lapack Working Note 176. */
                /* Computing 2nd power */
                r__2 = (r__1 = a[offpi + j * a_dim1], f2c_abs(r__1)) / vn1[j];
                temp = 1.f - r__2 * r__2;
                temp = max(temp,0.f);
                /* Computing 2nd power */
                r__1 = vn1[j] / vn2[j];
                temp2 = temp * (r__1 * r__1);
                if (temp2 <= tol3z)
                {
                    if (offpi < *m)
                    {
                        i__3 = *m - offpi;
                        vn1[j] = snrm2_(&i__3, &a[offpi + 1 + j * a_dim1], & c__1);
                        vn2[j] = vn1[j];
                    }
                    else
                    {
                        vn1[j] = 0.f;
                        vn2[j] = 0.f;
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
    /* End of SLAQP2 */
}
/* slaqp2_ */
