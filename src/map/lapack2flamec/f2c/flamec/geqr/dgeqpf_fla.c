/* dgeqpf.f -- translated by f2c (version 20190311). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b DGEQPF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGEQPF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqpf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqpf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqpf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPVT( * ) */
/* DOUBLE PRECISION A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine DGEQP3. */
/* > */
/* > DGEQPF computes a QR factorization with column pivoting of a */
/* > real M-by-N matrix A: A*P = Q*R. */
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
/* > The number of columns of the matrix A. N >= 0 */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the upper triangle of the array contains the */
/* > fla_min(M,N)-by-N upper triangular matrix R;
the elements */
/* > below the diagonal, together with the array TAU, */
/* > represent the orthogonal matrix Q as a product of */
/* > fla_min(m,n) elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
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
/* > TAU is DOUBLE PRECISION array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (3*N) */
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
/* > \ingroup doubleGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of elementary reflectors */
/* > */
/* > Q = H(1) H(2) . . . H(n) */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H = I - tau * v * v**T */
/* > */
/* > where tau is a real scalar, and v is a real vector with */
/* > v(1:i-1) = 0 and v(i) = 1;
v(i+1:m) is stored on exit in A(i+1:m,i). */
/* > */
/* > The matrix P is represented in jpvt as follows: If */
/* > jpvt(j) = i */
/* > then the jth column of P is the ith canonical unit vector. */
/* > */
/* > Partial column norm updating strategy modified by */
/* > Z. Drmac and Z. Bujanovic, Dept. of Mathematics, */
/* > University of Zagreb, Croatia. */
/* > -- April 2011 -- */
/* > For more details see LAPACK Working Note 176. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
int dgeqpf_fla(integer *m, integer *n, doublereal *a, integer * lda, integer *jpvt, doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, ma, mn;
    doublereal aii;
    integer pvt;
    doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    doublereal temp2, tol3z;
    extern /* Subroutine */
    int dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *);
    integer itemp;
    extern /* Subroutine */
    int dswap_(integer *, doublereal *, integer *, doublereal *, integer *), dgeqr2_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), dorm2r_(char *, char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
    int dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
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
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -4;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGEQPF", &i__1);
        return 0;
    }
    mn = fla_min(*m,*n);
    tol3z = sqrt(dlamch_("Epsilon"));
    /* Move initial columns up front */
    itemp = 1;
    i__1 = *n;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        if (jpvt[i__] != 0)
        {
            if (i__ != itemp)
            {
                dswap_(m, &a[i__ * a_dim1 + 1], &c__1, &a[itemp * a_dim1 + 1], &c__1);
                jpvt[i__] = jpvt[itemp];
                jpvt[itemp] = i__;
            }
            else
            {
                jpvt[i__] = i__;
            }
            ++itemp;
        }
        else
        {
            jpvt[i__] = i__;
        }
        /* L10: */
    }
    --itemp;
    /* Compute the QR factorization and update remaining columns */
    if (itemp > 0)
    {
        ma = fla_min(itemp,*m);
        dgeqr2_(m, &ma, &a[a_offset], lda, &tau[1], &work[1], info);
        if (ma < *n)
        {
            i__1 = *n - ma;
            dorm2r_("Left", "Transpose", m, &i__1, &ma, &a[a_offset], lda, & tau[1], &a[(ma + 1) * a_dim1 + 1], lda, &work[1], info);
        }
    }
    if (itemp < mn)
    {
        /* Initialize partial column norms. The first n elements of */
        /* work store the exact column norms. */
        i__1 = *n;
        for (i__ = itemp + 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = *m - itemp;
            work[i__] = dnrm2_(&i__2, &a[itemp + 1 + i__ * a_dim1], &c__1);
            work[*n + i__] = work[i__];
            /* L20: */
        }
        /* Compute factorization */
        i__1 = mn;
        for (i__ = itemp + 1;
                i__ <= i__1;
                ++i__)
        {
            /* Determine ith pivot column and swap if necessary */
            i__2 = *n - i__ + 1;
            pvt = i__ - 1 + idamax_(&i__2, &work[i__], &c__1);
            if (pvt != i__)
            {
                dswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], & c__1);
                itemp = jpvt[pvt];
                jpvt[pvt] = jpvt[i__];
                jpvt[i__] = itemp;
                work[pvt] = work[i__];
                work[*n + pvt] = work[*n + i__];
            }
            /* Generate elementary reflector H(i) */
            if (i__ < *m)
            {
                i__2 = *m - i__ + 1;
                dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[i__]);
            }
            else
            {
                dlarfg_(&c__1, &a[*m + *m * a_dim1], &a[*m + *m * a_dim1], & c__1, &tau[*m]);
            }
            if (i__ < *n)
            {
                /* Apply H(i) to A(i:m,i+1:n) from the left */
                aii = a[i__ + i__ * a_dim1];
                a[i__ + i__ * a_dim1] = 1.;
                i__2 = *m - i__ + 1;
                i__3 = *n - i__;
                dlarf_("LEFT", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, & tau[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[(* n << 1) + 1]);
                a[i__ + i__ * a_dim1] = aii;
            }
            /* Update partial column norms */
            i__2 = *n;
            for (j = i__ + 1;
                    j <= i__2;
                    ++j)
            {
                if (work[j] != 0.)
                {
                    /* NOTE: The following 4 lines follow from the analysis in */
                    /* Lapack Working Note 176. */
                    temp = (d__1 = a[i__ + j * a_dim1], f2c_abs(d__1)) / work[j];
                    /* Computing MAX */
                    d__1 = 0.;
                    d__2 = (temp + 1.) * (1. - temp); // , expr subst
                    temp = fla_max(d__1,d__2);
                    /* Computing 2nd power */
                    d__1 = work[j] / work[*n + j];
                    temp2 = temp * (d__1 * d__1);
                    if (temp2 <= tol3z)
                    {
                        if (*m - i__ > 0)
                        {
                            i__3 = *m - i__;
                            work[j] = dnrm2_(&i__3, &a[i__ + 1 + j * a_dim1], &c__1);
                            work[*n + j] = work[j];
                        }
                        else
                        {
                            work[j] = 0.;
                            work[*n + j] = 0.;
                        }
                    }
                    else
                    {
                        work[j] *= sqrt(temp);
                    }
                }
                /* L30: */
            }
            /* L40: */
        }
    }
    return 0;
    /* End of DGEQPF */
}
/* dgeqpf_ */
