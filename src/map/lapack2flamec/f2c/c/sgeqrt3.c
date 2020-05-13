/* ../netlib/sgeqrt3.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b8 = 1.f;
static real c_b20 = -1.f;
/* > \brief \b SGEQRT3 recursively computes a QR factorization of a general real or complex matrix using the c ompact WY representation of Q. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGEQRT3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqrt3 .f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqrt3 .f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqrt3 .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGEQRT3( M, N, A, LDA, T, LDT, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, M, N, LDT */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), T( LDT, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEQRT3 recursively computes a QR factorization of a real M-by-N */
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
/* > The number of rows of the matrix A. M >= N. */
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
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the real M-by-N matrix A. On exit, the elements on and */
/* > above the diagonal contain the N-by-N upper triangular matrix R;
the */
/* > elements below the diagonal are the columns of V. See below for */
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
/* > T is REAL array, dimension (LDT,N) */
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
/* > \date September 2012 */
/* > \ingroup realGEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix V stores the elementary reflectors H(i) in the i-th column */
/* > below the diagonal. For example, if M=5 and N=3, the matrix V is */
/* > */
/* > V = ( 1 ) */
/* > ( v1 1 ) */
/* > ( v1 v2 1 ) */
/* > ( v1 v2 v3 ) */
/* > ( v1 v2 v3 ) */
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
int sgeqrt3_(integer *m, integer *n, real *a, integer *lda, real *t, integer *ldt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2;
    /* Local variables */
    integer i__, j, i1, j1, n1, n2, iinfo;
    extern /* Subroutine */
    int sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *), strmm_(char *, char *, char *, char *, integer *, integer *, real *, real *, integer *, real *, integer *), xerbla_(char *, integer *), slarfg_(integer *, real *, real *, integer *, real *);
    /* -- LAPACK computational routine (version 3.4.2) -- */
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
    if (*n < 0)
    {
        *info = -2;
    }
    else if (*m < *n)
    {
        *info = -1;
    }
    else if (*lda < max(1,*m))
    {
        *info = -4;
    }
    else if (*ldt < max(1,*n))
    {
        *info = -6;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGEQRT3", &i__1);
        return 0;
    }
    if (*n == 1)
    {
        /* Compute Householder transform when N=1 */
        slarfg_(m, &a[a_offset], &a[min(2,*m) + a_dim1], &c__1, &t[t_offset]);
    }
    else
    {
        /* Otherwise, split A into blocks... */
        n1 = *n / 2;
        n2 = *n - n1;
        /* Computing MIN */
        i__1 = n1 + 1;
        j1 = min(i__1,*n);
        /* Computing MIN */
        i__1 = *n + 1;
        i1 = min(i__1,*m);
        /* Compute A(1:M,1:N1) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H */
        sgeqrt3_(m, &n1, &a[a_offset], lda, &t[t_offset], ldt, &iinfo);
        /* Compute A(1:M,J1:N) = Q1^H A(1:M,J1:N) [workspace: T(1:N1,J1:N)] */
        i__1 = n2;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = n1;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                t[i__ + (j + n1) * t_dim1] = a[i__ + (j + n1) * a_dim1];
            }
        }
        strmm_("L", "L", "T", "U", &n1, &n2, &c_b8, &a[a_offset], lda, &t[j1 * t_dim1 + 1], ldt) ;
        i__1 = *m - n1;
        sgemm_("T", "N", &n1, &n2, &i__1, &c_b8, &a[j1 + a_dim1], lda, &a[j1 + j1 * a_dim1], lda, &c_b8, &t[j1 * t_dim1 + 1], ldt);
        strmm_("L", "U", "T", "N", &n1, &n2, &c_b8, &t[t_offset], ldt, &t[j1 * t_dim1 + 1], ldt) ;
        i__1 = *m - n1;
        sgemm_("N", "N", &i__1, &n2, &n1, &c_b20, &a[j1 + a_dim1], lda, &t[j1 * t_dim1 + 1], ldt, &c_b8, &a[j1 + j1 * a_dim1], lda);
        strmm_("L", "L", "N", "U", &n1, &n2, &c_b8, &a[a_offset], lda, &t[j1 * t_dim1 + 1], ldt) ;
        i__1 = n2;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = n1;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                a[i__ + (j + n1) * a_dim1] -= t[i__ + (j + n1) * t_dim1];
            }
        }
        /* Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H */
        i__1 = *m - n1;
        sgeqrt3_(&i__1, &n2, &a[j1 + j1 * a_dim1], lda, &t[j1 + j1 * t_dim1], ldt, &iinfo);
        /* Compute T3 = T(1:N1,J1:N) = -T1 Y1^H Y2 T2 */
        i__1 = n1;
        for (i__ = 1;
                i__ <= i__1;
                ++i__)
        {
            i__2 = n2;
            for (j = 1;
                    j <= i__2;
                    ++j)
            {
                t[i__ + (j + n1) * t_dim1] = a[j + n1 + i__ * a_dim1];
            }
        }
        strmm_("R", "L", "N", "U", &n1, &n2, &c_b8, &a[j1 + j1 * a_dim1], lda, &t[j1 * t_dim1 + 1], ldt);
        i__1 = *m - *n;
        sgemm_("T", "N", &n1, &n2, &i__1, &c_b8, &a[i1 + a_dim1], lda, &a[i1 + j1 * a_dim1], lda, &c_b8, &t[j1 * t_dim1 + 1], ldt);
        strmm_("L", "U", "N", "N", &n1, &n2, &c_b20, &t[t_offset], ldt, &t[j1 * t_dim1 + 1], ldt);
        strmm_("R", "U", "N", "N", &n1, &n2, &c_b8, &t[j1 + j1 * t_dim1], ldt, &t[j1 * t_dim1 + 1], ldt);
        /* Y = (Y1,Y2);
        R = [ R1 A(1:N1,J1:N) ];
        T = [T1 T3] */
        /* [ 0 R2 ] [ 0 T2] */
    }
    return 0;
    /* End of SGEQRT3 */
}
/* sgeqrt3_ */
