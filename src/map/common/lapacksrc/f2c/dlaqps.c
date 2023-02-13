/* ../netlib/dlaqps.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static doublereal c_b8 = -1.;
static doublereal c_b9 = 1.;
static doublereal c_b16 = 0.;
/* > \brief \b DLAQPS computes a step of QR factorization with column pivoting of a real m-by-n matrix A by us ing BLAS level 3. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLAQPS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqps. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqps. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqps. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1, */
/* VN2, AUXV, F, LDF ) */
/* .. Scalar Arguments .. */
/* INTEGER KB, LDA, LDF, M, N, NB, OFFSET */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPVT( * ) */
/* DOUBLE PRECISION A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ), */
/* $ VN1( * ), VN2( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAQPS computes a step of QR factorization with column pivoting */
/* > of a real M-by-N matrix A by using Blas-3. It tries to factorize */
/* > NB columns from A starting from the row OFFSET+1, and updates all */
/* > of the matrix with Blas-3 xGEMM. */
/* > */
/* > In some cases, due to catastrophic cancellations, it cannot */
/* > factorize NB columns. Hence, the actual number of factorized */
/* > columns is returned in KB. */
/* > */
/* > Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized. */
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
/* > \param[in] OFFSET */
/* > \verbatim */
/* > OFFSET is INTEGER */
/* > The number of rows of A that have been factorized in */
/* > previous steps. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The number of columns to factorize. */
/* > \endverbatim */
/* > */
/* > \param[out] KB */
/* > \verbatim */
/* > KB is INTEGER */
/* > The number of columns actually factorized. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, block A(OFFSET+1:M,1:KB) is the triangular */
/* > factor obtained and block A(1:OFFSET,1:N) has been */
/* > accordingly pivoted, but no factorized. */
/* > The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has */
/* > been updated. */
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
/* > JPVT(I) = K <==> Column K of the full matrix A has been */
/* > permuted into position I in AP. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION array, dimension (KB) */
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
/* > \param[in,out] AUXV */
/* > \verbatim */
/* > AUXV is DOUBLE PRECISION array, dimension (NB) */
/* > Auxiliar vector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] F */
/* > \verbatim */
/* > F is DOUBLE PRECISION array, dimension (LDF,NB) */
/* > Matrix F**T = L*Y**T*A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDF */
/* > \verbatim */
/* > LDF is INTEGER */
/* > The leading dimension of the array F. LDF >= max(1,N). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
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
int dlaqps_(integer *m, integer *n, integer *offset, integer *nb, integer *kb, doublereal *a, integer *lda, integer *jpvt, doublereal *tau, doublereal *vn1, doublereal *vn2, doublereal *auxv, doublereal *f, integer *ldf)
{
    /* System generated locals */
    integer a_dim1, a_offset, f_dim1, f_offset, i__1, i__2;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal);
    integer i_dnnt(doublereal *);
    /* Local variables */
    integer j, k, rk;
    doublereal akk;
    integer pvt;
    doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    doublereal temp2, tol3z;
    extern /* Subroutine */
    int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *), dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *);
    integer itemp;
    extern /* Subroutine */
    int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
    int dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    integer lsticc, lastrk;
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
    --auxv;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    /* Function Body */
    /* Computing MIN */
    i__1 = *m;
    i__2 = *n + *offset; // , expr subst
    lastrk = min(i__1,i__2);
    lsticc = 0;
    k = 0;
    tol3z = sqrt(dlamch_("Epsilon"));
    /* Beginning of while loop. */
L10:
    if (k < *nb && lsticc == 0)
    {
        ++k;
        rk = *offset + k;
        /* Determine ith pivot column and swap if necessary */
        i__1 = *n - k + 1;
        pvt = k - 1 + idamax_(&i__1, &vn1[k], &c__1);
        if (pvt != k)
        {
            dswap_(m, &a[pvt * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
            i__1 = k - 1;
            dswap_(&i__1, &f[pvt + f_dim1], ldf, &f[k + f_dim1], ldf);
            itemp = jpvt[pvt];
            jpvt[pvt] = jpvt[k];
            jpvt[k] = itemp;
            vn1[pvt] = vn1[k];
            vn2[pvt] = vn2[k];
        }
        /* Apply previous Householder reflectors to column K: */
        /* A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**T. */
        if (k > 1)
        {
            i__1 = *m - rk + 1;
            i__2 = k - 1;
            dgemv_("No transpose", &i__1, &i__2, &c_b8, &a[rk + a_dim1], lda, &f[k + f_dim1], ldf, &c_b9, &a[rk + k * a_dim1], &c__1);
        }
        /* Generate elementary reflector H(k). */
        if (rk < *m)
        {
            i__1 = *m - rk + 1;
            dlarfg_(&i__1, &a[rk + k * a_dim1], &a[rk + 1 + k * a_dim1], & c__1, &tau[k]);
        }
        else
        {
            dlarfg_(&c__1, &a[rk + k * a_dim1], &a[rk + k * a_dim1], &c__1, & tau[k]);
        }
        akk = a[rk + k * a_dim1];
        a[rk + k * a_dim1] = 1.;
        /* Compute Kth column of F: */
        /* Compute F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**T*A(RK:M,K). */
        if (k < *n)
        {
            i__1 = *m - rk + 1;
            i__2 = *n - k;
            dgemv_("Transpose", &i__1, &i__2, &tau[k], &a[rk + (k + 1) * a_dim1], lda, &a[rk + k * a_dim1], &c__1, &c_b16, &f[k + 1 + k * f_dim1], &c__1);
        }
        /* Padding F(1:K,K) with zeros. */
        i__1 = k;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            f[j + k * f_dim1] = 0.;
            /* L20: */
        }
        /* Incremental updating of F: */
        /* F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**T */
        /* *A(RK:M,K). */
        if (k > 1)
        {
            i__1 = *m - rk + 1;
            i__2 = k - 1;
            d__1 = -tau[k];
            dgemv_("Transpose", &i__1, &i__2, &d__1, &a[rk + a_dim1], lda, &a[ rk + k * a_dim1], &c__1, &c_b16, &auxv[1], &c__1);
            i__1 = k - 1;
            dgemv_("No transpose", n, &i__1, &c_b9, &f[f_dim1 + 1], ldf, & auxv[1], &c__1, &c_b9, &f[k * f_dim1 + 1], &c__1);
        }
        /* Update the current row of A: */
        /* A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**T. */
        if (k < *n)
        {
            i__1 = *n - k;
            dgemv_("No transpose", &i__1, &k, &c_b8, &f[k + 1 + f_dim1], ldf, &a[rk + a_dim1], lda, &c_b9, &a[rk + (k + 1) * a_dim1], lda);
        }
        /* Update partial column norms. */
        if (rk < lastrk)
        {
            i__1 = *n;
            for (j = k + 1;
                    j <= i__1;
                    ++j)
            {
                if (vn1[j] != 0.)
                {
                    /* NOTE: The following 4 lines follow from the analysis in */
                    /* Lapack Working Note 176. */
                    temp = (d__1 = a[rk + j * a_dim1], f2c_abs(d__1)) / vn1[j];
                    /* Computing MAX */
                    d__1 = 0.;
                    d__2 = (temp + 1.) * (1. - temp); // , expr subst
                    temp = max(d__1,d__2);
                    /* Computing 2nd power */
                    d__1 = vn1[j] / vn2[j];
                    temp2 = temp * (d__1 * d__1);
                    if (temp2 <= tol3z)
                    {
                        vn2[j] = (doublereal) lsticc;
                        lsticc = j;
                    }
                    else
                    {
                        vn1[j] *= sqrt(temp);
                    }
                }
                /* L30: */
            }
        }
        a[rk + k * a_dim1] = akk;
        /* End of while loop. */
        goto L10;
    }
    *kb = k;
    rk = *offset + *kb;
    /* Apply the block reflector to the rest of the matrix: */
    /* A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) - */
    /* A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**T. */
    /* Computing MIN */
    i__1 = *n;
    i__2 = *m - *offset; // , expr subst
    if (*kb < min(i__1,i__2))
    {
        i__1 = *m - rk;
        i__2 = *n - *kb;
        dgemm_("No transpose", "Transpose", &i__1, &i__2, kb, &c_b8, &a[rk + 1 + a_dim1], lda, &f[*kb + 1 + f_dim1], ldf, &c_b9, &a[rk + 1 + (*kb + 1) * a_dim1], lda);
    }
    /* Recomputation of difficult columns. */
L40:
    if (lsticc > 0)
    {
        itemp = i_dnnt(&vn2[lsticc]);
        i__1 = *m - rk;
        vn1[lsticc] = dnrm2_(&i__1, &a[rk + 1 + lsticc * a_dim1], &c__1);
        /* NOTE: The computation of VN1( LSTICC ) relies on the fact that */
        /* SNRM2 does not fail on vectors with norm below the value of */
        /* SQRT(DLAMCH('S')) */
        vn2[lsticc] = vn1[lsticc];
        lsticc = itemp;
        goto L40;
    }
    return 0;
    /* End of DLAQPS */
}
/* dlaqps_ */
