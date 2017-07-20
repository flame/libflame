/* ../netlib/clasyf_rook.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 =
{
    1.f,0.f
}
;
static integer c__1 = 1;
/* > \brief \b CLASYF_ROOK computes a partial factorization of a complex symmetric matrix using the bounded Bu nch-Kaufman ("rook") diagonal pivoting method. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLASYF_ROOK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clasyf_ rook.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clasyf_ rook.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clasyf_ rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLASYF_ROOK( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, KB, LDA, LDW, N, NB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), W( LDW, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLASYF_ROOK computes a partial factorization of a complex symmetric */
/* > matrix A using the bounded Bunch-Kaufman ("rook") diagonal */
/* > pivoting method. The partial factorization has the form: */
/* > */
/* > A = ( I U12 ) ( A11 0 ) ( I 0 ) if UPLO = 'U', or: */
/* > ( 0 U22 ) ( 0 D ) ( U12**T U22**T ) */
/* > */
/* > A = ( L11 0 ) ( D 0 ) ( L11**T L21**T ) if UPLO = 'L' */
/* > ( L21 I ) ( 0 A22 ) ( 0 I ) */
/* > */
/* > where the order of D is at most NB. The actual order is returned in */
/* > the argument KB, and is either NB or NB-1, or N if N <= NB. */
/* > */
/* > CLASYF_ROOK is an auxiliary routine called by CSYTRF_ROOK. It uses */
/* > blocked code (calling Level 3 BLAS) to update the submatrix */
/* > A11 (if UPLO = 'U') or A22 (if UPLO = 'L'). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > symmetric matrix A is stored: */
/* > = 'U': Upper triangular */
/* > = 'L': Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The maximum number of columns of the matrix A that should be */
/* > factored. NB should be at least 2 to allow for 2-by-2 pivot */
/* > blocks. */
/* > \endverbatim */
/* > */
/* > \param[out] KB */
/* > \verbatim */
/* > KB is INTEGER */
/* > The number of columns of A that were actually factored. */
/* > KB is either NB-1 or NB, or N if N <= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the symmetric matrix A. If UPLO = 'U', the leading */
/* > n-by-n upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading n-by-n lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > On exit, A contains details of the partial factorization. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D. */
/* > */
/* > If UPLO = 'U': */
/* > Only the last KB elements of IPIV are set. */
/* > */
/* > If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* > interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* > If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and */
/* > columns k and -IPIV(k) were interchanged and rows and */
/* > columns k-1 and -IPIV(k-1) were inerchaged, */
/* > D(k-1:k,k-1:k) is a 2-by-2 diagonal block. */
/* > */
/* > If UPLO = 'L': */
/* > Only the first KB elements of IPIV are set. */
/* > */
/* > If IPIV(k) > 0, then rows and columns k and IPIV(k) */
/* > were interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* > If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and */
/* > columns k and -IPIV(k) were interchanged and rows and */
/* > columns k+1 and -IPIV(k+1) were inerchaged, */
/* > D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is COMPLEX array, dimension (LDW,NB) */
/* > \endverbatim */
/* > */
/* > \param[in] LDW */
/* > \verbatim */
/* > LDW is INTEGER */
/* > The leading dimension of the array W. LDW >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > > 0: if INFO = k, D(k,k) is exactly zero. The factorization */
/* > has been completed, but the block diagonal matrix D is */
/* > exactly singular. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup complexSYcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > November 2013, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* > School of Mathematics, */
/* > University of Manchester */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int clasyf_rook_(char *uplo, integer *n, integer *nb, integer *kb, complex *a, integer *lda, integer *ipiv, complex *w, integer *ldw, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4;
    /* Builtin functions */
    double sqrt(doublereal), r_imag(complex *);
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    integer j, k, p;
    complex t, r1, d11, d12, d21, d22;
    integer jb, ii, jj, kk, kp, kw, jp1, jp2, kkw;
    logical done;
    integer imax, jmax;
    real alpha;
    extern /* Subroutine */
    int cscal_(integer *, complex *, complex *, integer *), cgemm_(char *, char *, integer *, integer *, integer * , complex *, complex *, integer *, complex *, integer *, complex * , complex *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int cgemv_(char *, integer *, integer *, complex * , complex *, integer *, complex *, integer *, complex *, complex * , integer *);
    real sfmin;
    extern /* Subroutine */
    int ccopy_(integer *, complex *, integer *, complex *, integer *);
    integer itemp;
    extern /* Subroutine */
    int cswap_(integer *, complex *, integer *, complex *, integer *);
    integer kstep;
    real stemp, absakk;
    extern integer icamax_(integer *, complex *, integer *);
    extern real slamch_(char *);
    real colmax, rowmax;
    /* -- LAPACK computational routine (version 3.5.0) -- */
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    w_dim1 = *ldw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    /* Function Body */
    *info = 0;
    /* Initialize ALPHA for use in choosing pivot block size. */
    alpha = (sqrt(17.f) + 1.f) / 8.f;
    /* Compute machine safe minimum */
    sfmin = slamch_("S");
    if (lsame_(uplo, "U"))
    {
        /* Factorize the trailing columns of A using the upper triangle */
        /* of A and working backwards, and compute the matrix W = U12*D */
        /* for use in updating A11 */
        /* K is the main loop index, decreasing from N in steps of 1 or 2 */
        k = *n;
L10: /* KW is the column of W which corresponds to column K of A */
        kw = *nb + k - *n;
        /* Exit from loop */
        if (k <= *n - *nb + 1 && *nb < *n || k < 1)
        {
            goto L30;
        }
        kstep = 1;
        p = k;
        /* Copy column K of A to column KW of W and update it */
        ccopy_(&k, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
        if (k < *n)
        {
            i__1 = *n - k;
            q__1.r = -1.f;
            q__1.i = -0.f; // , expr subst
            cgemv_("No transpose", &k, &i__1, &q__1, &a[(k + 1) * a_dim1 + 1], lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b1, &w[kw * w_dim1 + 1], &c__1);
        }
        /* Determine rows and columns to be interchanged and whether */
        /* a 1-by-1 or 2-by-2 pivot block will be used */
        i__1 = k + kw * w_dim1;
        absakk = (r__1 = w[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&w[k + kw * w_dim1]), f2c_abs(r__2));
        /* IMAX is the row-index of the largest off-diagonal element in */
        /* column K, and COLMAX is its absolute value. */
        /* Determine both COLMAX and IMAX. */
        if (k > 1)
        {
            i__1 = k - 1;
            imax = icamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
            i__1 = imax + kw * w_dim1;
            colmax = (r__1 = w[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&w[imax + kw * w_dim1]), f2c_abs(r__2));
        }
        else
        {
            colmax = 0.f;
        }
        if (max(absakk,colmax) == 0.f)
        {
            /* Column K is zero or underflow: set INFO and continue */
            if (*info == 0)
            {
                *info = k;
            }
            kp = k;
            ccopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
        }
        else
        {
            /* ============================================================ */
            /* Test for interchange */
            /* Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
            /* (used to handle NaN and Inf) */
            if (! (absakk < alpha * colmax))
            {
                /* no interchange, use 1-by-1 pivot block */
                kp = k;
            }
            else
            {
                done = FALSE_;
                /* Loop until pivot found */
L12: /* Begin pivot search loop body */
                /* Copy column IMAX to column KW-1 of W and update it */
                ccopy_(&imax, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
                i__1 = k - imax;
                ccopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 1 + (kw - 1) * w_dim1], &c__1);
                if (k < *n)
                {
                    i__1 = *n - k;
                    q__1.r = -1.f;
                    q__1.i = -0.f; // , expr subst
                    cgemv_("No transpose", &k, &i__1, &q__1, &a[(k + 1) * a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1], &c__1);
                }
                /* JMAX is the column-index of the largest off-diagonal */
                /* element in row IMAX, and ROWMAX is its absolute value. */
                /* Determine both ROWMAX and JMAX. */
                if (imax != k)
                {
                    i__1 = k - imax;
                    jmax = imax + icamax_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1], &c__1);
                    i__1 = jmax + (kw - 1) * w_dim1;
                    rowmax = (r__1 = w[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(& w[jmax + (kw - 1) * w_dim1]), f2c_abs(r__2));
                }
                else
                {
                    rowmax = 0.f;
                }
                if (imax > 1)
                {
                    i__1 = imax - 1;
                    itemp = icamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
                    i__1 = itemp + (kw - 1) * w_dim1;
                    stemp = (r__1 = w[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&w[ itemp + (kw - 1) * w_dim1]), f2c_abs(r__2));
                    if (stemp > rowmax)
                    {
                        rowmax = stemp;
                        jmax = itemp;
                    }
                }
                /* Equivalent to testing for */
                /* CABS1( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX */
                /* (used to handle NaN and Inf) */
                i__1 = imax + (kw - 1) * w_dim1;
                if (! ((r__1 = w[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&w[imax + (kw - 1) * w_dim1]), f2c_abs(r__2)) < alpha * rowmax))
                {
                    /* interchange rows and columns K and IMAX, */
                    /* use 1-by-1 pivot block */
                    kp = imax;
                    /* copy column KW-1 of W to column KW of W */
                    ccopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
                    done = TRUE_;
                    /* Equivalent to testing for ROWMAX.EQ.COLMAX, */
                    /* (used to handle NaN and Inf) */
                }
                else if (p == jmax || rowmax <= colmax)
                {
                    /* interchange rows and columns K-1 and IMAX, */
                    /* use 2-by-2 pivot block */
                    kp = imax;
                    kstep = 2;
                    done = TRUE_;
                }
                else
                {
                    /* Pivot not found: set params and repeat */
                    p = imax;
                    colmax = rowmax;
                    imax = jmax;
                    /* Copy updated JMAXth (next IMAXth) column to Kth of W */
                    ccopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
                }
                /* End pivot search loop body */
                if (! done)
                {
                    goto L12;
                }
            }
            /* ============================================================ */
            kk = k - kstep + 1;
            /* KKW is the column of W which corresponds to column KK of A */
            kkw = *nb + kk - *n;
            if (kstep == 2 && p != k)
            {
                /* Copy non-updated column K to column P */
                i__1 = k - p;
                ccopy_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 1) * a_dim1], lda);
                ccopy_(&p, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], & c__1);
                /* Interchange rows K and P in last N-K+1 columns of A */
                /* and last N-K+2 columns of W */
                i__1 = *n - k + 1;
                cswap_(&i__1, &a[k + k * a_dim1], lda, &a[p + k * a_dim1], lda);
                i__1 = *n - kk + 1;
                cswap_(&i__1, &w[k + kkw * w_dim1], ldw, &w[p + kkw * w_dim1], ldw);
            }
            /* Updated column KP is already stored in column KKW of W */
            if (kp != kk)
            {
                /* Copy non-updated column KK to column KP */
                i__1 = kp + k * a_dim1;
                i__2 = kk + k * a_dim1;
                a[i__1].r = a[i__2].r;
                a[i__1].i = a[i__2].i; // , expr subst
                i__1 = k - 1 - kp;
                ccopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1], lda);
                ccopy_(&kp, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], & c__1);
                /* Interchange rows KK and KP in last N-KK+1 columns */
                /* of A and W */
                i__1 = *n - kk + 1;
                cswap_(&i__1, &a[kk + kk * a_dim1], lda, &a[kp + kk * a_dim1], lda);
                i__1 = *n - kk + 1;
                cswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * w_dim1], ldw);
            }
            if (kstep == 1)
            {
                /* 1-by-1 pivot block D(k): column KW of W now holds */
                /* W(k) = U(k)*D(k) */
                /* where U(k) is the k-th column of U */
                /* Store U(k) in column k of A */
                ccopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], & c__1);
                if (k > 1)
                {
                    i__1 = k + k * a_dim1;
                    if ((r__1 = a[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]), f2c_abs(r__2)) >= sfmin)
                    {
                        c_div(&q__1, &c_b1, &a[k + k * a_dim1]);
                        r1.r = q__1.r;
                        r1.i = q__1.i; // , expr subst
                        i__1 = k - 1;
                        cscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
                    }
                    else /* if(complicated condition) */
                    {
                        i__1 = k + k * a_dim1;
                        if (a[i__1].r != 0.f || a[i__1].i != 0.f)
                        {
                            i__1 = k - 1;
                            for (ii = 1;
                                    ii <= i__1;
                                    ++ii)
                            {
                                i__2 = ii + k * a_dim1;
                                c_div(&q__1, &a[ii + k * a_dim1], &a[k + k * a_dim1]);
                                a[i__2].r = q__1.r;
                                a[i__2].i = q__1.i; // , expr subst
                                /* L14: */
                            }
                        }
                    }
                }
            }
            else
            {
                /* 2-by-2 pivot block D(k): columns KW and KW-1 of W now */
                /* hold */
                /* ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */
                /* where U(k) and U(k-1) are the k-th and (k-1)-th columns */
                /* of U */
                if (k > 2)
                {
                    /* Store U(k) and U(k-1) in columns k and k-1 of A */
                    i__1 = k - 1 + kw * w_dim1;
                    d12.r = w[i__1].r;
                    d12.i = w[i__1].i; // , expr subst
                    c_div(&q__1, &w[k + kw * w_dim1], &d12);
                    d11.r = q__1.r;
                    d11.i = q__1.i; // , expr subst
                    c_div(&q__1, &w[k - 1 + (kw - 1) * w_dim1], &d12);
                    d22.r = q__1.r;
                    d22.i = q__1.i; // , expr subst
                    q__3.r = d11.r * d22.r - d11.i * d22.i;
                    q__3.i = d11.r * d22.i + d11.i * d22.r; // , expr subst
                    q__2.r = q__3.r - 1.f;
                    q__2.i = q__3.i - 0.f; // , expr subst
                    c_div(&q__1, &c_b1, &q__2);
                    t.r = q__1.r;
                    t.i = q__1.i; // , expr subst
                    i__1 = k - 2;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = j + (k - 1) * a_dim1;
                        i__3 = j + (kw - 1) * w_dim1;
                        q__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i;
                        q__4.i = d11.r * w[i__3].i + d11.i * w[i__3] .r; // , expr subst
                        i__4 = j + kw * w_dim1;
                        q__3.r = q__4.r - w[i__4].r;
                        q__3.i = q__4.i - w[i__4] .i; // , expr subst
                        c_div(&q__2, &q__3, &d12);
                        q__1.r = t.r * q__2.r - t.i * q__2.i;
                        q__1.i = t.r * q__2.i + t.i * q__2.r; // , expr subst
                        a[i__2].r = q__1.r;
                        a[i__2].i = q__1.i; // , expr subst
                        i__2 = j + k * a_dim1;
                        i__3 = j + kw * w_dim1;
                        q__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i;
                        q__4.i = d22.r * w[i__3].i + d22.i * w[i__3] .r; // , expr subst
                        i__4 = j + (kw - 1) * w_dim1;
                        q__3.r = q__4.r - w[i__4].r;
                        q__3.i = q__4.i - w[i__4] .i; // , expr subst
                        c_div(&q__2, &q__3, &d12);
                        q__1.r = t.r * q__2.r - t.i * q__2.i;
                        q__1.i = t.r * q__2.i + t.i * q__2.r; // , expr subst
                        a[i__2].r = q__1.r;
                        a[i__2].i = q__1.i; // , expr subst
                        /* L20: */
                    }
                }
                /* Copy D(k) to A */
                i__1 = k - 1 + (k - 1) * a_dim1;
                i__2 = k - 1 + (kw - 1) * w_dim1;
                a[i__1].r = w[i__2].r;
                a[i__1].i = w[i__2].i; // , expr subst
                i__1 = k - 1 + k * a_dim1;
                i__2 = k - 1 + kw * w_dim1;
                a[i__1].r = w[i__2].r;
                a[i__1].i = w[i__2].i; // , expr subst
                i__1 = k + k * a_dim1;
                i__2 = k + kw * w_dim1;
                a[i__1].r = w[i__2].r;
                a[i__1].i = w[i__2].i; // , expr subst
            }
        }
        /* Store details of the interchanges in IPIV */
        if (kstep == 1)
        {
            ipiv[k] = kp;
        }
        else
        {
            ipiv[k] = -p;
            ipiv[k - 1] = -kp;
        }
        /* Decrease K and return to the start of the main loop */
        k -= kstep;
        goto L10;
L30: /* Update the upper triangle of A11 (= A(1:k,1:k)) as */
        /* A11 := A11 - U12*D*U12**T = A11 - U12*W**T */
        /* computing blocks of NB columns at a time */
        i__1 = -(*nb);
        for (j = (k - 1) / *nb * *nb + 1;
                i__1 < 0 ? j >= 1 : j <= 1;
                j += i__1)
        {
            /* Computing MIN */
            i__2 = *nb;
            i__3 = k - j + 1; // , expr subst
            jb = min(i__2,i__3);
            /* Update the upper triangle of the diagonal block */
            i__2 = j + jb - 1;
            for (jj = j;
                    jj <= i__2;
                    ++jj)
            {
                i__3 = jj - j + 1;
                i__4 = *n - k;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgemv_("No transpose", &i__3, &i__4, &q__1, &a[j + (k + 1) * a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b1, &a[j + jj * a_dim1], &c__1);
                /* L40: */
            }
            /* Update the rectangular superdiagonal block */
            if (j >= 2)
            {
                i__2 = j - 1;
                i__3 = *n - k;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &q__1, &a[(k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * w_dim1], ldw, &c_b1, &a[j * a_dim1 + 1], lda);
            }
            /* L50: */
        }
        /* Put U12 in standard form by partially undoing the interchanges */
        /* in columns k+1:n */
        j = k + 1;
L60:
        kstep = 1;
        jp1 = 1;
        jj = j;
        jp2 = ipiv[j];
        if (jp2 < 0)
        {
            jp2 = -jp2;
            ++j;
            jp1 = -ipiv[j];
            kstep = 2;
        }
        ++j;
        if (jp2 != jj && j <= *n)
        {
            i__1 = *n - j + 1;
            cswap_(&i__1, &a[jp2 + j * a_dim1], lda, &a[jj + j * a_dim1], lda) ;
        }
        jj = j - 1;
        if (jp1 != jj && kstep == 2)
        {
            i__1 = *n - j + 1;
            cswap_(&i__1, &a[jp1 + j * a_dim1], lda, &a[jj + j * a_dim1], lda) ;
        }
        if (j <= *n)
        {
            goto L60;
        }
        /* Set KB to the number of columns factorized */
        *kb = *n - k;
    }
    else
    {
        /* Factorize the leading columns of A using the lower triangle */
        /* of A and working forwards, and compute the matrix W = L21*D */
        /* for use in updating A22 */
        /* K is the main loop index, increasing from 1 in steps of 1 or 2 */
        k = 1;
L70: /* Exit from loop */
        if (k >= *nb && *nb < *n || k > *n)
        {
            goto L90;
        }
        kstep = 1;
        p = k;
        /* Copy column K of A to column K of W and update it */
        i__1 = *n - k + 1;
        ccopy_(&i__1, &a[k + k * a_dim1], &c__1, &w[k + k * w_dim1], &c__1);
        if (k > 1)
        {
            i__1 = *n - k + 1;
            i__2 = k - 1;
            q__1.r = -1.f;
            q__1.i = -0.f; // , expr subst
            cgemv_("No transpose", &i__1, &i__2, &q__1, &a[k + a_dim1], lda, & w[k + w_dim1], ldw, &c_b1, &w[k + k * w_dim1], &c__1);
        }
        /* Determine rows and columns to be interchanged and whether */
        /* a 1-by-1 or 2-by-2 pivot block will be used */
        i__1 = k + k * w_dim1;
        absakk = (r__1 = w[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&w[k + k * w_dim1]), f2c_abs(r__2));
        /* IMAX is the row-index of the largest off-diagonal element in */
        /* column K, and COLMAX is its absolute value. */
        /* Determine both COLMAX and IMAX. */
        if (k < *n)
        {
            i__1 = *n - k;
            imax = k + icamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
            i__1 = imax + k * w_dim1;
            colmax = (r__1 = w[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&w[imax + k * w_dim1]), f2c_abs(r__2));
        }
        else
        {
            colmax = 0.f;
        }
        if (max(absakk,colmax) == 0.f)
        {
            /* Column K is zero or underflow: set INFO and continue */
            if (*info == 0)
            {
                *info = k;
            }
            kp = k;
            i__1 = *n - k + 1;
            ccopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], & c__1);
        }
        else
        {
            /* ============================================================ */
            /* Test for interchange */
            /* Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
            /* (used to handle NaN and Inf) */
            if (! (absakk < alpha * colmax))
            {
                /* no interchange, use 1-by-1 pivot block */
                kp = k;
            }
            else
            {
                done = FALSE_;
                /* Loop until pivot found */
L72: /* Begin pivot search loop body */
                /* Copy column IMAX to column K+1 of W and update it */
                i__1 = imax - k;
                ccopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * w_dim1], &c__1);
                i__1 = *n - imax + 1;
                ccopy_(&i__1, &a[imax + imax * a_dim1], &c__1, &w[imax + (k + 1) * w_dim1], &c__1);
                if (k > 1)
                {
                    i__1 = *n - k + 1;
                    i__2 = k - 1;
                    q__1.r = -1.f;
                    q__1.i = -0.f; // , expr subst
                    cgemv_("No transpose", &i__1, &i__2, &q__1, &a[k + a_dim1] , lda, &w[imax + w_dim1], ldw, &c_b1, &w[k + (k + 1) * w_dim1], &c__1);
                }
                /* JMAX is the column-index of the largest off-diagonal */
                /* element in row IMAX, and ROWMAX is its absolute value. */
                /* Determine both ROWMAX and JMAX. */
                if (imax != k)
                {
                    i__1 = imax - k;
                    jmax = k - 1 + icamax_(&i__1, &w[k + (k + 1) * w_dim1], & c__1);
                    i__1 = jmax + (k + 1) * w_dim1;
                    rowmax = (r__1 = w[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(& w[jmax + (k + 1) * w_dim1]), f2c_abs(r__2));
                }
                else
                {
                    rowmax = 0.f;
                }
                if (imax < *n)
                {
                    i__1 = *n - imax;
                    itemp = imax + icamax_(&i__1, &w[imax + 1 + (k + 1) * w_dim1], &c__1);
                    i__1 = itemp + (k + 1) * w_dim1;
                    stemp = (r__1 = w[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&w[ itemp + (k + 1) * w_dim1]), f2c_abs(r__2));
                    if (stemp > rowmax)
                    {
                        rowmax = stemp;
                        jmax = itemp;
                    }
                }
                /* Equivalent to testing for */
                /* CABS1( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX */
                /* (used to handle NaN and Inf) */
                i__1 = imax + (k + 1) * w_dim1;
                if (! ((r__1 = w[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&w[imax + (k + 1) * w_dim1]), f2c_abs(r__2)) < alpha * rowmax))
                {
                    /* interchange rows and columns K and IMAX, */
                    /* use 1-by-1 pivot block */
                    kp = imax;
                    /* copy column K+1 of W to column K of W */
                    i__1 = *n - k + 1;
                    ccopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * w_dim1], &c__1);
                    done = TRUE_;
                    /* Equivalent to testing for ROWMAX.EQ.COLMAX, */
                    /* (used to handle NaN and Inf) */
                }
                else if (p == jmax || rowmax <= colmax)
                {
                    /* interchange rows and columns K+1 and IMAX, */
                    /* use 2-by-2 pivot block */
                    kp = imax;
                    kstep = 2;
                    done = TRUE_;
                }
                else
                {
                    /* Pivot not found: set params and repeat */
                    p = imax;
                    colmax = rowmax;
                    imax = jmax;
                    /* Copy updated JMAXth (next IMAXth) column to Kth of W */
                    i__1 = *n - k + 1;
                    ccopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * w_dim1], &c__1);
                }
                /* End pivot search loop body */
                if (! done)
                {
                    goto L72;
                }
            }
            /* ============================================================ */
            kk = k + kstep - 1;
            if (kstep == 2 && p != k)
            {
                /* Copy non-updated column K to column P */
                i__1 = p - k;
                ccopy_(&i__1, &a[k + k * a_dim1], &c__1, &a[p + k * a_dim1], lda);
                i__1 = *n - p + 1;
                ccopy_(&i__1, &a[p + k * a_dim1], &c__1, &a[p + p * a_dim1], & c__1);
                /* Interchange rows K and P in first K columns of A */
                /* and first K+1 columns of W */
                cswap_(&k, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
                cswap_(&kk, &w[k + w_dim1], ldw, &w[p + w_dim1], ldw);
            }
            /* Updated column KP is already stored in column KK of W */
            if (kp != kk)
            {
                /* Copy non-updated column KK to column KP */
                i__1 = kp + k * a_dim1;
                i__2 = kk + k * a_dim1;
                a[i__1].r = a[i__2].r;
                a[i__1].i = a[i__2].i; // , expr subst
                i__1 = kp - k - 1;
                ccopy_(&i__1, &a[k + 1 + kk * a_dim1], &c__1, &a[kp + (k + 1) * a_dim1], lda);
                i__1 = *n - kp + 1;
                ccopy_(&i__1, &a[kp + kk * a_dim1], &c__1, &a[kp + kp * a_dim1], &c__1);
                /* Interchange rows KK and KP in first KK columns of A and W */
                cswap_(&kk, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
                cswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
            }
            if (kstep == 1)
            {
                /* 1-by-1 pivot block D(k): column k of W now holds */
                /* W(k) = L(k)*D(k) */
                /* where L(k) is the k-th column of L */
                /* Store L(k) in column k of A */
                i__1 = *n - k + 1;
                ccopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], & c__1);
                if (k < *n)
                {
                    i__1 = k + k * a_dim1;
                    if ((r__1 = a[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&a[k + k * a_dim1]), f2c_abs(r__2)) >= sfmin)
                    {
                        c_div(&q__1, &c_b1, &a[k + k * a_dim1]);
                        r1.r = q__1.r;
                        r1.i = q__1.i; // , expr subst
                        i__1 = *n - k;
                        cscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
                    }
                    else /* if(complicated condition) */
                    {
                        i__1 = k + k * a_dim1;
                        if (a[i__1].r != 0.f || a[i__1].i != 0.f)
                        {
                            i__1 = *n;
                            for (ii = k + 1;
                                    ii <= i__1;
                                    ++ii)
                            {
                                i__2 = ii + k * a_dim1;
                                c_div(&q__1, &a[ii + k * a_dim1], &a[k + k * a_dim1]);
                                a[i__2].r = q__1.r;
                                a[i__2].i = q__1.i; // , expr subst
                                /* L74: */
                            }
                        }
                    }
                }
            }
            else
            {
                /* 2-by-2 pivot block D(k): columns k and k+1 of W now hold */
                /* ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */
                /* where L(k) and L(k+1) are the k-th and (k+1)-th columns */
                /* of L */
                if (k < *n - 1)
                {
                    /* Store L(k) and L(k+1) in columns k and k+1 of A */
                    i__1 = k + 1 + k * w_dim1;
                    d21.r = w[i__1].r;
                    d21.i = w[i__1].i; // , expr subst
                    c_div(&q__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
                    d11.r = q__1.r;
                    d11.i = q__1.i; // , expr subst
                    c_div(&q__1, &w[k + k * w_dim1], &d21);
                    d22.r = q__1.r;
                    d22.i = q__1.i; // , expr subst
                    q__3.r = d11.r * d22.r - d11.i * d22.i;
                    q__3.i = d11.r * d22.i + d11.i * d22.r; // , expr subst
                    q__2.r = q__3.r - 1.f;
                    q__2.i = q__3.i - 0.f; // , expr subst
                    c_div(&q__1, &c_b1, &q__2);
                    t.r = q__1.r;
                    t.i = q__1.i; // , expr subst
                    i__1 = *n;
                    for (j = k + 2;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = j + k * a_dim1;
                        i__3 = j + k * w_dim1;
                        q__4.r = d11.r * w[i__3].r - d11.i * w[i__3].i;
                        q__4.i = d11.r * w[i__3].i + d11.i * w[i__3] .r; // , expr subst
                        i__4 = j + (k + 1) * w_dim1;
                        q__3.r = q__4.r - w[i__4].r;
                        q__3.i = q__4.i - w[i__4] .i; // , expr subst
                        c_div(&q__2, &q__3, &d21);
                        q__1.r = t.r * q__2.r - t.i * q__2.i;
                        q__1.i = t.r * q__2.i + t.i * q__2.r; // , expr subst
                        a[i__2].r = q__1.r;
                        a[i__2].i = q__1.i; // , expr subst
                        i__2 = j + (k + 1) * a_dim1;
                        i__3 = j + (k + 1) * w_dim1;
                        q__4.r = d22.r * w[i__3].r - d22.i * w[i__3].i;
                        q__4.i = d22.r * w[i__3].i + d22.i * w[i__3] .r; // , expr subst
                        i__4 = j + k * w_dim1;
                        q__3.r = q__4.r - w[i__4].r;
                        q__3.i = q__4.i - w[i__4] .i; // , expr subst
                        c_div(&q__2, &q__3, &d21);
                        q__1.r = t.r * q__2.r - t.i * q__2.i;
                        q__1.i = t.r * q__2.i + t.i * q__2.r; // , expr subst
                        a[i__2].r = q__1.r;
                        a[i__2].i = q__1.i; // , expr subst
                        /* L80: */
                    }
                }
                /* Copy D(k) to A */
                i__1 = k + k * a_dim1;
                i__2 = k + k * w_dim1;
                a[i__1].r = w[i__2].r;
                a[i__1].i = w[i__2].i; // , expr subst
                i__1 = k + 1 + k * a_dim1;
                i__2 = k + 1 + k * w_dim1;
                a[i__1].r = w[i__2].r;
                a[i__1].i = w[i__2].i; // , expr subst
                i__1 = k + 1 + (k + 1) * a_dim1;
                i__2 = k + 1 + (k + 1) * w_dim1;
                a[i__1].r = w[i__2].r;
                a[i__1].i = w[i__2].i; // , expr subst
            }
        }
        /* Store details of the interchanges in IPIV */
        if (kstep == 1)
        {
            ipiv[k] = kp;
        }
        else
        {
            ipiv[k] = -p;
            ipiv[k + 1] = -kp;
        }
        /* Increase K and return to the start of the main loop */
        k += kstep;
        goto L70;
L90: /* Update the lower triangle of A22 (= A(k:n,k:n)) as */
        /* A22 := A22 - L21*D*L21**T = A22 - L21*W**T */
        /* computing blocks of NB columns at a time */
        i__1 = *n;
        i__2 = *nb;
        for (j = k;
                i__2 < 0 ? j >= i__1 : j <= i__1;
                j += i__2)
        {
            /* Computing MIN */
            i__3 = *nb;
            i__4 = *n - j + 1; // , expr subst
            jb = min(i__3,i__4);
            /* Update the lower triangle of the diagonal block */
            i__3 = j + jb - 1;
            for (jj = j;
                    jj <= i__3;
                    ++jj)
            {
                i__4 = j + jb - jj;
                i__5 = k - 1;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgemv_("No transpose", &i__4, &i__5, &q__1, &a[jj + a_dim1], lda, &w[jj + w_dim1], ldw, &c_b1, &a[jj + jj * a_dim1] , &c__1);
                /* L100: */
            }
            /* Update the rectangular subdiagonal block */
            if (j + jb <= *n)
            {
                i__3 = *n - j - jb + 1;
                i__4 = k - 1;
                q__1.r = -1.f;
                q__1.i = -0.f; // , expr subst
                cgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &q__1, &a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b1, &a[j + jb + j * a_dim1], lda);
            }
            /* L110: */
        }
        /* Put L21 in standard form by partially undoing the interchanges */
        /* in columns 1:k-1 */
        j = k - 1;
L120:
        kstep = 1;
        jp1 = 1;
        jj = j;
        jp2 = ipiv[j];
        if (jp2 < 0)
        {
            jp2 = -jp2;
            --j;
            jp1 = -ipiv[j];
            kstep = 2;
        }
        --j;
        if (jp2 != jj && j >= 1)
        {
            cswap_(&j, &a[jp2 + a_dim1], lda, &a[jj + a_dim1], lda);
        }
        jj = j + 1;
        if (jp1 != jj && kstep == 2)
        {
            cswap_(&j, &a[jp1 + a_dim1], lda, &a[jj + a_dim1], lda);
        }
        if (j >= 1)
        {
            goto L120;
        }
        /* Set KB to the number of columns factorized */
        *kb = k - 1;
    }
    return 0;
    /* End of CLASYF_ROOK */
}
/* clasyf_rook__ */
