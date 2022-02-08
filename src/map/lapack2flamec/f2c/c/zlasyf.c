/* ../netlib/zlasyf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 =
{
    1.,0.
}
;
static integer c__1 = 1;
/* > \brief \b ZLASYF computes a partial factorization of a complex symmetric matrix using the Bunch-Kaufman d iagonal pivoting method. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLASYF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlasyf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlasyf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlasyf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, KB, LDA, LDW, N, NB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), W( LDW, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLASYF computes a partial factorization of a complex symmetric matrix */
/* > A using the Bunch-Kaufman diagonal pivoting method. The partial */
/* > factorization has the form: */
/* > */
/* > A = ( I U12 ) ( A11 0 ) ( I 0 ) if UPLO = 'U', or: */
/* > ( 0 U22 ) ( 0 D ) ( U12**T U22**T ) */
/* > */
/* > A = ( L11 0 ) ( D 0 ) ( L11**T L21**T ) if UPLO = 'L' */
/* > ( L21 I ) ( 0 A22 ) ( 0 I ) */
/* > */
/* > where the order of D is at most NB. The actual order is returned in */
/* > the argument KB, and is either NB or NB-1, or N if N <= NB. */
/* > Note that U**T denotes the transpose of U. */
/* > */
/* > ZLASYF is an auxiliary routine called by ZSYTRF. It uses blocked code */
/* > (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or */
/* > A22 (if UPLO = 'L'). */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
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
/* > If IPIV(k) = IPIV(k-1) < 0, then rows and columns */
/* > k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
/* > is a 2-by-2 diagonal block. */
/* > */
/* > If UPLO = 'L': */
/* > Only the first KB elements of IPIV are set. */
/* > */
/* > If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* > interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* > If IPIV(k) = IPIV(k+1) < 0, then rows and columns */
/* > k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1) */
/* > is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* > W is COMPLEX*16 array, dimension (LDW,NB) */
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
/* > \ingroup complex16SYcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > November 2013, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int zlasyf_(char *uplo, integer *n, integer *nb, integer *kb, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w, integer *ldw, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    snprintf(buffer, 256,"zlasyf inputs: uplo %c, n %d, nb %d, lda %d, ldw %d",*uplo, *n, *nb, *lda, *ldw);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;
    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    /* Local variables */
    integer j, k;
    doublecomplex t, r1, d11, d21, d22;
    integer jb, jj, kk, jp, kp, kw, kkw, imax, jmax;
    doublereal alpha;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int zscal_(integer *, doublecomplex *, doublecomplex *, integer *), zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *);
    integer kstep;
    extern /* Subroutine */
    int zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *), zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *), zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    doublereal absakk, colmax;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    doublereal rowmax;
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
    alpha = (sqrt(17.) + 1.) / 8.;
    if (lsame_(uplo, "U"))
    {
        /* Factorize the trailing columns of A using the upper triangle */
        /* of A and working backwards, and compute the matrix W = U12*D */
        /* for use in updating A11 */
        /* K is the main loop index, decreasing from N in steps of 1 or 2 */
        /* KW is the column of W which corresponds to column K of A */
        k = *n;
L10:
        kw = *nb + k - *n;
        /* Exit from loop */
        if (k <= *n - *nb + 1 && *nb < *n || k < 1)
        {
            goto L30;
        }
        /* Copy column K of A to column KW of W and update it */
        zcopy_(&k, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
        if (k < *n)
        {
            i__1 = *n - k;
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1], lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b1, &w[kw * w_dim1 + 1], &c__1);
        }
        kstep = 1;
        /* Determine rows and columns to be interchanged and whether */
        /* a 1-by-1 or 2-by-2 pivot block will be used */
        i__1 = k + kw * w_dim1;
        absakk = (d__1 = w[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&w[k + kw * w_dim1]), f2c_dabs(d__2));
        /* IMAX is the row-index of the largest off-diagonal element in */
        if (k > 1)
        {
            i__1 = k - 1;
            imax = izamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
            i__1 = imax + kw * w_dim1;
            colmax = (d__1 = w[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&w[imax + kw * w_dim1]), f2c_dabs(d__2));
        }
        else
        {
            colmax = 0.;
        }
        if (max(absakk,colmax) == 0.)
        {
            /* Column K is zero or underflow: set INFO and continue */
            if (*info == 0)
            {
                *info = k;
            }
            kp = k;
        }
        else
        {
            if (absakk >= alpha * colmax)
            {
                /* no interchange, use 1-by-1 pivot block */
                kp = k;
            }
            else
            {
                /* Copy column IMAX to column KW-1 of W and update it */
                zcopy_(&imax, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
                i__1 = k - imax;
                zcopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 1 + (kw - 1) * w_dim1], &c__1);
                if (k < *n)
                {
                    i__1 = *n - k;
                    z__1.r = -1.;
                    z__1.i = -0.; // , expr subst
                    zgemv_("No transpose", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1], &c__1);
                }
                /* JMAX is the column-index of the largest off-diagonal */
                /* element in row IMAX, and ROWMAX is its absolute value */
                i__1 = k - imax;
                jmax = imax + izamax_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1], &c__1);
                i__1 = jmax + (kw - 1) * w_dim1;
                rowmax = (d__1 = w[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&w[ jmax + (kw - 1) * w_dim1]), f2c_dabs(d__2));
                if (imax > 1)
                {
                    i__1 = imax - 1;
                    jmax = izamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
                    /* Computing MAX */
                    i__1 = jmax + (kw - 1) * w_dim1;
                    d__3 = rowmax;
                    d__4 = (d__1 = w[i__1].r, f2c_dabs(d__1)) + ( d__2 = d_imag(&w[jmax + (kw - 1) * w_dim1]), f2c_dabs( d__2)); // , expr subst
                    rowmax = max(d__3,d__4);
                }
                if (absakk >= alpha * colmax * (colmax / rowmax))
                {
                    /* no interchange, use 1-by-1 pivot block */
                    kp = k;
                }
                else /* if(complicated condition) */
                {
                    i__1 = imax + (kw - 1) * w_dim1;
                    if ((d__1 = w[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&w[ imax + (kw - 1) * w_dim1]), f2c_dabs(d__2)) >= alpha * rowmax)
                    {
                        /* interchange rows and columns K and IMAX, use 1-by-1 */
                        /* pivot block */
                        kp = imax;
                        /* copy column KW-1 of W to column KW of W */
                        zcopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
                    }
                    else
                    {
                        /* interchange rows and columns K-1 and IMAX, use 2-by-2 */
                        /* pivot block */
                        kp = imax;
                        kstep = 2;
                    }
                }
            }
            /* ============================================================ */
            /* KK is the column of A where pivoting step stopped */
            kk = k - kstep + 1;
            /* KKW is the column of W which corresponds to column KK of A */
            kkw = *nb + kk - *n;
            /* Interchange rows and columns KP and KK. */
            /* Updated column KP is already stored in column KKW of W. */
            if (kp != kk)
            {
                /* Copy non-updated column KK to column KP of submatrix A */
                /* at step K. No need to copy element into column K */
                /* (or K and K-1 for 2-by-2 pivot) of A, since these columns */
                /* will be later overwritten. */
                i__1 = kp + kp * a_dim1;
                i__2 = kk + kk * a_dim1;
                a[i__1].r = a[i__2].r;
                a[i__1].i = a[i__2].i; // , expr subst
                i__1 = kk - 1 - kp;
                zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1], lda);
                if (kp > 1)
                {
                    i__1 = kp - 1;
                    zcopy_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
                }
                /* Interchange rows KK and KP in last K+1 to N columns of A */
                /* (columns K (or K and K-1 for 2-by-2 pivot) of A will be */
                /* later overwritten). Interchange rows KK and KP */
                /* in last KKW to NB columns of W. */
                if (k < *n)
                {
                    i__1 = *n - k;
                    zswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k + 1) * a_dim1], lda);
                }
                i__1 = *n - kk + 1;
                zswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * w_dim1], ldw);
            }
            if (kstep == 1)
            {
                /* 1-by-1 pivot block D(k): column kw of W now holds */
                /* W(kw) = U(k)*D(k), */
                /* where U(k) is the k-th column of U */
                /* Store subdiag. elements of column U(k) */
                /* and 1-by-1 block D(k) in column k of A. */
                /* NOTE: Diagonal element U(k,k) is a UNIT element */
                /* and not stored. */
                /* A(k,k) := D(k,k) = W(k,kw) */
                /* A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k) */
                zcopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], & c__1);
                z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
                r1.r = z__1.r;
                r1.i = z__1.i; // , expr subst
                i__1 = k - 1;
                zscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
            }
            else
            {
                /* 2-by-2 pivot block D(k): columns kw and kw-1 of W now hold */
                /* ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k) */
                /* where U(k) and U(k-1) are the k-th and (k-1)-th columns */
                /* of U */
                /* Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2 */
                /* block D(k-1:k,k-1:k) in columns k-1 and k of A. */
                /* NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT */
                /* block and not stored. */
                /* A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw) */
                /* A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) = */
                /* = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) ) */
                if (k > 2)
                {
                    /* Compose the columns of the inverse of 2-by-2 pivot */
                    /* block D in the following way to reduce the number */
                    /* of FLOPS when we myltiply panel ( W(kw-1) W(kw) ) by */
                    /* this inverse */
                    /* D**(-1) = ( d11 d21 )**(-1) = */
                    /* ( d21 d22 ) */
                    /* = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) = */
                    /* ( (-d21 ) ( d11 ) ) */
                    /* = 1/d21 * 1/((d11/d21)*(d22/d21)-1) * */
                    /* * ( ( d22/d21 ) ( -1 ) ) = */
                    /* ( ( -1 ) ( d11/d21 ) ) */
                    /* = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) ( -1 ) ) = */
                    /* ( ( -1 ) ( D22 ) ) */
                    /* = 1/d21 * T * ( ( D11 ) ( -1 ) ) */
                    /* ( ( -1 ) ( D22 ) ) */
                    /* = D21 * ( ( D11 ) ( -1 ) ) */
                    /* ( ( -1 ) ( D22 ) ) */
                    i__1 = k - 1 + kw * w_dim1;
                    d21.r = w[i__1].r;
                    d21.i = w[i__1].i; // , expr subst
                    z_div(&z__1, &w[k + kw * w_dim1], &d21);
                    d11.r = z__1.r;
                    d11.i = z__1.i; // , expr subst
                    z_div(&z__1, &w[k - 1 + (kw - 1) * w_dim1], &d21);
                    d22.r = z__1.r;
                    d22.i = z__1.i; // , expr subst
                    z__3.r = d11.r * d22.r - d11.i * d22.i;
                    z__3.i = d11.r * d22.i + d11.i * d22.r; // , expr subst
                    z__2.r = z__3.r - 1.;
                    z__2.i = z__3.i - 0.; // , expr subst
                    z_div(&z__1, &c_b1, &z__2);
                    t.r = z__1.r;
                    t.i = z__1.i; // , expr subst
                    z_div(&z__1, &t, &d21);
                    d21.r = z__1.r;
                    d21.i = z__1.i; // , expr subst
                    /* Update elements in columns A(k-1) and A(k) as */
                    /* dot products of rows of ( W(kw-1) W(kw) ) and columns */
                    /* of D**(-1) */
                    i__1 = k - 2;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = j + (k - 1) * a_dim1;
                        i__3 = j + (kw - 1) * w_dim1;
                        z__3.r = d11.r * w[i__3].r - d11.i * w[i__3].i;
                        z__3.i = d11.r * w[i__3].i + d11.i * w[i__3] .r; // , expr subst
                        i__4 = j + kw * w_dim1;
                        z__2.r = z__3.r - w[i__4].r;
                        z__2.i = z__3.i - w[i__4] .i; // , expr subst
                        z__1.r = d21.r * z__2.r - d21.i * z__2.i;
                        z__1.i = d21.r * z__2.i + d21.i * z__2.r; // , expr subst
                        a[i__2].r = z__1.r;
                        a[i__2].i = z__1.i; // , expr subst
                        i__2 = j + k * a_dim1;
                        i__3 = j + kw * w_dim1;
                        z__3.r = d22.r * w[i__3].r - d22.i * w[i__3].i;
                        z__3.i = d22.r * w[i__3].i + d22.i * w[i__3] .r; // , expr subst
                        i__4 = j + (kw - 1) * w_dim1;
                        z__2.r = z__3.r - w[i__4].r;
                        z__2.i = z__3.i - w[i__4] .i; // , expr subst
                        z__1.r = d21.r * z__2.r - d21.i * z__2.i;
                        z__1.i = d21.r * z__2.i + d21.i * z__2.r; // , expr subst
                        a[i__2].r = z__1.r;
                        a[i__2].i = z__1.i; // , expr subst
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
            ipiv[k] = -kp;
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
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemv_("No transpose", &i__3, &i__4, &z__1, &a[j + (k + 1) * a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b1, &a[j + jj * a_dim1], &c__1);
                /* L40: */
            }
            /* Update the rectangular superdiagonal block */
            i__2 = j - 1;
            i__3 = *n - k;
            z__1.r = -1.;
            z__1.i = -0.; // , expr subst
            zgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &z__1, &a[( k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * w_dim1], ldw, &c_b1, &a[j * a_dim1 + 1], lda);
            /* L50: */
        }
        /* Put U12 in standard form by partially undoing the interchanges */
        /* in columns k+1:n looping backwards from k+1 to n */
        j = k + 1;
L60: /* Undo the interchanges (if any) of rows JJ and JP at each */
        /* step J */
        /* (Here, J is a diagonal index) */
        jj = j;
        jp = ipiv[j];
        if (jp < 0)
        {
            jp = -jp;
            /* (Here, J is a diagonal index) */
            ++j;
        }
        /* (NOTE: Here, J is used to determine row length. Length N-J+1 */
        /* of the rows to swap back doesn't include diagonal element) */
        ++j;
        if (jp != jj && j <= *n)
        {
            i__1 = *n - j + 1;
            zswap_(&i__1, &a[jp + j * a_dim1], lda, &a[jj + j * a_dim1], lda);
        }
        if (j < *n)
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
        /* Copy column K of A to column K of W and update it */
        i__1 = *n - k + 1;
        zcopy_(&i__1, &a[k + k * a_dim1], &c__1, &w[k + k * w_dim1], &c__1);
        i__1 = *n - k + 1;
        i__2 = k - 1;
        z__1.r = -1.;
        z__1.i = -0.; // , expr subst
        zgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &w[k + w_dim1], ldw, &c_b1, &w[k + k * w_dim1], &c__1);
        kstep = 1;
        /* Determine rows and columns to be interchanged and whether */
        /* a 1-by-1 or 2-by-2 pivot block will be used */
        i__1 = k + k * w_dim1;
        absakk = (d__1 = w[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&w[k + k * w_dim1]), f2c_dabs(d__2));
        /* IMAX is the row-index of the largest off-diagonal element in */
        if (k < *n)
        {
            i__1 = *n - k;
            imax = k + izamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
            i__1 = imax + k * w_dim1;
            colmax = (d__1 = w[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&w[imax + k * w_dim1]), f2c_dabs(d__2));
        }
        else
        {
            colmax = 0.;
        }
        if (max(absakk,colmax) == 0.)
        {
            /* Column K is zero or underflow: set INFO and continue */
            if (*info == 0)
            {
                *info = k;
            }
            kp = k;
        }
        else
        {
            if (absakk >= alpha * colmax)
            {
                /* no interchange, use 1-by-1 pivot block */
                kp = k;
            }
            else
            {
                /* Copy column IMAX to column K+1 of W and update it */
                i__1 = imax - k;
                zcopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * w_dim1], &c__1);
                i__1 = *n - imax + 1;
                zcopy_(&i__1, &a[imax + imax * a_dim1], &c__1, &w[imax + (k + 1) * w_dim1], &c__1);
                i__1 = *n - k + 1;
                i__2 = k - 1;
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemv_("No transpose", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &w[imax + w_dim1], ldw, &c_b1, &w[k + (k + 1) * w_dim1], &c__1);
                /* JMAX is the column-index of the largest off-diagonal */
                /* element in row IMAX, and ROWMAX is its absolute value */
                i__1 = imax - k;
                jmax = k - 1 + izamax_(&i__1, &w[k + (k + 1) * w_dim1], &c__1) ;
                i__1 = jmax + (k + 1) * w_dim1;
                rowmax = (d__1 = w[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&w[ jmax + (k + 1) * w_dim1]), f2c_dabs(d__2));
                if (imax < *n)
                {
                    i__1 = *n - imax;
                    jmax = imax + izamax_(&i__1, &w[imax + 1 + (k + 1) * w_dim1], &c__1);
                    /* Computing MAX */
                    i__1 = jmax + (k + 1) * w_dim1;
                    d__3 = rowmax;
                    d__4 = (d__1 = w[i__1].r, f2c_dabs(d__1)) + ( d__2 = d_imag(&w[jmax + (k + 1) * w_dim1]), f2c_dabs( d__2)); // , expr subst
                    rowmax = max(d__3,d__4);
                }
                if (absakk >= alpha * colmax * (colmax / rowmax))
                {
                    /* no interchange, use 1-by-1 pivot block */
                    kp = k;
                }
                else /* if(complicated condition) */
                {
                    i__1 = imax + (k + 1) * w_dim1;
                    if ((d__1 = w[i__1].r, f2c_dabs(d__1)) + (d__2 = d_imag(&w[ imax + (k + 1) * w_dim1]), f2c_dabs(d__2)) >= alpha * rowmax)
                    {
                        /* interchange rows and columns K and IMAX, use 1-by-1 */
                        /* pivot block */
                        kp = imax;
                        /* copy column K+1 of W to column K of W */
                        i__1 = *n - k + 1;
                        zcopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * w_dim1], &c__1);
                    }
                    else
                    {
                        /* interchange rows and columns K+1 and IMAX, use 2-by-2 */
                        /* pivot block */
                        kp = imax;
                        kstep = 2;
                    }
                }
            }
            /* ============================================================ */
            /* KK is the column of A where pivoting step stopped */
            kk = k + kstep - 1;
            /* Interchange rows and columns KP and KK. */
            /* Updated column KP is already stored in column KK of W. */
            if (kp != kk)
            {
                /* Copy non-updated column KK to column KP of submatrix A */
                /* at step K. No need to copy element into column K */
                /* (or K and K+1 for 2-by-2 pivot) of A, since these columns */
                /* will be later overwritten. */
                i__1 = kp + kp * a_dim1;
                i__2 = kk + kk * a_dim1;
                a[i__1].r = a[i__2].r;
                a[i__1].i = a[i__2].i; // , expr subst
                i__1 = kp - kk - 1;
                zcopy_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 1) * a_dim1], lda);
                if (kp < *n)
                {
                    i__1 = *n - kp;
                    zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1], &c__1);
                }
                /* Interchange rows KK and KP in first K-1 columns of A */
                /* (columns K (or K and K+1 for 2-by-2 pivot) of A will be */
                /* later overwritten). Interchange rows KK and KP */
                /* in first KK columns of W. */
                if (k > 1)
                {
                    i__1 = k - 1;
                    zswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
                }
                zswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
            }
            if (kstep == 1)
            {
                /* 1-by-1 pivot block D(k): column k of W now holds */
                /* W(k) = L(k)*D(k), */
                /* where L(k) is the k-th column of L */
                /* Store subdiag. elements of column L(k) */
                /* and 1-by-1 block D(k) in column k of A. */
                /* (NOTE: Diagonal element L(k,k) is a UNIT element */
                /* and not stored) */
                /* A(k,k) := D(k,k) = W(k,k) */
                /* A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k) */
                i__1 = *n - k + 1;
                zcopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], & c__1);
                if (k < *n)
                {
                    z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
                    r1.r = z__1.r;
                    r1.i = z__1.i; // , expr subst
                    i__1 = *n - k;
                    zscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
                }
            }
            else
            {
                /* 2-by-2 pivot block D(k): columns k and k+1 of W now hold */
                /* ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */
                /* where L(k) and L(k+1) are the k-th and (k+1)-th columns */
                /* of L */
                /* Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2 */
                /* block D(k:k+1,k:k+1) in columns k and k+1 of A. */
                /* (NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT */
                /* block and not stored) */
                /* A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1) */
                /* A(k+2:N,k:k+1) := L(k+2:N,k:k+1) = */
                /* = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) ) */
                if (k < *n - 1)
                {
                    /* Compose the columns of the inverse of 2-by-2 pivot */
                    /* block D in the following way to reduce the number */
                    /* of FLOPS when we myltiply panel ( W(k) W(k+1) ) by */
                    /* this inverse */
                    /* D**(-1) = ( d11 d21 )**(-1) = */
                    /* ( d21 d22 ) */
                    /* = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) = */
                    /* ( (-d21 ) ( d11 ) ) */
                    /* = 1/d21 * 1/((d11/d21)*(d22/d21)-1) * */
                    /* * ( ( d22/d21 ) ( -1 ) ) = */
                    /* ( ( -1 ) ( d11/d21 ) ) */
                    /* = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) ( -1 ) ) = */
                    /* ( ( -1 ) ( D22 ) ) */
                    /* = 1/d21 * T * ( ( D11 ) ( -1 ) ) */
                    /* ( ( -1 ) ( D22 ) ) */
                    /* = D21 * ( ( D11 ) ( -1 ) ) */
                    /* ( ( -1 ) ( D22 ) ) */
                    i__1 = k + 1 + k * w_dim1;
                    d21.r = w[i__1].r;
                    d21.i = w[i__1].i; // , expr subst
                    z_div(&z__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
                    d11.r = z__1.r;
                    d11.i = z__1.i; // , expr subst
                    z_div(&z__1, &w[k + k * w_dim1], &d21);
                    d22.r = z__1.r;
                    d22.i = z__1.i; // , expr subst
                    z__3.r = d11.r * d22.r - d11.i * d22.i;
                    z__3.i = d11.r * d22.i + d11.i * d22.r; // , expr subst
                    z__2.r = z__3.r - 1.;
                    z__2.i = z__3.i - 0.; // , expr subst
                    z_div(&z__1, &c_b1, &z__2);
                    t.r = z__1.r;
                    t.i = z__1.i; // , expr subst
                    z_div(&z__1, &t, &d21);
                    d21.r = z__1.r;
                    d21.i = z__1.i; // , expr subst
                    /* Update elements in columns A(k) and A(k+1) as */
                    /* dot products of rows of ( W(k) W(k+1) ) and columns */
                    /* of D**(-1) */
                    i__1 = *n;
                    for (j = k + 2;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = j + k * a_dim1;
                        i__3 = j + k * w_dim1;
                        z__3.r = d11.r * w[i__3].r - d11.i * w[i__3].i;
                        z__3.i = d11.r * w[i__3].i + d11.i * w[i__3] .r; // , expr subst
                        i__4 = j + (k + 1) * w_dim1;
                        z__2.r = z__3.r - w[i__4].r;
                        z__2.i = z__3.i - w[i__4] .i; // , expr subst
                        z__1.r = d21.r * z__2.r - d21.i * z__2.i;
                        z__1.i = d21.r * z__2.i + d21.i * z__2.r; // , expr subst
                        a[i__2].r = z__1.r;
                        a[i__2].i = z__1.i; // , expr subst
                        i__2 = j + (k + 1) * a_dim1;
                        i__3 = j + (k + 1) * w_dim1;
                        z__3.r = d22.r * w[i__3].r - d22.i * w[i__3].i;
                        z__3.i = d22.r * w[i__3].i + d22.i * w[i__3] .r; // , expr subst
                        i__4 = j + k * w_dim1;
                        z__2.r = z__3.r - w[i__4].r;
                        z__2.i = z__3.i - w[i__4] .i; // , expr subst
                        z__1.r = d21.r * z__2.r - d21.i * z__2.i;
                        z__1.i = d21.r * z__2.i + d21.i * z__2.r; // , expr subst
                        a[i__2].r = z__1.r;
                        a[i__2].i = z__1.i; // , expr subst
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
            ipiv[k] = -kp;
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
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemv_("No transpose", &i__4, &i__5, &z__1, &a[jj + a_dim1], lda, &w[jj + w_dim1], ldw, &c_b1, &a[jj + jj * a_dim1], &c__1);
                /* L100: */
            }
            /* Update the rectangular subdiagonal block */
            if (j + jb <= *n)
            {
                i__3 = *n - j - jb + 1;
                i__4 = k - 1;
                z__1.r = -1.;
                z__1.i = -0.; // , expr subst
                zgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &z__1, &a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b1, &a[j + jb + j * a_dim1], lda);
            }
            /* L110: */
        }
        /* Put L21 in standard form by partially undoing the interchanges */
        /* of rows in columns 1:k-1 looping backwards from k-1 to 1 */
        j = k - 1;
L120: /* Undo the interchanges (if any) of rows JJ and JP at each */
        /* step J */
        /* (Here, J is a diagonal index) */
        jj = j;
        jp = ipiv[j];
        if (jp < 0)
        {
            jp = -jp;
            /* (Here, J is a diagonal index) */
            --j;
        }
        /* (NOTE: Here, J is used to determine row length. Length J */
        /* of the rows to swap back doesn't include diagonal element) */
        --j;
        if (jp != jj && j >= 1)
        {
            zswap_(&j, &a[jp + a_dim1], lda, &a[jj + a_dim1], lda);
        }
        if (j > 1)
        {
            goto L120;
        }
        /* Set KB to the number of columns factorized */
        *kb = k - 1;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of ZLASYF */
}
/* zlasyf_ */
