/* ../netlib/zhetf2.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZHETF2 computes the factorization of a complex Hermitian matrix, using the diagonal pivoting me thod (unblocked algorithm, calling Level 2 BLAS). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZHETF2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetf2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetf2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetf2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZHETF2( UPLO, N, A, LDA, IPIV, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETF2 computes the factorization of a complex Hermitian matrix A */
/* > using the Bunch-Kaufman diagonal pivoting method: */
/* > */
/* > A = U*D*U**H or A = L*D*L**H */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, U**H is the conjugate transpose of U, and D is */
/* > Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > This is the unblocked version of the algorithm, calling Level 2 BLAS. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the upper or lower triangular part of the */
/* > Hermitian matrix A is stored: */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the Hermitian matrix A. If UPLO = 'U', the leading */
/* > n-by-n upper triangular part of A contains the upper */
/* > triangular part of the matrix A, and the strictly lower */
/* > triangular part of A is not referenced. If UPLO = 'L', the */
/* > leading n-by-n lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, the block diagonal matrix D and the multipliers used */
/* > to obtain the factor U or L (see below for further details). */
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
/* > If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* > interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* > If IPIV(k) = IPIV(k-1) < 0, then rows and columns */
/* > k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
/* > is a 2-by-2 diagonal block. */
/* > */
/* > If UPLO = 'L': */
/* > If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* > interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > */
/* > If IPIV(k) = IPIV(k+1) < 0, then rows and columns */
/* > k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1) */
/* > is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -k, the k-th argument had an illegal value */
/* > > 0: if INFO = k, D(k,k) is exactly zero. The factorization */
/* > has been completed, but the block diagonal matrix D is */
/* > exactly singular, and division by zero will occur if it */
/* > is used to solve a system of equations. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2013 */
/* > \ingroup complex16HEcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > If UPLO = 'U', then A = U*D*U**H, where */
/* > U = P(n)*U(n)* ... *P(k)U(k)* ..., */
/* > i.e., U is a product of terms P(k)*U(k), where k decreases from n to */
/* > 1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1 */
/* > and 2-by-2 diagonal blocks D(k). P(k) is a permutation matrix as */
/* > defined by IPIV(k), and U(k) is a unit upper triangular matrix, such */
/* > that if the diagonal block D(k) is of order s (s = 1 or 2), then */
/* > */
/* > ( I v 0 ) k-s */
/* > U(k) = ( 0 I 0 ) s */
/* > ( 0 0 I ) n-k */
/* > k-s s n-k */
/* > */
/* > If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k). */
/* > If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k), */
/* > and A(k,k), and v overwrites A(1:k-2,k-1:k). */
/* > */
/* > If UPLO = 'L', then A = L*D*L**H, where */
/* > L = P(1)*L(1)* ... *P(k)*L(k)* ..., */
/* > i.e., L is a product of terms P(k)*L(k), where k increases from 1 to */
/* > n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1 */
/* > and 2-by-2 diagonal blocks D(k). P(k) is a permutation matrix as */
/* > defined by IPIV(k), and L(k) is a unit lower triangular matrix, such */
/* > that if the diagonal block D(k) is of order s (s = 1 or 2), then */
/* > */
/* > ( I 0 0 ) k-1 */
/* > L(k) = ( 0 I 0 ) s */
/* > ( 0 v I ) n-k-s+1 */
/* > k-1 s n-k-s+1 */
/* > */
/* > If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k). */
/* > If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k), */
/* > and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1). */
/* > \endverbatim */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > 09-29-06 - patch from */
/* > Bobby Cheng, MathWorks */
/* > */
/* > Replace l.210 and l.393 */
/* > IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN */
/* > by */
/* > IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN */
/* > */
/* > 01-01-96 - Based on modifications by */
/* > J. Lewis, Boeing Computer Services Company */
/* > A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
int zhetf2_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;
    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    doublereal d__;
    integer i__, j, k;
    doublecomplex t;
    doublereal r1, d11;
    doublecomplex d12;
    doublereal d22;
    doublecomplex d21;
    integer kk, kp;
    doublecomplex wk;
    doublereal tt;
    doublecomplex wkm1, wkp1;
    integer imax, jmax;
    extern /* Subroutine */
    int zher_(char *, integer *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *);
    doublereal alpha;
    extern logical lsame_(char *, char *);
    integer kstep;
    logical upper;
    extern /* Subroutine */
    int zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    doublereal absakk;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */
    int xerbla_(char *, integer *), zdscal_( integer *, doublereal *, doublecomplex *, integer *);
    doublereal colmax;
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
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
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
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZHETF2", &i__1);
        return 0;
    }
    /* Initialize ALPHA for use in choosing pivot block size. */
    alpha = (sqrt(17.) + 1.) / 8.;
    if (upper)
    {
        /* Factorize A as U*D*U**H using the upper triangle of A */
        /* K is the main loop index, decreasing from N to 1 in steps of */
        /* 1 or 2 */
        k = *n;
L10: /* If K < 1, exit from loop */
        if (k < 1)
        {
            goto L90;
        }
        kstep = 1;
        /* Determine rows and columns to be interchanged and whether */
        /* a 1-by-1 or 2-by-2 pivot block will be used */
        i__1 = k + k * a_dim1;
        absakk = (d__1 = a[i__1].r, f2c_abs(d__1));
        /* IMAX is the row-index of the largest off-diagonal element in */
        /* column K, and COLMAX is its absolute value. */
        /* Determine both COLMAX and IMAX. */
        if (k > 1)
        {
            i__1 = k - 1;
            imax = izamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
            i__1 = imax + k * a_dim1;
            colmax = (d__1 = a[i__1].r, f2c_abs(d__1)) + (d__2 = d_imag(&a[imax + k * a_dim1]), f2c_abs(d__2));
        }
        else
        {
            colmax = 0.;
        }
        if (max(absakk,colmax) == 0. || disnan_(&absakk))
        {
            /* Column K is zero or underflow, or contains a NaN: */
            /* set INFO and continue */
            if (*info == 0)
            {
                *info = k;
            }
            kp = k;
            i__1 = k + k * a_dim1;
            i__2 = k + k * a_dim1;
            d__1 = a[i__2].r;
            a[i__1].r = d__1;
            a[i__1].i = 0.; // , expr subst
        }
        else
        {
            /* ============================================================ */
            /* Test for interchange */
            if (absakk >= alpha * colmax)
            {
                /* no interchange, use 1-by-1 pivot block */
                kp = k;
            }
            else
            {
                /* JMAX is the column-index of the largest off-diagonal */
                /* element in row IMAX, and ROWMAX is its absolute value. */
                /* Determine only ROWMAX. */
                i__1 = k - imax;
                jmax = imax + izamax_(&i__1, &a[imax + (imax + 1) * a_dim1], lda);
                i__1 = imax + jmax * a_dim1;
                rowmax = (d__1 = a[i__1].r, f2c_abs(d__1)) + (d__2 = d_imag(&a[ imax + jmax * a_dim1]), f2c_abs(d__2));
                if (imax > 1)
                {
                    i__1 = imax - 1;
                    jmax = izamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
                    /* Computing MAX */
                    i__1 = jmax + imax * a_dim1;
                    d__3 = rowmax;
                    d__4 = (d__1 = a[i__1].r, f2c_abs(d__1)) + ( d__2 = d_imag(&a[jmax + imax * a_dim1]), f2c_abs(d__2) ); // , expr subst
                    rowmax = max(d__3,d__4);
                }
                if (absakk >= alpha * colmax * (colmax / rowmax))
                {
                    /* no interchange, use 1-by-1 pivot block */
                    kp = k;
                }
                else /* if(complicated condition) */
                {
                    i__1 = imax + imax * a_dim1;
                    if ((d__1 = a[i__1].r, f2c_abs(d__1)) >= alpha * rowmax)
                    {
                        /* interchange rows and columns K and IMAX, use 1-by-1 */
                        /* pivot block */
                        kp = imax;
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
            kk = k - kstep + 1;
            if (kp != kk)
            {
                /* Interchange rows and columns KK and KP in the leading */
                /* submatrix A(1:k,1:k) */
                i__1 = kp - 1;
                zswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
                i__1 = kk - 1;
                for (j = kp + 1;
                        j <= i__1;
                        ++j)
                {
                    d_cnjg(&z__1, &a[j + kk * a_dim1]);
                    t.r = z__1.r;
                    t.i = z__1.i; // , expr subst
                    i__2 = j + kk * a_dim1;
                    d_cnjg(&z__1, &a[kp + j * a_dim1]);
                    a[i__2].r = z__1.r;
                    a[i__2].i = z__1.i; // , expr subst
                    i__2 = kp + j * a_dim1;
                    a[i__2].r = t.r;
                    a[i__2].i = t.i; // , expr subst
                    /* L20: */
                }
                i__1 = kp + kk * a_dim1;
                d_cnjg(&z__1, &a[kp + kk * a_dim1]);
                a[i__1].r = z__1.r;
                a[i__1].i = z__1.i; // , expr subst
                i__1 = kk + kk * a_dim1;
                r1 = a[i__1].r;
                i__1 = kk + kk * a_dim1;
                i__2 = kp + kp * a_dim1;
                d__1 = a[i__2].r;
                a[i__1].r = d__1;
                a[i__1].i = 0.; // , expr subst
                i__1 = kp + kp * a_dim1;
                a[i__1].r = r1;
                a[i__1].i = 0.; // , expr subst
                if (kstep == 2)
                {
                    i__1 = k + k * a_dim1;
                    i__2 = k + k * a_dim1;
                    d__1 = a[i__2].r;
                    a[i__1].r = d__1;
                    a[i__1].i = 0.; // , expr subst
                    i__1 = k - 1 + k * a_dim1;
                    t.r = a[i__1].r;
                    t.i = a[i__1].i; // , expr subst
                    i__1 = k - 1 + k * a_dim1;
                    i__2 = kp + k * a_dim1;
                    a[i__1].r = a[i__2].r;
                    a[i__1].i = a[i__2].i; // , expr subst
                    i__1 = kp + k * a_dim1;
                    a[i__1].r = t.r;
                    a[i__1].i = t.i; // , expr subst
                }
            }
            else
            {
                i__1 = k + k * a_dim1;
                i__2 = k + k * a_dim1;
                d__1 = a[i__2].r;
                a[i__1].r = d__1;
                a[i__1].i = 0.; // , expr subst
                if (kstep == 2)
                {
                    i__1 = k - 1 + (k - 1) * a_dim1;
                    i__2 = k - 1 + (k - 1) * a_dim1;
                    d__1 = a[i__2].r;
                    a[i__1].r = d__1;
                    a[i__1].i = 0.; // , expr subst
                }
            }
            /* Update the leading submatrix */
            if (kstep == 1)
            {
                /* 1-by-1 pivot block D(k): column k now holds */
                /* W(k) = U(k)*D(k) */
                /* where U(k) is the k-th column of U */
                /* Perform a rank-1 update of A(1:k-1,1:k-1) as */
                /* A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H */
                i__1 = k + k * a_dim1;
                r1 = 1. / a[i__1].r;
                i__1 = k - 1;
                d__1 = -r1;
                zher_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &a[ a_offset], lda);
                /* Store U(k) in column k */
                i__1 = k - 1;
                zdscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
            }
            else
            {
                /* 2-by-2 pivot block D(k): columns k and k-1 now hold */
                /* ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */
                /* where U(k) and U(k-1) are the k-th and (k-1)-th columns */
                /* of U */
                /* Perform a rank-2 update of A(1:k-2,1:k-2) as */
                /* A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H */
                /* = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H */
                if (k > 2)
                {
                    i__1 = k - 1 + k * a_dim1;
                    d__1 = a[i__1].r;
                    d__2 = d_imag(&a[k - 1 + k * a_dim1]);
                    d__ = dlapy2_(&d__1, &d__2);
                    i__1 = k - 1 + (k - 1) * a_dim1;
                    d22 = a[i__1].r / d__;
                    i__1 = k + k * a_dim1;
                    d11 = a[i__1].r / d__;
                    tt = 1. / (d11 * d22 - 1.);
                    i__1 = k - 1 + k * a_dim1;
                    z__1.r = a[i__1].r / d__;
                    z__1.i = a[i__1].i / d__; // , expr subst
                    d12.r = z__1.r;
                    d12.i = z__1.i; // , expr subst
                    d__ = tt / d__;
                    for (j = k - 2;
                            j >= 1;
                            --j)
                    {
                        i__1 = j + (k - 1) * a_dim1;
                        z__3.r = d11 * a[i__1].r;
                        z__3.i = d11 * a[i__1].i; // , expr subst
                        d_cnjg(&z__5, &d12);
                        i__2 = j + k * a_dim1;
                        z__4.r = z__5.r * a[i__2].r - z__5.i * a[i__2].i;
                        z__4.i = z__5.r * a[i__2].i + z__5.i * a[i__2] .r; // , expr subst
                        z__2.r = z__3.r - z__4.r;
                        z__2.i = z__3.i - z__4.i; // , expr subst
                        z__1.r = d__ * z__2.r;
                        z__1.i = d__ * z__2.i; // , expr subst
                        wkm1.r = z__1.r;
                        wkm1.i = z__1.i; // , expr subst
                        i__1 = j + k * a_dim1;
                        z__3.r = d22 * a[i__1].r;
                        z__3.i = d22 * a[i__1].i; // , expr subst
                        i__2 = j + (k - 1) * a_dim1;
                        z__4.r = d12.r * a[i__2].r - d12.i * a[i__2].i;
                        z__4.i = d12.r * a[i__2].i + d12.i * a[i__2] .r; // , expr subst
                        z__2.r = z__3.r - z__4.r;
                        z__2.i = z__3.i - z__4.i; // , expr subst
                        z__1.r = d__ * z__2.r;
                        z__1.i = d__ * z__2.i; // , expr subst
                        wk.r = z__1.r;
                        wk.i = z__1.i; // , expr subst
                        for (i__ = j;
                                i__ >= 1;
                                --i__)
                        {
                            i__1 = i__ + j * a_dim1;
                            i__2 = i__ + j * a_dim1;
                            i__3 = i__ + k * a_dim1;
                            d_cnjg(&z__4, &wk);
                            z__3.r = a[i__3].r * z__4.r - a[i__3].i * z__4.i;
                            z__3.i = a[i__3].r * z__4.i + a[i__3].i * z__4.r; // , expr subst
                            z__2.r = a[i__2].r - z__3.r;
                            z__2.i = a[i__2].i - z__3.i; // , expr subst
                            i__4 = i__ + (k - 1) * a_dim1;
                            d_cnjg(&z__6, &wkm1);
                            z__5.r = a[i__4].r * z__6.r - a[i__4].i * z__6.i;
                            z__5.i = a[i__4].r * z__6.i + a[i__4].i * z__6.r; // , expr subst
                            z__1.r = z__2.r - z__5.r;
                            z__1.i = z__2.i - z__5.i; // , expr subst
                            a[i__1].r = z__1.r;
                            a[i__1].i = z__1.i; // , expr subst
                            /* L30: */
                        }
                        i__1 = j + k * a_dim1;
                        a[i__1].r = wk.r;
                        a[i__1].i = wk.i; // , expr subst
                        i__1 = j + (k - 1) * a_dim1;
                        a[i__1].r = wkm1.r;
                        a[i__1].i = wkm1.i; // , expr subst
                        i__1 = j + j * a_dim1;
                        i__2 = j + j * a_dim1;
                        d__1 = a[i__2].r;
                        z__1.r = d__1;
                        z__1.i = 0.; // , expr subst
                        a[i__1].r = z__1.r;
                        a[i__1].i = z__1.i; // , expr subst
                        /* L40: */
                    }
                }
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
    }
    else
    {
        /* Factorize A as L*D*L**H using the lower triangle of A */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2 */
        k = 1;
L50: /* If K > N, exit from loop */
        if (k > *n)
        {
            goto L90;
        }
        kstep = 1;
        /* Determine rows and columns to be interchanged and whether */
        /* a 1-by-1 or 2-by-2 pivot block will be used */
        i__1 = k + k * a_dim1;
        absakk = (d__1 = a[i__1].r, f2c_abs(d__1));
        /* IMAX is the row-index of the largest off-diagonal element in */
        /* column K, and COLMAX is its absolute value. */
        /* Determine both COLMAX and IMAX. */
        if (k < *n)
        {
            i__1 = *n - k;
            imax = k + izamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
            i__1 = imax + k * a_dim1;
            colmax = (d__1 = a[i__1].r, f2c_abs(d__1)) + (d__2 = d_imag(&a[imax + k * a_dim1]), f2c_abs(d__2));
        }
        else
        {
            colmax = 0.;
        }
        if (max(absakk,colmax) == 0. || disnan_(&absakk))
        {
            /* Column K is zero or underflow, or contains a NaN: */
            /* set INFO and continue */
            if (*info == 0)
            {
                *info = k;
            }
            kp = k;
            i__1 = k + k * a_dim1;
            i__2 = k + k * a_dim1;
            d__1 = a[i__2].r;
            a[i__1].r = d__1;
            a[i__1].i = 0.; // , expr subst
        }
        else
        {
            /* ============================================================ */
            /* Test for interchange */
            if (absakk >= alpha * colmax)
            {
                /* no interchange, use 1-by-1 pivot block */
                kp = k;
            }
            else
            {
                /* JMAX is the column-index of the largest off-diagonal */
                /* element in row IMAX, and ROWMAX is its absolute value. */
                /* Determine only ROWMAX. */
                i__1 = imax - k;
                jmax = k - 1 + izamax_(&i__1, &a[imax + k * a_dim1], lda);
                i__1 = imax + jmax * a_dim1;
                rowmax = (d__1 = a[i__1].r, f2c_abs(d__1)) + (d__2 = d_imag(&a[ imax + jmax * a_dim1]), f2c_abs(d__2));
                if (imax < *n)
                {
                    i__1 = *n - imax;
                    jmax = imax + izamax_(&i__1, &a[imax + 1 + imax * a_dim1], &c__1);
                    /* Computing MAX */
                    i__1 = jmax + imax * a_dim1;
                    d__3 = rowmax;
                    d__4 = (d__1 = a[i__1].r, f2c_abs(d__1)) + ( d__2 = d_imag(&a[jmax + imax * a_dim1]), f2c_abs(d__2) ); // , expr subst
                    rowmax = max(d__3,d__4);
                }
                if (absakk >= alpha * colmax * (colmax / rowmax))
                {
                    /* no interchange, use 1-by-1 pivot block */
                    kp = k;
                }
                else /* if(complicated condition) */
                {
                    i__1 = imax + imax * a_dim1;
                    if ((d__1 = a[i__1].r, f2c_abs(d__1)) >= alpha * rowmax)
                    {
                        /* interchange rows and columns K and IMAX, use 1-by-1 */
                        /* pivot block */
                        kp = imax;
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
            kk = k + kstep - 1;
            if (kp != kk)
            {
                /* Interchange rows and columns KK and KP in the trailing */
                /* submatrix A(k:n,k:n) */
                if (kp < *n)
                {
                    i__1 = *n - kp;
                    zswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1], &c__1);
                }
                i__1 = kp - 1;
                for (j = kk + 1;
                        j <= i__1;
                        ++j)
                {
                    d_cnjg(&z__1, &a[j + kk * a_dim1]);
                    t.r = z__1.r;
                    t.i = z__1.i; // , expr subst
                    i__2 = j + kk * a_dim1;
                    d_cnjg(&z__1, &a[kp + j * a_dim1]);
                    a[i__2].r = z__1.r;
                    a[i__2].i = z__1.i; // , expr subst
                    i__2 = kp + j * a_dim1;
                    a[i__2].r = t.r;
                    a[i__2].i = t.i; // , expr subst
                    /* L60: */
                }
                i__1 = kp + kk * a_dim1;
                d_cnjg(&z__1, &a[kp + kk * a_dim1]);
                a[i__1].r = z__1.r;
                a[i__1].i = z__1.i; // , expr subst
                i__1 = kk + kk * a_dim1;
                r1 = a[i__1].r;
                i__1 = kk + kk * a_dim1;
                i__2 = kp + kp * a_dim1;
                d__1 = a[i__2].r;
                a[i__1].r = d__1;
                a[i__1].i = 0.; // , expr subst
                i__1 = kp + kp * a_dim1;
                a[i__1].r = r1;
                a[i__1].i = 0.; // , expr subst
                if (kstep == 2)
                {
                    i__1 = k + k * a_dim1;
                    i__2 = k + k * a_dim1;
                    d__1 = a[i__2].r;
                    a[i__1].r = d__1;
                    a[i__1].i = 0.; // , expr subst
                    i__1 = k + 1 + k * a_dim1;
                    t.r = a[i__1].r;
                    t.i = a[i__1].i; // , expr subst
                    i__1 = k + 1 + k * a_dim1;
                    i__2 = kp + k * a_dim1;
                    a[i__1].r = a[i__2].r;
                    a[i__1].i = a[i__2].i; // , expr subst
                    i__1 = kp + k * a_dim1;
                    a[i__1].r = t.r;
                    a[i__1].i = t.i; // , expr subst
                }
            }
            else
            {
                i__1 = k + k * a_dim1;
                i__2 = k + k * a_dim1;
                d__1 = a[i__2].r;
                a[i__1].r = d__1;
                a[i__1].i = 0.; // , expr subst
                if (kstep == 2)
                {
                    i__1 = k + 1 + (k + 1) * a_dim1;
                    i__2 = k + 1 + (k + 1) * a_dim1;
                    d__1 = a[i__2].r;
                    a[i__1].r = d__1;
                    a[i__1].i = 0.; // , expr subst
                }
            }
            /* Update the trailing submatrix */
            if (kstep == 1)
            {
                /* 1-by-1 pivot block D(k): column k now holds */
                /* W(k) = L(k)*D(k) */
                /* where L(k) is the k-th column of L */
                if (k < *n)
                {
                    /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
                    /* A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H */
                    i__1 = k + k * a_dim1;
                    r1 = 1. / a[i__1].r;
                    i__1 = *n - k;
                    d__1 = -r1;
                    zher_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &c__1, & a[k + 1 + (k + 1) * a_dim1], lda);
                    /* Store L(k) in column K */
                    i__1 = *n - k;
                    zdscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
                }
            }
            else
            {
                /* 2-by-2 pivot block D(k) */
                if (k < *n - 1)
                {
                    /* Perform a rank-2 update of A(k+2:n,k+2:n) as */
                    /* A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H */
                    /* = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H */
                    /* where L(k) and L(k+1) are the k-th and (k+1)-th */
                    /* columns of L */
                    i__1 = k + 1 + k * a_dim1;
                    d__1 = a[i__1].r;
                    d__2 = d_imag(&a[k + 1 + k * a_dim1]);
                    d__ = dlapy2_(&d__1, &d__2);
                    i__1 = k + 1 + (k + 1) * a_dim1;
                    d11 = a[i__1].r / d__;
                    i__1 = k + k * a_dim1;
                    d22 = a[i__1].r / d__;
                    tt = 1. / (d11 * d22 - 1.);
                    i__1 = k + 1 + k * a_dim1;
                    z__1.r = a[i__1].r / d__;
                    z__1.i = a[i__1].i / d__; // , expr subst
                    d21.r = z__1.r;
                    d21.i = z__1.i; // , expr subst
                    d__ = tt / d__;
                    i__1 = *n;
                    for (j = k + 2;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = j + k * a_dim1;
                        z__3.r = d11 * a[i__2].r;
                        z__3.i = d11 * a[i__2].i; // , expr subst
                        i__3 = j + (k + 1) * a_dim1;
                        z__4.r = d21.r * a[i__3].r - d21.i * a[i__3].i;
                        z__4.i = d21.r * a[i__3].i + d21.i * a[i__3] .r; // , expr subst
                        z__2.r = z__3.r - z__4.r;
                        z__2.i = z__3.i - z__4.i; // , expr subst
                        z__1.r = d__ * z__2.r;
                        z__1.i = d__ * z__2.i; // , expr subst
                        wk.r = z__1.r;
                        wk.i = z__1.i; // , expr subst
                        i__2 = j + (k + 1) * a_dim1;
                        z__3.r = d22 * a[i__2].r;
                        z__3.i = d22 * a[i__2].i; // , expr subst
                        d_cnjg(&z__5, &d21);
                        i__3 = j + k * a_dim1;
                        z__4.r = z__5.r * a[i__3].r - z__5.i * a[i__3].i;
                        z__4.i = z__5.r * a[i__3].i + z__5.i * a[i__3] .r; // , expr subst
                        z__2.r = z__3.r - z__4.r;
                        z__2.i = z__3.i - z__4.i; // , expr subst
                        z__1.r = d__ * z__2.r;
                        z__1.i = d__ * z__2.i; // , expr subst
                        wkp1.r = z__1.r;
                        wkp1.i = z__1.i; // , expr subst
                        i__2 = *n;
                        for (i__ = j;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + j * a_dim1;
                            i__4 = i__ + j * a_dim1;
                            i__5 = i__ + k * a_dim1;
                            d_cnjg(&z__4, &wk);
                            z__3.r = a[i__5].r * z__4.r - a[i__5].i * z__4.i;
                            z__3.i = a[i__5].r * z__4.i + a[i__5].i * z__4.r; // , expr subst
                            z__2.r = a[i__4].r - z__3.r;
                            z__2.i = a[i__4].i - z__3.i; // , expr subst
                            i__6 = i__ + (k + 1) * a_dim1;
                            d_cnjg(&z__6, &wkp1);
                            z__5.r = a[i__6].r * z__6.r - a[i__6].i * z__6.i;
                            z__5.i = a[i__6].r * z__6.i + a[i__6].i * z__6.r; // , expr subst
                            z__1.r = z__2.r - z__5.r;
                            z__1.i = z__2.i - z__5.i; // , expr subst
                            a[i__3].r = z__1.r;
                            a[i__3].i = z__1.i; // , expr subst
                            /* L70: */
                        }
                        i__2 = j + k * a_dim1;
                        a[i__2].r = wk.r;
                        a[i__2].i = wk.i; // , expr subst
                        i__2 = j + (k + 1) * a_dim1;
                        a[i__2].r = wkp1.r;
                        a[i__2].i = wkp1.i; // , expr subst
                        i__2 = j + j * a_dim1;
                        i__3 = j + j * a_dim1;
                        d__1 = a[i__3].r;
                        z__1.r = d__1;
                        z__1.i = 0.; // , expr subst
                        a[i__2].r = z__1.r;
                        a[i__2].i = z__1.i; // , expr subst
                        /* L80: */
                    }
                }
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
        goto L50;
    }
L90:
    return 0;
    /* End of ZHETF2 */
}
/* zhetf2_ */
