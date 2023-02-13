/* ../netlib/ssptrf.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SSPTRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSPTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssptrf. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssptrf. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssptrf. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSPTRF( UPLO, N, AP, IPIV, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL AP( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPTRF computes the factorization of a real symmetric matrix A stored */
/* > in packed format using the Bunch-Kaufman diagonal pivoting method: */
/* > */
/* > A = U*D*U**T or A = L*D*L**T */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is symmetric and block diagonal with */
/* > 1-by-1 and 2-by-2 diagonal blocks. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
*/
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is REAL array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the symmetric matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*/
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* > On exit, the block diagonal matrix D and the multipliers used */
/* > to obtain the factor U or L, stored as a packed triangular */
/* > matrix overwriting A (see below for further details). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D. */
/* > If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* > interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and */
/* > columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
/* > is a 2-by-2 diagonal block. If UPLO = 'L' and IPIV(k) = */
/* > IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were */
/* > interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, D(i,i) is exactly zero. The factorization */
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
/* > \date November 2011 */
/* > \ingroup realOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > 5-96 - Based on modifications by J. Lewis, Boeing Computer Services */
/* > Company */
/* > */
/* > If UPLO = 'U', then A = U*D*U**T, where */
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
/* > If UPLO = 'L', then A = L*D*L**T, where */
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
/* > */
/* ===================================================================== */
/* Subroutine */
int ssptrf_(char *uplo, integer *n, real *ap, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j, k;
    real t, r1, d11, d12, d21, d22;
    integer kc, kk, kp;
    real wk;
    integer kx, knc, kpc, npp;
    real wkm1, wkp1;
    integer imax, jmax;
    extern /* Subroutine */
    int sspr_(char *, integer *, real *, real *, integer *, real *);
    real alpha;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */
    int sscal_(integer *, real *, real *, integer *);
    integer kstep;
    logical upper;
    extern /* Subroutine */
    int sswap_(integer *, real *, integer *, real *, integer *);
    real absakk;
    extern /* Subroutine */
    int xerbla_(char *, integer *);
    extern integer isamax_(integer *, real *, integer *);
    real colmax, rowmax;
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
    /* Test the input parameters. */
    /* Parameter adjustments */
    --ipiv;
    --ap;
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
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSPTRF", &i__1);
        return 0;
    }
    /* Initialize ALPHA for use in choosing pivot block size. */
    alpha = (sqrt(17.f) + 1.f) / 8.f;
    if (upper)
    {
        /* Factorize A as U*D*U**T using the upper triangle of A */
        /* K is the main loop index, decreasing from N to 1 in steps of */
        /* 1 or 2 */
        k = *n;
        kc = (*n - 1) * *n / 2 + 1;
L10:
        knc = kc;
        /* If K < 1, exit from loop */
        if (k < 1)
        {
            goto L110;
        }
        kstep = 1;
        /* Determine rows and columns to be interchanged and whether */
        /* a 1-by-1 or 2-by-2 pivot block will be used */
        absakk = (r__1 = ap[kc + k - 1], f2c_abs(r__1));
        /* IMAX is the row-index of the largest off-diagonal element in */
        /* column K, and COLMAX is its absolute value */
        if (k > 1)
        {
            i__1 = k - 1;
            imax = isamax_(&i__1, &ap[kc], &c__1);
            colmax = (r__1 = ap[kc + imax - 1], f2c_abs(r__1));
        }
        else
        {
            colmax = 0.f;
        }
        if (max(absakk,colmax) == 0.f)
        {
            /* Column K is zero: set INFO and continue */
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
                rowmax = 0.f;
                jmax = imax;
                kx = imax * (imax + 1) / 2 + imax;
                i__1 = k;
                for (j = imax + 1;
                        j <= i__1;
                        ++j)
                {
                    if ((r__1 = ap[kx], f2c_abs(r__1)) > rowmax)
                    {
                        rowmax = (r__1 = ap[kx], f2c_abs(r__1));
                        jmax = j;
                    }
                    kx += j;
                    /* L20: */
                }
                kpc = (imax - 1) * imax / 2 + 1;
                if (imax > 1)
                {
                    i__1 = imax - 1;
                    jmax = isamax_(&i__1, &ap[kpc], &c__1);
                    /* Computing MAX */
                    r__2 = rowmax;
                    r__3 = (r__1 = ap[kpc + jmax - 1], f2c_abs( r__1)); // , expr subst
                    rowmax = max(r__2,r__3);
                }
                if (absakk >= alpha * colmax * (colmax / rowmax))
                {
                    /* no interchange, use 1-by-1 pivot block */
                    kp = k;
                }
                else if ((r__1 = ap[kpc + imax - 1], f2c_abs(r__1)) >= alpha * rowmax)
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
            kk = k - kstep + 1;
            if (kstep == 2)
            {
                knc = knc - k + 1;
            }
            if (kp != kk)
            {
                /* Interchange rows and columns KK and KP in the leading */
                /* submatrix A(1:k,1:k) */
                i__1 = kp - 1;
                sswap_(&i__1, &ap[knc], &c__1, &ap[kpc], &c__1);
                kx = kpc + kp - 1;
                i__1 = kk - 1;
                for (j = kp + 1;
                        j <= i__1;
                        ++j)
                {
                    kx = kx + j - 1;
                    t = ap[knc + j - 1];
                    ap[knc + j - 1] = ap[kx];
                    ap[kx] = t;
                    /* L30: */
                }
                t = ap[knc + kk - 1];
                ap[knc + kk - 1] = ap[kpc + kp - 1];
                ap[kpc + kp - 1] = t;
                if (kstep == 2)
                {
                    t = ap[kc + k - 2];
                    ap[kc + k - 2] = ap[kc + kp - 1];
                    ap[kc + kp - 1] = t;
                }
            }
            /* Update the leading submatrix */
            if (kstep == 1)
            {
                /* 1-by-1 pivot block D(k): column k now holds */
                /* W(k) = U(k)*D(k) */
                /* where U(k) is the k-th column of U */
                /* Perform a rank-1 update of A(1:k-1,1:k-1) as */
                /* A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T */
                r1 = 1.f / ap[kc + k - 1];
                i__1 = k - 1;
                r__1 = -r1;
                sspr_(uplo, &i__1, &r__1, &ap[kc], &c__1, &ap[1]);
                /* Store U(k) in column k */
                i__1 = k - 1;
                sscal_(&i__1, &r1, &ap[kc], &c__1);
            }
            else
            {
                /* 2-by-2 pivot block D(k): columns k and k-1 now hold */
                /* ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */
                /* where U(k) and U(k-1) are the k-th and (k-1)-th columns */
                /* of U */
                /* Perform a rank-2 update of A(1:k-2,1:k-2) as */
                /* A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
                /* = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T */
                if (k > 2)
                {
                    d12 = ap[k - 1 + (k - 1) * k / 2];
                    d22 = ap[k - 1 + (k - 2) * (k - 1) / 2] / d12;
                    d11 = ap[k + (k - 1) * k / 2] / d12;
                    t = 1.f / (d11 * d22 - 1.f);
                    d12 = t / d12;
                    for (j = k - 2;
                            j >= 1;
                            --j)
                    {
                        wkm1 = d12 * (d11 * ap[j + (k - 2) * (k - 1) / 2] - ap[j + (k - 1) * k / 2]);
                        wk = d12 * (d22 * ap[j + (k - 1) * k / 2] - ap[j + (k - 2) * (k - 1) / 2]);
                        for (i__ = j;
                                i__ >= 1;
                                --i__)
                        {
                            ap[i__ + (j - 1) * j / 2] = ap[i__ + (j - 1) * j / 2] - ap[i__ + (k - 1) * k / 2] * wk - ap[ i__ + (k - 2) * (k - 1) / 2] * wkm1;
                            /* L40: */
                        }
                        ap[j + (k - 1) * k / 2] = wk;
                        ap[j + (k - 2) * (k - 1) / 2] = wkm1;
                        /* L50: */
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
        kc = knc - k;
        goto L10;
    }
    else
    {
        /* Factorize A as L*D*L**T using the lower triangle of A */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2 */
        k = 1;
        kc = 1;
        npp = *n * (*n + 1) / 2;
L60:
        knc = kc;
        /* If K > N, exit from loop */
        if (k > *n)
        {
            goto L110;
        }
        kstep = 1;
        /* Determine rows and columns to be interchanged and whether */
        /* a 1-by-1 or 2-by-2 pivot block will be used */
        absakk = (r__1 = ap[kc], f2c_abs(r__1));
        /* IMAX is the row-index of the largest off-diagonal element in */
        /* column K, and COLMAX is its absolute value */
        if (k < *n)
        {
            i__1 = *n - k;
            imax = k + isamax_(&i__1, &ap[kc + 1], &c__1);
            colmax = (r__1 = ap[kc + imax - k], f2c_abs(r__1));
        }
        else
        {
            colmax = 0.f;
        }
        if (max(absakk,colmax) == 0.f)
        {
            /* Column K is zero: set INFO and continue */
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
                /* JMAX is the column-index of the largest off-diagonal */
                /* element in row IMAX, and ROWMAX is its absolute value */
                rowmax = 0.f;
                kx = kc + imax - k;
                i__1 = imax - 1;
                for (j = k;
                        j <= i__1;
                        ++j)
                {
                    if ((r__1 = ap[kx], f2c_abs(r__1)) > rowmax)
                    {
                        rowmax = (r__1 = ap[kx], f2c_abs(r__1));
                        jmax = j;
                    }
                    kx = kx + *n - j;
                    /* L70: */
                }
                kpc = npp - (*n - imax + 1) * (*n - imax + 2) / 2 + 1;
                if (imax < *n)
                {
                    i__1 = *n - imax;
                    jmax = imax + isamax_(&i__1, &ap[kpc + 1], &c__1);
                    /* Computing MAX */
                    r__2 = rowmax;
                    r__3 = (r__1 = ap[kpc + jmax - imax], f2c_abs( r__1)); // , expr subst
                    rowmax = max(r__2,r__3);
                }
                if (absakk >= alpha * colmax * (colmax / rowmax))
                {
                    /* no interchange, use 1-by-1 pivot block */
                    kp = k;
                }
                else if ((r__1 = ap[kpc], f2c_abs(r__1)) >= alpha * rowmax)
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
            kk = k + kstep - 1;
            if (kstep == 2)
            {
                knc = knc + *n - k + 1;
            }
            if (kp != kk)
            {
                /* Interchange rows and columns KK and KP in the trailing */
                /* submatrix A(k:n,k:n) */
                if (kp < *n)
                {
                    i__1 = *n - kp;
                    sswap_(&i__1, &ap[knc + kp - kk + 1], &c__1, &ap[kpc + 1], &c__1);
                }
                kx = knc + kp - kk;
                i__1 = kp - 1;
                for (j = kk + 1;
                        j <= i__1;
                        ++j)
                {
                    kx = kx + *n - j + 1;
                    t = ap[knc + j - kk];
                    ap[knc + j - kk] = ap[kx];
                    ap[kx] = t;
                    /* L80: */
                }
                t = ap[knc];
                ap[knc] = ap[kpc];
                ap[kpc] = t;
                if (kstep == 2)
                {
                    t = ap[kc + 1];
                    ap[kc + 1] = ap[kc + kp - k];
                    ap[kc + kp - k] = t;
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
                    /* A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */
                    r1 = 1.f / ap[kc];
                    i__1 = *n - k;
                    r__1 = -r1;
                    sspr_(uplo, &i__1, &r__1, &ap[kc + 1], &c__1, &ap[kc + *n - k + 1]);
                    /* Store L(k) in column K */
                    i__1 = *n - k;
                    sscal_(&i__1, &r1, &ap[kc + 1], &c__1);
                }
            }
            else
            {
                /* 2-by-2 pivot block D(k): columns K and K+1 now hold */
                /* ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */
                /* where L(k) and L(k+1) are the k-th and (k+1)-th columns */
                /* of L */
                if (k < *n - 1)
                {
                    /* Perform a rank-2 update of A(k+2:n,k+2:n) as */
                    /* A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T */
                    /* = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T */
                    /* where L(k) and L(k+1) are the k-th and (k+1)-th */
                    /* columns of L */
                    d21 = ap[k + 1 + (k - 1) * ((*n << 1) - k) / 2];
                    d11 = ap[k + 1 + k * ((*n << 1) - k - 1) / 2] / d21;
                    d22 = ap[k + (k - 1) * ((*n << 1) - k) / 2] / d21;
                    t = 1.f / (d11 * d22 - 1.f);
                    d21 = t / d21;
                    i__1 = *n;
                    for (j = k + 2;
                            j <= i__1;
                            ++j)
                    {
                        wk = d21 * (d11 * ap[j + (k - 1) * ((*n << 1) - k) / 2] - ap[j + k * ((*n << 1) - k - 1) / 2]);
                        wkp1 = d21 * (d22 * ap[j + k * ((*n << 1) - k - 1) / 2] - ap[j + (k - 1) * ((*n << 1) - k) / 2]);
                        i__2 = *n;
                        for (i__ = j;
                                i__ <= i__2;
                                ++i__)
                        {
                            ap[i__ + (j - 1) * ((*n << 1) - j) / 2] = ap[i__ + (j - 1) * ((*n << 1) - j) / 2] - ap[i__ + (k - 1) * ((*n << 1) - k) / 2] * wk - ap[i__ + k * ((*n << 1) - k - 1) / 2] * wkp1;
                            /* L90: */
                        }
                        ap[j + (k - 1) * ((*n << 1) - k) / 2] = wk;
                        ap[j + k * ((*n << 1) - k - 1) / 2] = wkp1;
                        /* L100: */
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
        kc = knc + *n - k + 2;
        goto L60;
    }
L110:
    return 0;
    /* End of SSPTRF */
}
/* ssptrf_ */
