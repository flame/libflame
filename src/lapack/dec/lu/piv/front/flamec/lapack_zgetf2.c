/*
    Copyright (c) 2021-2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"

/* Subroutine */ integer lapack_zgetf2(integer *m, integer *n, dcomplex *a,
	integer *lda, integer *ipiv, integer *info)
{

/*
    Purpose
    =======

    ZGETF2 computes an LU factorization of a general m-by-n matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 2 BLAS version of the algorithm.

    Arguments
    =========

    M       (input) INTEGER
            The number of rows of the matrix A.  M >= 0.

    N       (input) INTEGER
            The number of columns of the matrix A.  N >= 0.

    A       (input/output) COMPLEX*16 array, dimension (LDA,N)
            On entry, the m by n matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    LDA     (input) INTEGER
            The leading dimension of the array A.  LDA >= fla_max(1,M).

    IPIV    (output) INTEGER array, dimension (fla_min(M,N))
            The pivot indices; for 1 <= i <= fla_min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    INFO    (output) INTEGER
            = 0: successful exit
            < 0: if INFO = -k, the k-th argument had an illegal value
            > 0: if INFO = k, U(k,k) is exactly zero. The factorization
                 has been completed, but the factor U is exactly
                 singular, and division by zero will occur if it is used
                 to solve a system of equations.

    =====================================================================


       Test the input parameters.

       Parameter adjustments */
    /* Table of constant values */
    static TLS_CLASS_SPEC dcomplex c_b1 = {1.,0.};
    static TLS_CLASS_SPEC integer c__1 = 1;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    dcomplex z__1;
    /* Builtin functions */
    void z_div(dcomplex *, dcomplex *, dcomplex *);
    /* Local variables */
    static TLS_CLASS_SPEC integer j;

    static TLS_CLASS_SPEC integer jp;
    extern /* Subroutine */ int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
#define a_subscr(a_1,a_2) (a_2)*a_dim1 + a_1
#define a_ref(a_1,a_2) a[a_subscr(a_1,a_2)]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < fla_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("LAPACK_ZGETF2", &i__1, (ftnlen)13);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = fla_min(*m,*n);
    for (j = 1; j <= i__1; ++j) {

/*        Find pivot and test for singularity. */

	i__2 = *m - j + 1;
	jp = j - 1 + izamax_(&i__2, &a_ref(j, j), &c__1);
	ipiv[j] = jp;
	i__2 = a_subscr(jp, j);
	if (a[i__2].real != 0. || a[i__2].imag != 0.) {

/*           Apply the interchange to columns 1:N. */

	    if (jp != j) {
		zswap_(n, &a_ref(j, 1), lda, &a_ref(jp, 1), lda);
	    }

/*           Compute elements J+1:M of J-th column. */

	    if (j < *m) {
		i__2 = *m - j;
		z_div(&z__1, &c_b1, &a_ref(j, j));
		zscal_(&i__2, &z__1, &a_ref(j + 1, j), &c__1);
	    }

	} else if (*info == 0) {

	    *info = j;
	}

	if (j < fla_min(*m,*n)) {

/*           Update trailing submatrix. */

	    i__2 = *m - j;
	    i__3 = *n - j;
	    z__1.real = -1., z__1.imag = 0.;
	    zgeru_(&i__2, &i__3, &z__1, &a_ref(j + 1, j), &c__1, &a_ref(j, j
		    + 1), lda, &a_ref(j + 1, j + 1), lda);
	}
/* L10: */
    }
    return 0;

/*     End of ZGETF2 */

} /* zgetf2_ */

#undef a_ref
#undef a_subscr

